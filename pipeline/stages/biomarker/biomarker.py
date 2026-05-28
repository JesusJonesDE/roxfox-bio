from __future__ import annotations

import io
from datetime import datetime
from pathlib import Path

import httpx
import pandas as pd
from scipy.stats import fisher_exact
from rich.console import Console
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from pipeline.cache import CacheManager
from pipeline.config import Settings

_STRONGLY_DEP = -0.5
_MIN_MUT_LINES = 3  # minimum lines with mutation for inclusion


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _http_get(url: str, timeout: int = 120) -> httpx.Response:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(url)
        r.raise_for_status()
        return r


def _get_manifest() -> dict[str, str]:
    r = _http_get("https://depmap.org/portal/api/download/files")
    try:
        data = r.json()
        files = data if isinstance(data, list) else data.get("files", data.get("table", []))
        result: dict[str, str] = {}
        for entry in (files if isinstance(files, list) else []):
            name = entry.get("fileName") or entry.get("name") or entry.get("file_name")
            url = entry.get("downloadUrl") or entry.get("download_url") or entry.get("url")
            if name and url:
                result[str(name)] = str(url)
        if result:
            return result
    except Exception:
        pass
    # CSV fallback (manifest is a CSV: release, release_date, filename, url, md5_hash)
    df = pd.read_csv(io.StringIO(r.text))
    col_lower = {c.lower().replace(" ", "_"): c for c in df.columns}
    name_col = next((col_lower[k] for k in col_lower if "file_name" in k or "filename" in k),
                    next((col_lower[k] for k in col_lower if "name" in k), None))
    url_col = next((col_lower[k] for k in col_lower if "url" in k or "link" in k), None)
    if name_col and url_col:
        # Use setdefault so first occurrence (latest release) wins for duplicate filenames
        result: dict[str, str] = {}
        for name, url in zip(df[name_col], df[url_col]):
            if pd.notna(name) and pd.notna(url):
                result.setdefault(str(name), str(url))
        return result
    raise ValueError("Cannot parse DepMap manifest")


def _download_mutations(manifest: dict[str, str], cache_dir: Path) -> pd.DataFrame:
    """Download OmicsSomaticMutations.csv, filter to non-silent deleterious mutations,
    deduplicate to one row per (ModelID, Hugo_Symbol)."""
    cache_file = cache_dir / "OmicsSomaticMutations_filtered.parquet"
    if cache_file.exists():
        return pd.read_parquet(cache_file)

    url = manifest.get("OmicsSomaticMutations.csv") or next(
        (v for k, v in manifest.items()
         if "OmicsSomaticMutations" in k and k.endswith(".csv")),
        None,
    )
    if not url:
        available = ", ".join(list(manifest.keys())[:10])
        raise ValueError(
            f"OmicsSomaticMutations.csv not found in DepMap manifest. "
            f"Available (first 10): {available}"
        )

    # Peek at headers to handle old (24Q4-) vs new (26Q1+) schema
    header_df = pd.read_csv(url, nrows=0)
    cols = set(header_df.columns)

    # New schema (26Q1+): HugoSymbol, MolecularConsequence, LikelyLoF
    # Old schema (-24Q4): Hugo_Symbol, Variant_Classification, isDeleterious
    if "HugoSymbol" in cols:
        df = pd.read_csv(url, usecols=["ModelID", "HugoSymbol", "MolecularConsequence",
                                        "LikelyLoF"])
        df = df.rename(columns={
            "HugoSymbol": "Hugo_Symbol",
            "MolecularConsequence": "Variant_Classification",
            "LikelyLoF": "isDeleterious",
        })
        # MolecularConsequence "synonymous_variant" ≈ old "Silent"
        df = df[
            (df["Variant_Classification"] != "synonymous_variant") &
            (df["isDeleterious"] == True)  # noqa: E712
        ].copy()
    else:
        df = pd.read_csv(url, usecols=["ModelID", "Hugo_Symbol", "Variant_Classification",
                                        "isDeleterious"])
        df = df[
            (df["Variant_Classification"] != "Silent") &
            (df["isDeleterious"] == True)  # noqa: E712
        ].copy()

    # One row per (ModelID, Hugo_Symbol) — any deleterious mutation in gene counts
    df = df.drop_duplicates(subset=["ModelID", "Hugo_Symbol"])

    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache_file, index=False)
    return df


def _compute_enrichment(
    dep_records: list[dict],
    mutations: pd.DataFrame,
    lineage: str,
    min_lines: int,
) -> pd.DataFrame:
    """Compute Fisher's exact enrichment per gene for a lineage.

    Cell lines in dep_records are split into:
      - strongly dependent: gene_effect <= _STRONGLY_DEP AND OncotreeLineage == lineage
      - non-dependent: OncotreeLineage == lineage AND gene_effect > _STRONGLY_DEP
    """
    dep_df = pd.DataFrame(dep_records)
    lin = dep_df[dep_df["OncotreeLineage"] == lineage]
    if len(lin) == 0:
        available = sorted(dep_df["OncotreeLineage"].unique())
        raise ValueError(
            f"Lineage '{lineage}' not found in DepMap cache. "
            f"Available lineages: {available}"
        )

    dependent_ids = set(lin.loc[lin["gene_effect"] <= _STRONGLY_DEP, "ModelID"])
    nondependent_ids = set(lin.loc[lin["gene_effect"] > _STRONGLY_DEP, "ModelID"])
    all_ids = dependent_ids | nondependent_ids

    # Restrict mutations to this lineage's cell lines
    muts = mutations[mutations["ModelID"].isin(all_ids)]
    genes = muts["Hugo_Symbol"].unique()

    rows = []
    for gene in genes:
        gene_ids = set(muts.loc[muts["Hugo_Symbol"] == gene, "ModelID"])
        n_dep_mut = len(gene_ids & dependent_ids)
        n_dep_nomut = len(dependent_ids) - n_dep_mut
        n_nondep_mut = len(gene_ids & nondependent_ids)
        n_nondep_nomut = len(nondependent_ids) - n_nondep_mut

        if n_dep_mut < min_lines:
            continue

        table = [[n_dep_mut, n_dep_nomut], [n_nondep_mut, n_nondep_nomut]]
        odds_ratio, p_value = fisher_exact(table, alternative="greater")

        if odds_ratio == float("inf"):
            direction = "enriched_in_dependent"
        elif p_value < 0.05 and odds_ratio > 1:
            direction = "enriched_in_dependent"
        elif p_value < 0.05 and odds_ratio < 1:
            direction = "depleted_in_dependent"
        else:
            direction = "ns"

        rows.append({
            "gene": gene,
            "lineage": lineage,
            "variant_classification": "any_deleterious",
            "n_dependent_with_mut": n_dep_mut,
            "n_dependent_without_mut": n_dep_nomut,
            "n_nondependent_with_mut": n_nondep_mut,
            "n_nondependent_without_mut": n_nondep_nomut,
            "odds_ratio": round(odds_ratio, 3),
            "p_value": round(p_value, 6),
            "enrichment_direction": direction,
        })

    if not rows:
        return pd.DataFrame(columns=[
            "gene", "lineage", "variant_classification",
            "n_dependent_with_mut", "n_dependent_without_mut",
            "n_nondependent_with_mut", "n_nondependent_without_mut",
            "odds_ratio", "p_value", "enrichment_direction",
        ])

    return pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)


def _write_report(
    gene_symbol: str,
    lineage: str,
    results: pd.DataFrame,
    n_dependent: int,
    n_nondependent: int,
    results_dir: Path,
) -> Path:
    significant = results[results["p_value"] < 0.05]
    enriched = significant[significant["enrichment_direction"] == "enriched_in_dependent"]

    if len(enriched) == 0:
        top_section = "_No co-mutations significantly enriched in strongly-dependent lines (p < 0.05)._"
    else:
        top_rows = [
            f"- **{r.gene}**: OR={r.odds_ratio:.2f}, p={r.p_value:.4f}, "
            f"{r.n_dependent_with_mut} dependent lines vs {r.n_nondependent_with_mut} non-dependent"
            for _, r in enriched.head(5).iterrows()
        ]
        top_section = "\n".join(top_rows)

    table_rows = [
        f"| {r.gene} | {r.odds_ratio:.2f} | {r.p_value:.4f} "
        f"| {r.n_dependent_with_mut} | {r.n_nondependent_with_mut} "
        f"| {r.enrichment_direction.replace('_', ' ').title()} |"
        for _, r in results.head(30).iterrows()
    ]

    sig_count = len(significant)
    total_tested = len(results)

    report = f"""# {gene_symbol} Biomarker Analysis — {lineage}

**Source**: Broad Institute DepMap + CCLE OmicsSomaticMutations
**Generated**: {datetime.now().strftime("%Y-%m-%d")}
**Lineage**: {lineage}
**Cell lines**: {n_dependent} strongly dependent | {n_nondependent} non-dependent
**Genes tested**: {total_tested:,} | **Significant (p < 0.05)**: {sig_count}

---

## Top Candidate Biomarkers

{top_section}

---

## Full Ranked Table (top 30)

*Sorted by p-value ascending. Fisher's exact test, alternative="greater" (enrichment in dependent lines).*

| Gene | Odds Ratio | p-value | Dep. lines w/ mut | Non-dep. lines w/ mut | Direction |
|------|-----------|---------|------------------|----------------------|-----------|
{chr(10).join(table_rows) if table_rows else "| — | No genes passed the minimum cell-line filter | — | — | — | — |"}

---

## Interpretation

{"A co-mutation biomarker was identified. Patients with mutations in the top candidate gene(s) may show higher VRK1 dependency in this lineage, supporting use as a patient stratification marker in clinical trial design." if len(enriched) > 0 else "No significant co-mutation biomarker was identified in this lineage at p < 0.05. VRK1 dependency in " + lineage + " may be driven by non-mutational mechanisms (epigenetic, expression-level) or the lineage may be too small for statistical power."}

---

## Notes

- Strongly dependent: Chronos score ≤ −0.5
- Mutations: non-silent, isDeleterious == True; deduplicated to one entry per cell line per gene
- Minimum {_MIN_MUT_LINES} dependent lines with mutation required for inclusion
"""
    path = results_dir / "biomarker_report.md"
    path.write_text(report)
    return path


def update_research_report(
    gene_symbol: str,
    lineage: str,
    results: pd.DataFrame,
    results_dir: Path,
) -> bool:
    report_path = results_dir.parent / "research_report.md"
    if not report_path.exists():
        return False
    text = report_path.read_text()
    if "Biomarker Analysis" in text:
        return False

    sig = results[results["p_value"] < 0.05]
    top = sig.head(3)
    rows = [
        f"| {r.gene} | {r.odds_ratio:.2f} | {r.p_value:.4f} | {r.enrichment_direction} |"
        for _, r in top.iterrows()
    ]
    section = (
        f"\n\n### Biomarker Analysis ({lineage})\n\n"
        f"Significant co-mutations: {len(sig)} | Genes tested: {len(results)}\n\n"
        f"| Gene | OR | p-value | Direction |\n"
        f"|------|-----|---------|----------|\n"
        + "\n".join(rows)
        + f"\n\n*Full report: [{gene_symbol}/biomarker_report.md]({gene_symbol}/biomarker_report.md)*\n"
    )

    marker = f"## {gene_symbol}"
    if marker in text:
        idx = text.index(marker)
        end_of_line = text.find("\n", idx)
        text = text[:end_of_line + 1] + section + text[end_of_line + 1:]
        report_path.write_text(text)
        return True

    report_path.write_text(text + section)
    return True


def run_biomarker(
    gene_symbol: str,
    lineage: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    min_lines: int,
    console: Console,
) -> None:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    cache_key = f"biomarker_{lineage.replace(' ', '_').replace('/', '_')}"

    if not force:
        cached = cache.load(gene_symbol, cache_key)
        if cached is not None:
            count = len(cached) if isinstance(cached, list) else "?"
            console.print(
                f"  [dim]{gene_symbol:10}[/dim] biomarker            "
                f"[yellow]SKIP[/yellow]  (cached {count} results)"
            )
            return

    # Load DepMap dependency cache
    dep_cached = cache.load(gene_symbol, "depmap")
    if dep_cached is None:
        raise ValueError(
            f"DepMap cache not found for {gene_symbol}. "
            f"Run `pipeline depmap --target {gene_symbol}` first."
        )

    console.print(f"  [dim]{gene_symbol}:[/dim] loading dependency cache ({len(dep_cached):,} lines)...")
    console.print(f"  [dim]{gene_symbol}:[/dim] fetching DepMap manifest...")
    manifest = _get_manifest()

    console.print(f"  [dim]{gene_symbol}:[/dim] downloading OmicsSomaticMutations.csv (~150 MB)...")
    mut_cache_dir = settings.cache_dir / gene_symbol
    mutations = _download_mutations(manifest, mut_cache_dir)

    console.print(f"  [dim]{gene_symbol}:[/dim] analysing {lineage}...")
    results = _compute_enrichment(dep_cached, mutations, lineage, min_lines)

    dep_df = pd.DataFrame(dep_cached)
    lin = dep_df[dep_df["OncotreeLineage"] == lineage]
    n_dependent = int((lin["gene_effect"] <= _STRONGLY_DEP).sum())
    n_nondependent = int((lin["gene_effect"] > _STRONGLY_DEP).sum())

    sig_count = int((results["p_value"] < 0.05).sum())
    top_gene = results.iloc[0]["gene"] if len(results) > 0 else "none"
    top_or = results.iloc[0]["odds_ratio"] if len(results) > 0 else float("nan")
    top_p = results.iloc[0]["p_value"] if len(results) > 0 else float("nan")

    # Save cache
    cache.save(gene_symbol, cache_key, results.to_dict("records"), len(results))

    # Write CSV and report
    results.to_csv(results_dir / "biomarker_results.csv", index=False)
    _write_report(gene_symbol, lineage, results, n_dependent, n_nondependent, results_dir)
    updated = update_research_report(gene_symbol, lineage, results, results_dir)

    console.print(
        f"  [dim]{gene_symbol}:[/dim] {len(results):,} genes tested | "
        f"{sig_count} significant (p < 0.05) | "
        f"top: {top_gene} (OR={top_or:.2f}, p={top_p:.4f})"
    )
    console.print(f"  [dim]{gene_symbol}:[/dim] [green]biomarker_report.md written[/green]")
    if updated:
        console.print(f"  [dim]{gene_symbol}:[/dim] [dim]research_report.md updated[/dim]")
