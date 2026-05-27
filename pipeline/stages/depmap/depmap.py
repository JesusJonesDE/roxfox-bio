from __future__ import annotations

import io
from datetime import datetime
from pathlib import Path

import httpx
import pandas as pd
from rich.console import Console
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from pipeline.cache import CacheManager
from pipeline.config import Settings

MANIFEST_URL = "https://depmap.org/portal/api/download/files"
_STRONGLY_DEP = -0.5
_MODERATELY_DEP = -0.3
_WEAKLY_DEP = -0.1
_PAN_ESSENTIAL_FRACTION = 0.70
_MIN_LINEAGE_LINES = 3


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_manifest() -> httpx.Response:
    with httpx.Client(timeout=60) as client:
        r = client.get(MANIFEST_URL, headers={"Accept": "application/json, text/csv"})
        r.raise_for_status()
        return r


def _parse_manifest(r: httpx.Response) -> dict[str, str]:
    """Return {filename: signed_url} from DepMap manifest (JSON or CSV format)."""
    # Try JSON (newer API)
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

    # CSV fallback
    df = pd.read_csv(io.StringIO(r.text))
    col_lower = {c.lower().replace(" ", "_"): c for c in df.columns}
    name_col = next(
        (col_lower[k] for k in col_lower if "file_name" in k or "filename" in k),
        next((col_lower[k] for k in col_lower if "name" in k), None),
    )
    url_col = next((col_lower[k] for k in col_lower if "url" in k or "link" in k), None)
    if name_col and url_col:
        return {str(k): str(v) for k, v in zip(df[name_col], df[url_col]) if pd.notna(v)}

    raise ValueError(
        f"Cannot parse DepMap manifest. Columns found: {list(df.columns if 'df' in dir() else [])}"
    )


def _classify_tier(effect: float) -> str:
    if effect <= _STRONGLY_DEP:
        return "strongly_dependent"
    if effect <= _MODERATELY_DEP:
        return "moderately_dependent"
    if effect <= _WEAKLY_DEP:
        return "weakly_dependent"
    return "not_essential"


def _download_gene_effects(url: str, gene_symbol: str) -> pd.DataFrame:
    """Stream CRISPRGeneEffect.csv keeping only the target gene column + ModelID.

    DepMap's first column is unnamed in the CSV header (empty string); ModelID
    values (ACH-XXXXXX) are the row labels.  We use a two-pass approach:
    first read headers only to find the gene column's integer position, then
    read the full file using integer usecols so the unnamed column 0 is captured.
    """
    prefix = f"{gene_symbol} ("

    # Pass 1: header only to find gene column integer position
    header = pd.read_csv(url, nrows=0)
    col_list = list(header.columns)
    gene_col_candidates = [(i, c) for i, c in enumerate(col_list) if c.startswith(prefix)]
    if not gene_col_candidates:
        raise ValueError(
            f"{gene_symbol} column not found in CRISPRGeneEffect.csv. "
            f"Expected a column like '{gene_symbol} (7443)'. "
            "Check that the gene symbol matches DepMap annotation."
        )
    gene_col_idx, gene_col_name = gene_col_candidates[0]

    # Pass 2: full data using integer positions (col 0 = unnamed ModelID)
    df = pd.read_csv(url, usecols=[0, gene_col_idx], header=0)
    df.columns = ["ModelID", "gene_effect"]
    df["gene_effect"] = pd.to_numeric(df["gene_effect"], errors="coerce")
    return df.dropna(subset=["gene_effect"]).reset_index(drop=True)


def _download_model(url: str) -> pd.DataFrame:
    """Stream Model.csv keeping only lineage columns."""
    df = pd.read_csv(
        url,
        usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease",
                 "OncotreeSubtype", "CCLEName"],
    )
    df["OncotreeLineage"] = df["OncotreeLineage"].fillna("Unknown").astype(str)
    return df


def _aggregate_lineages(merged: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    """Compute per-lineage LineageSummary records and global pan_essential_flag."""
    merged = merged.copy()
    merged["is_strongly_dep"] = merged["gene_effect"] <= _STRONGLY_DEP

    total = len(merged)
    pan_essential = bool(
        (merged["is_strongly_dep"].sum() / total) > _PAN_ESSENTIAL_FRACTION
        if total > 0 else False
    )

    agg = (
        merged.groupby("OncotreeLineage")
        .agg(
            n_lines=("gene_effect", "count"),
            median_effect=("gene_effect", "median"),
            mean_effect=("gene_effect", "mean"),
            n_strongly_dependent=("is_strongly_dep", "sum"),
        )
        .reset_index()
    )
    agg["pct_strongly_dependent"] = (
        (agg["n_strongly_dependent"] / agg["n_lines"]) * 100
    ).round(1)
    agg["dependency_tier"] = agg["median_effect"].apply(_classify_tier)
    agg["pan_essential_flag"] = pan_essential

    ranked = (
        agg[agg["n_lines"] >= _MIN_LINEAGE_LINES]
        .sort_values("median_effect")
        .reset_index(drop=True)
    )
    return ranked, pan_essential


def _write_report(
    gene_symbol: str,
    ranked: pd.DataFrame,
    pan_essential: bool,
    total_lines: int,
    results_dir: Path,
) -> Path:
    rows = []
    for _, row in ranked.iterrows():
        tier_label = row.dependency_tier.replace("_", " ").title()
        rows.append(
            f"| {row.OncotreeLineage} | {int(row.n_lines)} | {row.median_effect:+.3f} "
            f"| {int(row.pct_strongly_dependent)}% | {tier_label} |"
        )

    top3 = ranked.head(3)
    interp = [
        f"- **{r.OncotreeLineage}**: {int(r.n_lines)} lines, "
        f"median effect {r.median_effect:+.2f}, "
        f"{int(r.pct_strongly_dependent)}% strongly dependent"
        for _, r in top3.iterrows()
    ]

    pan_text = (
        "**Yes** — VRK1 is essential in >70% of all screened lines. "
        "Therapeutic window may be narrow; lineage-selective profiling recommended."
        if pan_essential else
        "**No** — VRK1 dependency is lineage-selective, supporting a targeted indication strategy."
    )

    report = f"""# {gene_symbol} DepMap Cancer Dependency Analysis

**Source**: Broad Institute DepMap Public (CRISPR Chronos method)
**Generated**: {datetime.now().strftime("%Y-%m-%d")}
**Cell lines screened**: {total_lines:,}
**Lineages with ≥{_MIN_LINEAGE_LINES} lines**: {len(ranked)}

---

## Pan-Essential Assessment

Pan-essential (>70% of all lines, gene effect ≤ −0.5): {pan_text}

---

## Cancer Lineage Dependency Ranking

*Sorted by median gene effect (most negative = most dependent). Lineages with <{_MIN_LINEAGE_LINES} cell lines excluded.*

| Lineage | Cell Lines | Median Effect | % Strongly Dep. | Dependency Tier |
|---------|-----------|--------------|----------------|-----------------|
{chr(10).join(rows)}

---

## Top Indication Candidates

{chr(10).join(interp) if interp else "No lineages meet the strongly dependent threshold."}

---

## Notes

- Chronos score ≤ −0.5: strongly dependent (essential for cell survival)
- Chronos score −0.5 to −0.3: moderately dependent
- Chronos score −0.3 to −0.1: weakly dependent
- Chronos score > −0.1: not essential
"""
    path = results_dir / "depmap_report.md"
    path.write_text(report)
    return path


def update_research_report(
    gene_symbol: str,
    ranked: pd.DataFrame,
    pan_essential: bool,
    results_dir: Path,
) -> bool:
    """Inject DepMap findings into data/results/research_report.md. Returns True if updated."""
    report_path = results_dir.parent / "research_report.md"
    if not report_path.exists():
        return False

    text = report_path.read_text()
    if "DepMap Cancer Dependency" in text:
        return False  # already injected

    top = ranked.head(3)
    rows = [
        f"| {r.OncotreeLineage} | {r.median_effect:+.2f} "
        f"| {int(r.pct_strongly_dependent)}% | {r.dependency_tier.replace('_', ' ').title()} |"
        for _, r in top.iterrows()
    ]

    section = (
        f"\n\n### DepMap Cancer Dependency (CRISPR Chronos)\n\n"
        f"Pan-essential: {'Yes' if pan_essential else 'No'} | "
        f"Lineages with ≥3 lines: {len(ranked)}\n\n"
        f"| Lineage | Median Effect | % Strongly Dep. | Tier |\n"
        f"|---------|--------------|----------------|------|\n"
        + "\n".join(rows)
        + f"\n\n*Full report: [{gene_symbol}/depmap_report.md]({gene_symbol}/depmap_report.md)*\n"
    )

    # Inject after the gene's primary section header
    marker = f"## {gene_symbol}"
    if marker in text:
        idx = text.index(marker)
        end_of_line = text.find("\n", idx)
        text = text[: end_of_line + 1] + section + text[end_of_line + 1 :]
        report_path.write_text(text)
        return True

    # Fallback: append to end
    report_path.write_text(text + section)
    return True


def run_depmap(
    gene_symbol: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> None:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    # Cache hit
    if not force:
        cached = cache.load(gene_symbol, "depmap")
        if cached is not None:
            count = len(cached) if isinstance(cached, list) else "?"
            console.print(
                f"  [dim]{gene_symbol:10}[/dim] depmap               "
                f"[yellow]SKIP[/yellow]  (cached {count} records)"
            )
            return

    console.print(f"  [dim]{gene_symbol}:[/dim] fetching DepMap manifest...")
    manifest = _parse_manifest(_get_manifest())

    # Exact match first (avoids CRISPRGeneEffectUncorrected.csv and older releases)
    crispr_url = next(
        (v for k, v in manifest.items() if k == "CRISPRGeneEffect.csv"),
        next(
            (v for k, v in manifest.items()
             if k.startswith("CRISPRGeneEffect") and "Uncorrected" not in k),
            None,
        ),
    )
    model_url = next(
        (v for k, v in manifest.items() if k == "Model.csv"),
        next(
            (v for k, v in manifest.items()
             if k.endswith(".csv") and "Model" in k and "Gene" not in k
             and "Mutation" not in k and "Copy" not in k),
            None,
        ),
    )

    if not crispr_url:
        available = ", ".join(list(manifest.keys())[:10])
        raise ValueError(
            f"CRISPRGeneEffect.csv not found in DepMap manifest. "
            f"Available files (first 10): {available}"
        )
    if not model_url:
        raise ValueError("Model.csv not found in DepMap manifest.")

    console.print(
        f"  [dim]{gene_symbol}:[/dim] downloading gene effect data (~200 MB, 2–4 min)..."
    )
    gene_effects = _download_gene_effects(crispr_url, gene_symbol)

    console.print(f"  [dim]{gene_symbol}:[/dim] downloading cell line metadata...")
    model_df = _download_model(model_url)

    merged = gene_effects.merge(model_df, on="ModelID", how="left")
    merged["OncotreeLineage"] = merged["OncotreeLineage"].fillna("Unknown")
    total_lines = len(merged)

    ranked, pan_essential = _aggregate_lineages(merged)

    # Cache raw records
    records = (
        merged[
            ["ModelID", "gene_effect", "OncotreeLineage",
             "OncotreePrimaryDisease", "OncotreeSubtype", "CCLEName"]
        ]
        .to_dict("records")
    )
    cache.save(gene_symbol, "depmap", records, len(records))

    # Write CSV and report
    ranked.to_csv(results_dir / "depmap_lineage_summary.csv", index=False)
    _write_report(gene_symbol, ranked, pan_essential, total_lines, results_dir)

    # Inject into master research report (US4)
    updated = update_research_report(gene_symbol, ranked, pan_essential, results_dir)

    # Console summary
    top = ranked.iloc[0] if len(ranked) > 0 else None
    top_str = (
        f"{top.OncotreeLineage} (median {top.median_effect:+.2f}, "
        f"{int(top.pct_strongly_dependent)}% strongly dependent)"
        if top is not None else "none"
    )
    n_lineages = merged["OncotreeLineage"].nunique()

    console.print(
        f"  [dim]{gene_symbol}:[/dim] {total_lines:,} cell lines across {n_lineages} lineages"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] pan-essential: {'Yes' if pan_essential else 'No'}"
    )
    console.print(f"  [dim]{gene_symbol}:[/dim] top lineage: {top_str}")
    console.print(f"  [dim]{gene_symbol}:[/dim] [green]depmap_report.md written[/green]")
    if updated:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [dim]research_report.md updated[/dim]"
        )
