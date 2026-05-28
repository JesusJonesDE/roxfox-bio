from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
from rich.console import Console

from pipeline.config import Settings


def _load_docking_results(results_dir: Path) -> pd.DataFrame:
    """Collect all docking_results_SCF-*.csv files into one DataFrame."""
    frames = []
    for csv in sorted(results_dir.glob("docking_results_*.csv")):
        df = pd.read_csv(csv)
        if df.empty:
            continue
        scaffold_id = csv.stem.replace("docking_results_", "")
        top = df.sort_values("affinity_kcal_mol").iloc[0]
        frames.append({
            "scaffold_id": scaffold_id,
            "top_affinity_kcal_mol": round(float(top["affinity_kcal_mol"]), 2),
            "n_poses": len(df),
        })
    return pd.DataFrame(frames) if frames else pd.DataFrame(
        columns=["scaffold_id", "top_affinity_kcal_mol", "n_poses"]
    )


def _load_cocrystal_flags(results_dir: Path) -> pd.DataFrame:
    """Extract space group and flagged atom count from cocrystal briefs."""
    rows = []
    for md in sorted(results_dir.glob("cocrystal_brief_*.md")):
        scaffold_id = md.stem.replace("cocrystal_brief_", "")
        text = md.read_text()
        space_group = "?"
        flagged = "?"
        for line in text.splitlines():
            if "| Target co-crystal" in line:
                parts = line.split("|")
                if len(parts) >= 4:
                    space_group = parts[3].strip().replace("Aim to match ", "")
            if "atoms flagged" in line.lower():
                import re
                m = re.search(r"(\d+) atoms? flagged", line)
                if m:
                    flagged = m.group(1)
            if "_No problematic atoms" in line:
                flagged = "0"
        rows.append({
            "scaffold_id": scaffold_id,
            "space_group": space_group,
            "flagged_atoms": flagged,
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["scaffold_id", "space_group", "flagged_atoms"]
    )


def _load_bioactivity(results_dir: Path) -> pd.DataFrame:
    """Load best bioactivity per scaffold from compounds_filtered.csv."""
    csv = results_dir / "compounds_filtered.csv"
    if not csv.exists():
        return pd.DataFrame(columns=["scaffold_id", "best_value_nm", "best_assay_type",
                                     "molecular_weight", "logp", "passes_ro5"])
    df = pd.read_csv(csv)
    if "scaffold_id" not in df.columns:
        return pd.DataFrame()

    # Best compound per scaffold: sort by best_value_nm, keep top row
    df = df.sort_values("best_value_nm", ascending=False)
    df = df.drop_duplicates(subset=["scaffold_id"])
    keep = ["scaffold_id", "best_value_nm", "best_assay_type", "molecular_weight", "logp", "passes_ro5"]
    return df[[c for c in keep if c in df.columns]]


def _write_rank_report(
    gene_symbol: str,
    ranked: pd.DataFrame,
    results_dir: Path,
) -> Path:
    if ranked.empty:
        rows_md = "| — | No docked scaffolds found | — | — | — | — | — |"
    else:
        rows_md = "\n".join(
            f"| **{r.scaffold_id}** "
            f"| {r.get('top_affinity_kcal_mol', '—')} "
            f"| {r.get('best_value_nm', '—')} "
            f"| {r.get('best_assay_type', '—')} "
            f"| {r.get('molecular_weight', '—')} "
            f"| {r.get('logp', '—')} "
            f"| {r.get('flagged_atoms', '—')} |"
            for _, r in ranked.iterrows()
        )

    report = f"""# {gene_symbol} Scaffold Ranking

**Generated**: {datetime.now().strftime("%Y-%m-%d")}
**Target**: {gene_symbol}
**Scaffolds docked**: {len(ranked)}

---

## Ranked by Docking Affinity

*Lower affinity (more negative kcal/mol) = stronger predicted binding.*
*Flagged atoms = atoms > 8 Å from pocket centroid adjacent to rotatable bonds (cocrystal clash risk).*

| Scaffold | Affinity (kcal/mol) | Best activity | Assay type | MW | cLogP | Flagged atoms |
|----------|---------------------|--------------|------------|-----|-------|---------------|
{rows_md}

---

## Interpretation

{"**Top scaffold: " + str(ranked.iloc[0]["scaffold_id"]) + "** — affinity " + str(ranked.iloc[0].get("top_affinity_kcal_mol", "?")) + " kcal/mol." if not ranked.empty else "No docked scaffolds to rank."}

Scaffolds with affinity < −8.0 kcal/mol and 0 flagged atoms are prioritised for:
1. Co-crystallisation experiment (brief already generated if cocrystal was run)
2. Kinase IC50 assay at 1 µM and 10 µM
3. Analog enumeration around the top scaffold core
"""
    path = results_dir / "scaffold_ranking.md"
    path.write_text(report)
    return path


def run_rank(
    gene_symbol: str,
    settings: Settings,
    console: Console,
) -> None:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    dock_df = _load_docking_results(results_dir)
    if dock_df.empty:
        console.print(
            f"  [yellow]{gene_symbol}:[/yellow] no docking results found. "
            f"Run `pipeline dock --target {gene_symbol} --all-scaffolds` first."
        )
        return

    cocrystal_df = _load_cocrystal_flags(results_dir)
    bioactivity_df = _load_bioactivity(results_dir)

    ranked = dock_df.merge(bioactivity_df, on="scaffold_id", how="left")
    if not cocrystal_df.empty:
        ranked = ranked.merge(cocrystal_df, on="scaffold_id", how="left")

    ranked = ranked.sort_values("top_affinity_kcal_mol")

    # Save CSV
    ranked.to_csv(results_dir / "scaffold_ranking.csv", index=False)

    report_path = _write_rank_report(gene_symbol, ranked, results_dir)

    console.print(f"  [dim]{gene_symbol}:[/dim] ranked {len(ranked)} scaffolds")
    for _, r in ranked.head(5).iterrows():
        affinity = r.get("top_affinity_kcal_mol", "?")
        activity = r.get("best_value_nm", "?")
        assay = r.get("best_assay_type", "?")
        console.print(
            f"    [bold]{r['scaffold_id']:12}[/bold] "
            f"affinity {affinity:>6} kcal/mol | "
            f"activity {activity} ({assay})"
        )
    console.print(f"  [dim]{gene_symbol}:[/dim] [green]scaffold_ranking.md written[/green]")
