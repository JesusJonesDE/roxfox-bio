"""VRK1 scaffold triage: scores and ranks scaffolds from clean compounds."""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
from rich.console import Console

from pipeline.config import Settings


def run_scaffold_triage(gene: str, settings: Settings, console: Console) -> Path:
    results_dir = settings.results_dir / gene
    compounds_path = results_dir / "compounds_filtered.csv"
    scaffolds_path = results_dir / "scaffolds.csv"

    if not compounds_path.exists():
        raise FileNotFoundError(f"{compounds_path} not found — run analyze first")
    if not scaffolds_path.exists():
        raise FileNotFoundError(f"{scaffolds_path} not found — run analyze first")

    compounds = pd.read_csv(compounds_path)
    scaffolds = pd.read_csv(scaffolds_path)

    # Clean compounds: Ro5-passing, selective, positive potency only
    clean = compounds[
        (compounds["passes_ro5"] == True) &
        (compounds["selectivity_flag"] == False) &
        (compounds["best_value_nm"] > 0)
    ].copy()

    if clean.empty:
        console.print(f"  [yellow]{gene}: no clean compounds after filtering[/yellow]")
        return results_dir / "scaffold_triage.csv"

    console.print(f"  {gene}: {len(clean)} clean compounds across {clean['scaffold_id'].nunique()} scaffolds")

    # Per-scaffold metrics computed from clean compounds only
    def _pIC50(nm: float) -> float:
        return -math.log10(max(nm, 0.001) * 1e-9)

    clean["pIC50"] = clean["best_value_nm"].apply(_pIC50)

    scaffold_stats = (
        clean.groupby("scaffold_id")
        .agg(
            clean_count=("molecule_chembl_id", "count"),
            best_potency_nm=("best_value_nm", "min"),
            median_potency_nm=("best_value_nm", "median"),
            best_pIC50=("pIC50", "max"),
        )
        .reset_index()
    )

    # Total compound count per scaffold (including non-clean)
    total_per_scaffold = compounds.groupby("scaffold_id")["molecule_chembl_id"].count().rename("total_count")
    scaffold_stats = scaffold_stats.join(total_per_scaffold, on="scaffold_id")
    scaffold_stats["selectivity_rate"] = scaffold_stats["clean_count"] / scaffold_stats["total_count"]

    # Merge scaffold SMILES
    scaffold_stats = scaffold_stats.merge(
        scaffolds[["scaffold_id", "scaffold_smiles"]], on="scaffold_id", how="left"
    )

    # Normalise sub-scores to [0, 1] within this set
    def _minmax(col: pd.Series) -> pd.Series:
        lo, hi = col.min(), col.max()
        return (col - lo) / (hi - lo) if hi > lo else pd.Series([1.0] * len(col), index=col.index)

    scaffold_stats["score_potency"] = _minmax(scaffold_stats["best_pIC50"])
    scaffold_stats["score_selectivity"] = _minmax(scaffold_stats["selectivity_rate"])
    scaffold_stats["score_size"] = _minmax(scaffold_stats["clean_count"].apply(math.log1p))

    scaffold_stats["composite_score"] = (
        0.40 * scaffold_stats["score_potency"] +
        0.40 * scaffold_stats["score_selectivity"] +
        0.20 * scaffold_stats["score_size"]
    ).round(4)

    scaffold_stats = scaffold_stats.sort_values("composite_score", ascending=False).reset_index(drop=True)
    scaffold_stats.insert(0, "rank", scaffold_stats.index + 1)

    # Series rank: scaffolds with ≥2 clean compounds (actionable SAR)
    series = scaffold_stats[scaffold_stats["clean_count"] >= 2].copy().reset_index(drop=True)
    series.insert(0, "series_rank", series.index + 1)
    scaffold_stats = scaffold_stats.merge(series[["scaffold_id", "series_rank"]], on="scaffold_id", how="left")
    scaffold_stats["is_series"] = scaffold_stats["series_rank"].notna()

    out_cols = [
        "rank", "series_rank", "scaffold_id", "scaffold_smiles", "composite_score",
        "score_potency", "score_selectivity", "score_size",
        "best_potency_nm", "median_potency_nm", "best_pIC50",
        "clean_count", "total_count", "selectivity_rate", "is_series",
    ]
    triage_path = results_dir / "scaffold_triage.csv"
    scaffold_stats[out_cols].to_csv(triage_path, index=False)

    # Leads from top-5 series (multi-compound scaffolds only)
    top_series_ids = series.head(5)["scaffold_id"].tolist()
    leads = clean[clean["scaffold_id"].isin(top_series_ids)].copy()
    leads = leads.merge(
        series[["scaffold_id", "series_rank", "composite_score"]], on="scaffold_id", how="left"
    )
    leads = leads.sort_values(["series_rank", "best_value_nm"]).reset_index(drop=True)

    leads_path = results_dir / "triage_leads.csv"
    leads.to_csv(leads_path, index=False)

    n_series = len(series)
    n_singletons = len(scaffold_stats) - n_series
    console.print(f"  {gene}: {n_series} series scaffolds (≥2 compounds), {n_singletons} singletons")
    console.print(f"  scaffold_triage.csv ({len(scaffold_stats)} scaffolds), triage_leads.csv ({len(leads)} leads)")
    if not series.empty:
        top = series.iloc[0]
        console.print(f"  Top series: {top['scaffold_id']}  "
                      f"score={top['composite_score']}  "
                      f"best={top['best_potency_nm']:.1f} nM  "
                      f"clean={int(top['clean_count'])}/{int(top['total_count'])}")
    return triage_path
