"""
SAR analysis for a lead scaffold series using RDKit descriptors.

Computes molecular descriptors for all compounds in a series,
correlates with potency, and identifies the key structural features
that drive activity.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from rich.console import Console

from pipeline.config import Settings

if TYPE_CHECKING:
    from pipeline.cache import CacheManager

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

DESCRIPTOR_FUNCS = {
    "mol_wt": Descriptors.MolWt if RDKIT_AVAILABLE else None,
    "logp": Descriptors.MolLogP if RDKIT_AVAILABLE else None,
    "tpsa": Descriptors.TPSA if RDKIT_AVAILABLE else None,
    "hbd": rdMolDescriptors.CalcNumHBD if RDKIT_AVAILABLE else None,
    "hba": rdMolDescriptors.CalcNumHBA if RDKIT_AVAILABLE else None,
    "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds if RDKIT_AVAILABLE else None,
    "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings if RDKIT_AVAILABLE else None,
    "n_rings": rdMolDescriptors.CalcNumRings if RDKIT_AVAILABLE else None,
    "fraction_csp3": rdMolDescriptors.CalcFractionCSP3 if RDKIT_AVAILABLE else None,
    "n_stereocenters": rdMolDescriptors.CalcNumAtomStereoCenters if RDKIT_AVAILABLE else None,
    "qed": None,  # computed separately
}


def _pIC50(nm: float) -> float:
    return -math.log10(max(nm, 0.001) * 1e-9)


def _compute_descriptors(smiles: str) -> dict | None:
    if not RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    result = {}
    for name, func in DESCRIPTOR_FUNCS.items():
        if func is None:
            continue
        try:
            result[name] = round(func(mol), 3)
        except Exception:
            result[name] = None
    try:
        from rdkit.Chem import QED
        result["qed"] = round(QED.qed(mol), 3)
    except Exception:
        result["qed"] = None
    return result


def _pearson(x: pd.Series, y: pd.Series) -> float:
    """Simple Pearson correlation, returns 0 on error."""
    try:
        return round(float(x.corr(y)), 3)
    except Exception:
        return 0.0


def run_sar_analysis(gene: str, series_id: str, settings: Settings, console: Console) -> Path:
    results_dir = settings.results_dir / gene
    compounds_path = results_dir / "compounds_filtered.csv"
    if not compounds_path.exists():
        raise FileNotFoundError(f"{compounds_path} not found")

    df = pd.read_csv(compounds_path)
    series = df[df["scaffold_id"] == series_id].copy()

    if series.empty:
        raise ValueError(f"No compounds found for scaffold {series_id} in {gene}")

    # Add pIC50
    series["pIC50"] = series["best_value_nm"].apply(
        lambda v: _pIC50(v) if v and v > 0 else None
    )
    series = series.dropna(subset=["pIC50", "smiles"])

    console.print(f"  {gene}: SAR analysis on {series_id} ({len(series)} compounds with valid pIC50 + SMILES)")

    if RDKIT_AVAILABLE:
        desc_rows = []
        for _, row in series.iterrows():
            desc = _compute_descriptors(row["smiles"])
            if desc:
                desc["molecule_chembl_id"] = row["molecule_chembl_id"]
                desc_rows.append(desc)
        desc_df = pd.DataFrame(desc_rows)
        series = series.merge(desc_df, on="molecule_chembl_id", how="left", suffixes=("", "_rdkit"))
    else:
        console.print(f"  [yellow]{gene}: RDKit not available — using pre-computed descriptors only[/yellow]")

    # Correlate all numeric columns with pIC50
    # Exclude best_value_nm — it's the source of pIC50, so r=-1 trivially
    EXCLUDE = {"molecule_chembl_id", "scaffold_id", "smiles", "best_assay_type",
               "pIC50", "selectivity_flag", "passes_ro5", "best_value_nm",
               "off_target_flags", "ro5_violations"}
    numeric_cols = [c for c in series.columns
                    if c not in EXCLUDE
                    and pd.api.types.is_numeric_dtype(series[c])]

    correlations = []
    for col in numeric_cols:
        valid = series[[col, "pIC50"]].dropna()
        if len(valid) >= 3:
            r = _pearson(valid[col], valid["pIC50"])
            correlations.append({"descriptor": col, "pearson_r": r, "n": len(valid)})

    corr_df = pd.DataFrame(correlations).sort_values("pearson_r", key=abs, ascending=False)

    # Summary stats per activity tier
    series["potency_tier"] = pd.cut(
        series["best_value_nm"],
        bins=[0, 10, 100, 1000, float("inf")],
        labels=["<10 nM", "10–100 nM", "100–1000 nM", ">1000 nM"],
    )
    tier_stats = series.groupby("potency_tier", observed=True).agg(
        count=("molecule_chembl_id", "count"),
        mean_pIC50=("pIC50", "mean"),
        mean_mw=("molecular_weight", "mean"),
        mean_logp=("logp", "mean"),
        mean_tpsa=("tpsa", "mean") if "tpsa" in series.columns else ("logp", "mean"),
    ).round(2)

    # Save
    sar_csv = results_dir / f"sar_{series_id}.csv"
    series[["molecule_chembl_id", "scaffold_id", "best_value_nm", "pIC50",
            "molecular_weight", "logp", "hbd", "hba", "rotatable_bonds",
            "passes_ro5", "selectivity_flag"] +
           (["tpsa", "aromatic_rings", "fraction_csp3", "qed"] if RDKIT_AVAILABLE else [])
           ].to_csv(sar_csv, index=False)

    _write_sar_report(gene, series_id, series, corr_df, tier_stats, results_dir)

    top_corr = corr_df.iloc[0] if not corr_df.empty else None
    msg = f"  {gene}: SAR {series_id} — strongest driver: {top_corr['descriptor']} (r={top_corr['pearson_r']})" if top_corr is not None else f"  {gene}: SAR {series_id} — no correlations found"
    console.print(msg)
    return sar_csv


def _write_sar_report(gene, series_id, series, corr_df, tier_stats, results_dir):
    lines = [
        f"# {gene} — SAR Analysis: {series_id}",
        "",
        "## Dataset Overview",
        "",
        f"- Compounds analysed: **{len(series)}**",
        f"- Potency range: **{series['best_value_nm'].min():.1f} – {series['best_value_nm'].max():.0f} nM**",
        f"- pIC50 range: **{series['pIC50'].min():.2f} – {series['pIC50'].max():.2f}**",
        "",
        "## Potency Distribution",
        "",
        "| Tier | Count | Mean pIC50 | Mean MW | Mean LogP |",
        "|---|---|---|---|---|",
    ]
    for tier, row in tier_stats.iterrows():
        lines.append(
            f"| {tier} | {int(row['count'])} | {row['mean_pIC50']:.2f} | "
            f"{row['mean_mw']:.0f} | {row['mean_logp']:.2f} |"
        )

    lines += [
        "",
        "## Descriptor–Activity Correlations",
        "",
        "Pearson r between each descriptor and pIC50. Positive r = descriptor increases potency.",
        "",
        "| Descriptor | Pearson r | n | Interpretation |",
        "|---|---|---|---|",
    ]

    interpretations = {
        "logp": "higher lipophilicity → more potent (hydrophobic pocket)",
        "molecular_weight": "larger molecules → more potent (more contact)",
        "tpsa": "lower TPSA → more potent (fewer polar groups in binding)",
        "hbd": "fewer H-bond donors → more potent",
        "hba": "more H-bond acceptors → more interactions",
        "aromatic_rings": "more aromatic rings → π-stacking in pocket",
        "fraction_csp3": "more sp3 character → better 3D fit",
        "qed": "drug-likeness correlates with potency",
        "ro5_violations": "fewer violations → more potent",
        "rotatable_bonds": "more flexibility → better induced fit",
    }

    for _, row in corr_df.head(12).iterrows():
        desc = row["descriptor"]
        r = row["pearson_r"]
        n = int(row["n"])
        direction = "↑ increases" if r > 0 else "↓ decreases"
        interp = interpretations.get(desc, f"{direction} potency")
        lines.append(f"| {desc} | **{r}** | {n} | {interp} |")

    # Top and bottom compounds
    best = series.nsmallest(3, "best_value_nm")
    worst = series.nlargest(3, "best_value_nm")

    lines += [
        "",
        "## Most Potent Compounds",
        "",
        "| ChEMBL ID | Potency (nM) | pIC50 | MW | LogP | Selective |",
        "|---|---|---|---|---|---|",
    ]
    for _, row in best.iterrows():
        sel = "Yes" if not row["selectivity_flag"] else "No"
        lines.append(
            f"| {row['molecule_chembl_id']} | {row['best_value_nm']:.1f} | {row['pIC50']:.2f} | "
            f"{row['molecular_weight']:.0f} | {row['logp']:.2f} | {sel} |"
        )

    lines += [
        "",
        "## Least Potent Compounds",
        "",
        "| ChEMBL ID | Potency (nM) | pIC50 | MW | LogP | Selective |",
        "|---|---|---|---|---|---|",
    ]
    for _, row in worst.iterrows():
        sel = "Yes" if not row["selectivity_flag"] else "No"
        lines.append(
            f"| {row['molecule_chembl_id']} | {row['best_value_nm']:.0f} | {row['pIC50']:.2f} | "
            f"{row['molecular_weight']:.0f} | {row['logp']:.2f} | {sel} |"
        )

    lines += [
        "",
        "## SAR Summary",
        "",
    ]

    # Narrative from correlations
    strong = corr_df[corr_df["pearson_r"].abs() >= 0.3]
    if strong.empty:
        lines.append("No strong (|r| ≥ 0.3) single-descriptor correlations with potency found.")
        lines.append("This may indicate multi-parameter optimisation is required, or that the series")
        lines.append("has limited structural variation.")
    else:
        lines.append(f"**{len(strong)} descriptors show meaningful correlation (|r| ≥ 0.3) with potency:**")
        lines.append("")
        for _, row in strong.iterrows():
            direction = "increases" if row["pearson_r"] > 0 else "decreases"
            lines.append(f"- **{row['descriptor']}** (r={row['pearson_r']}): {direction} potency — "
                         f"{interpretations.get(row['descriptor'], '')}")

    out_md = results_dir / f"sar_report_{series_id}.md"
    out_md.write_text("\n".join(lines) + "\n")
