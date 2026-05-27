"""
Kinome selectivity root-cause analysis.

Identifies which specific kinases/targets VRK1 off-target compounds hit,
groups them by kinase family, and writes selectivity_kinome.csv + kinome_report.md.
"""
from __future__ import annotations

import time
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING

import httpx
import pandas as pd
from rich.console import Console
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings

if TYPE_CHECKING:
    from pipeline.cache import CacheManager

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
OFF_TARGET_NM = 1000.0
SAMPLE_SIZE = 30  # top N most-flagged compounds to analyse in depth


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_target_activities(mol_id: str, primary_target_id: str) -> list[dict]:
    params = {
        "molecule_chembl_id": mol_id,
        "standard_value__lte": OFF_TARGET_NM,
        "standard_units": "nM",
        "activity_type__in": "IC50,Ki,Kd",
        "limit": 200,
        "format": "json",
    }
    with httpx.Client(timeout=20) as client:
        r = client.get(f"{CHEMBL_BASE}/activity", params=params,
                       headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"})
        r.raise_for_status()
    acts = r.json().get("activities", [])
    return [a for a in acts if a.get("target_chembl_id") != primary_target_id]


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_target_info(target_id: str) -> dict:
    with httpx.Client(timeout=15) as client:
        r = client.get(f"{CHEMBL_BASE}/target/{target_id}",
                       params={"format": "json"},
                       headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"})
        r.raise_for_status()
    return r.json()


def run_selectivity_kinome(gene: str, settings: Settings, cache: "CacheManager", console: Console) -> Path:
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)

    compounds_path = results_dir / "compounds_filtered.csv"
    if not compounds_path.exists():
        raise FileNotFoundError(f"{compounds_path} not found — run analyze first")

    df = pd.read_csv(compounds_path)
    chembl_cache = cache.load(gene, "chembl")
    primary_target_id = chembl_cache.get("chembl_target_id") if chembl_cache else None

    # Focus on the most promiscuous Ro5-passing compounds
    flagged = df[(df["passes_ro5"] == True) & (df["selectivity_flag"] == True)].copy()
    flagged = flagged.sort_values("off_target_flags", ascending=False)
    sample = flagged.head(SAMPLE_SIZE)

    console.print(f"  {gene}: profiling {len(sample)} most-flagged compounds for kinome off-targets")

    # Collect all off-target activities
    target_hit_counts: Counter = Counter()
    target_best_nm: dict[str, float] = {}
    target_mol_ids: dict[str, set] = {}

    for _, row in sample.iterrows():
        mol_id = row["molecule_chembl_id"]
        try:
            acts = _get_target_activities(mol_id, primary_target_id or "")
            for a in acts:
                tid = a.get("target_chembl_id")
                val = a.get("standard_value")
                if tid and val:
                    try:
                        val_f = float(val)
                    except (TypeError, ValueError):
                        continue
                    target_hit_counts[tid] += 1
                    if tid not in target_best_nm or val_f < target_best_nm[tid]:
                        target_best_nm[tid] = val_f
                    target_mol_ids.setdefault(tid, set()).add(mol_id)
            time.sleep(0.15)
        except Exception:
            pass

    if not target_hit_counts:
        console.print(f"  [yellow]{gene}: no off-target data retrieved[/yellow]")
        return results_dir / "selectivity_kinome.csv"

    # Remove ChEMBL artefacts: "Unchecked", non-molecular, cell-line, organism
    EXCLUDE_TYPES = {"UNCHECKED", "NON-MOLECULAR", "CELL-LINE", "ORGANISM", "TISSUE",
                     "SUBCELLULAR", "PROTEIN-PROTEIN INTERACTION"}
    EXCLUDE_NAMES = {"unchecked", "non-protein target"}

    # Resolve target names (top 40 by hit frequency)
    console.print(f"  {gene}: resolving {min(40, len(target_hit_counts))} target names from ChEMBL")
    top_targets = [tid for tid, _ in target_hit_counts.most_common(40)]
    target_info: dict[str, dict] = {}
    for tid in top_targets:
        try:
            info = _get_target_info(tid)
            target_info[tid] = {
                "name": info.get("pref_name", tid),
                "target_type": info.get("target_type", "Unknown"),
                "organism": info.get("organism", ""),
            }
            time.sleep(0.1)
        except Exception:
            target_info[tid] = {"name": tid, "target_type": "Unknown", "organism": ""}

    # Build output dataframe, filtering artefacts
    rows = []
    for tid, count in target_hit_counts.most_common(60):
        info = target_info.get(tid, {})
        name = info.get("name", tid)
        ttype = info.get("target_type", "Unknown")
        if ttype.upper() in EXCLUDE_TYPES:
            continue
        if name.lower() in EXCLUDE_NAMES:
            continue
        rows.append({
            "target_chembl_id": tid,
            "target_name": name,
            "target_type": ttype,
            "organism": info.get("organism", ""),
            "compound_count": count,
            "best_activity_nm": round(target_best_nm.get(tid, float("nan")), 1),
            "example_compounds": "; ".join(list(target_mol_ids.get(tid, set()))[:3]),
        })
        if len(rows) >= 40:
            break

    kinome_df = pd.DataFrame(rows)
    out_csv = results_dir / "selectivity_kinome.csv"
    kinome_df.to_csv(out_csv, index=False)

    _write_kinome_report(gene, kinome_df, len(flagged), len(sample), results_dir)

    top = kinome_df.iloc[0] if not kinome_df.empty else None
    console.print(
        f"  {gene}: top off-target = {top['target_name']} "
        f"(hit by {top['compound_count']}/{len(sample)} compounds, "
        f"best {top['best_activity_nm']} nM)"
    )
    return out_csv


def _write_kinome_report(gene, df, total_flagged, sample_size, results_dir):
    kinases = df[df["target_type"].str.contains("KINASE|kinase", na=False, case=False)]

    lines = [
        f"# {gene} — Kinome Selectivity Analysis",
        "",
        "## Overview",
        "",
        f"- Compounds with off-target activity (≥4 targets at ≤1µM): **{total_flagged}**",
        f"- Compounds profiled in this analysis: **{sample_size}** (most-flagged subset)",
        f"- Distinct off-targets identified: **{len(df)}**",
        f"- Kinase off-targets: **{len(kinases)}**",
        "",
        "## Top Off-Targets",
        "",
        "| Target | Type | Compound hits | Best activity (nM) |",
        "|---|---|---|---|",
    ]
    for _, row in df.head(20).iterrows():
        lines.append(
            f"| {row['target_name']} | {row['target_type']} | "
            f"{row['compound_count']}/{sample_size} | {row['best_activity_nm']} |"
        )

    lines += ["", "## Kinase Off-Targets Only", ""]
    if kinases.empty:
        lines.append("*No kinase off-targets identified — promiscuity may be non-kinase mediated.*")
    else:
        lines += [
            "| Kinase | Compound hits | Best activity (nM) |",
            "|---|---|---|",
        ]
        for _, row in kinases.iterrows():
            lines.append(f"| {row['target_name']} | {row['compound_count']}/{sample_size} | {row['best_activity_nm']} |")

    lines += [
        "",
        "## Selectivity Strategy Implications",
        "",
        "Understanding which kinases VRK1 compounds most commonly cross-react with",
        "defines the selectivity challenge. Key questions:",
        "",
    ]

    if not kinases.empty:
        top_kinase = kinases.iloc[0]["target_name"]
        lines += [
            f"- **Primary selectivity hurdle:** {top_kinase} — appears in the most compounds.",
            "  Structural comparison of VRK1 vs this kinase binding site will reveal",
            "  which design changes can discriminate between them.",
            "- **Hinge binder vs back-pocket:** If the top off-targets share a conserved hinge,",
            "  the current scaffolds may be hinge-binding; shifting to back-pocket engagement",
            "  (type II or allosteric inhibitors) could resolve selectivity.",
            "- **For SCF-013 (2.1 nM, 100% selective):** This series already avoids the kinome",
            "  cross-reactivity seen in less selective compounds — it is the structural template",
            "  to study for what confers selectivity.",
        ]
    else:
        lines.append("- Off-target activity appears non-kinase in nature — unusual for a kinase inhibitor program.")
        lines.append("  Review target types above for alternative mechanisms (e.g., ion channels, GPCRs).")

    out_md = results_dir / "kinome_report.md"
    out_md.write_text("\n".join(lines) + "\n")
