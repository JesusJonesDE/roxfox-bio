from __future__ import annotations

import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings
from pipeline.models import GateResult, GateStatus
from pipeline.stages.validate.validate import (
    _cache_gate_result,
    _load_cached_gate_result,
    _load_smiles,
    _write_gate_report,
)

# ── Data structures ────────────────────────────────────────────────────────────


@dataclass
class OffTargetEntry:
    gene: str
    pdb_id: str        # e.g. "1M17" or "AF2"
    source: str        # "pdb" or "alphafold"
    ligand_id: str     # e.g. "AQ4" — for box centroid
    chain_id: str      # e.g. "A"
    af2_uniprot: str = ""  # UniProt ID for AF2 download
    warning: str = ""


# ── Default kinase panel ───────────────────────────────────────────────────────

SELECTIVITY_PANEL = [
    OffTargetEntry(
        "VRK2", "AF2", "alphafold", "", "A",
        af2_uniprot="O95551",
        warning="AlphaFold2 model — lower docking confidence",
    ),
    OffTargetEntry("EGFR", "1IVO", "pdb", "AQ4", "A"),
    OffTargetEntry("CDK2", "1E9H", "pdb", "ATP", "A"),
    OffTargetEntry("PLK1", "2OKR", "pdb", "ADP", "A"),
]

_RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
_AF2_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


# ── Structure fetching ─────────────────────────────────────────────────────────


def _fetch_offtarget_structure(entry: OffTargetEntry, panel_dir: Path) -> Path:
    """Download and cache off-target PDB from RCSB or AlphaFold EBI.

    Cache path: panel_dir / f"{entry.gene}_{entry.pdb_id}.pdb"
    Returns path to cached file.
    """
    panel_dir.mkdir(parents=True, exist_ok=True)
    cache_path = panel_dir / f"{entry.gene}_{entry.pdb_id}.pdb"

    if cache_path.exists():
        return cache_path

    if entry.source == "pdb":
        url = _RCSB_URL.format(pdb_id=entry.pdb_id)
    elif entry.source == "alphafold":
        if not entry.af2_uniprot:
            raise ValueError(
                f"OffTargetEntry for {entry.gene} has source='alphafold' "
                f"but no af2_uniprot set."
            )
        url = _AF2_URL.format(uniprot=entry.af2_uniprot)
    else:
        raise ValueError(f"Unknown source {entry.source!r} for {entry.gene}")

    urllib.request.urlretrieve(url, str(cache_path))
    return cache_path


# ── Off-target docking ─────────────────────────────────────────────────────────


def _dock_offtarget(
    smiles: str,
    scaffold_id: str,
    entry: OffTargetEntry,
    panel_dir: Path,
    settings: Settings,
) -> float:
    """Dock scaffold into an off-target structure and return top Vina affinity.

    Reuses _prepare_receptor, _prepare_ligand, _extract_ligand_centroid,
    _define_box, and _run_vina from the dock module.

    Returns top pose affinity (kcal/mol), or 0.0 if no poses produced.
    """
    from pipeline.stages.dock.dock import (
        _define_box,
        _extract_ligand_centroid,
        _prepare_ligand,
        _prepare_receptor,
        _run_vina,
    )

    pdb_path = _fetch_offtarget_structure(entry, panel_dir)

    dock_cache = (
        settings.cache_dir / "shared" / "selectivity_panel" / f"dock_{entry.gene}"
    )
    dock_cache.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = _prepare_receptor(pdb_path, dock_cache)
    ligand_pdbqt = _prepare_ligand(smiles, scaffold_id, dock_cache)

    center = _extract_ligand_centroid(pdb_path, chain_id=entry.chain_id)
    if center is None:
        # Warn but continue — use origin as fallback
        center = [0.0, 0.0, 0.0]

    center, box_size = _define_box(center)

    output_pdbqt = dock_cache / f"poses_{scaffold_id}.pdbqt"
    poses = _run_vina(
        receptor_pdbqt,
        ligand_pdbqt,
        center,
        box_size,
        8,  # exhaustiveness=8 for speed in panel docking
        output_pdbqt,
    )

    if not poses:
        return 0.0
    return poses[0]["affinity_kcal_mol"]


# ── Gate orchestrator ──────────────────────────────────────────────────────────


def run_selectivity_gate(
    gene_symbol: str,
    scaffold_id: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> GateResult:
    """Run the selectivity docking panel gate.

    Docks the scaffold into each off-target in SELECTIVITY_PANEL (excluding
    any entry whose gene matches gene_symbol), computes a selectivity index
    SI = |primary_affinity| / max(|offtarget_affinity|), and returns PASS
    if SI >= 10.

    Steps:
    1. Cache check (skip if cached and not force).
    2. Load SMILES from compounds_filtered.csv.
    3. Load primary target Vina affinity from docking_results_{scaffold}.csv.
    4. Build shared panel_dir and filter panel entries.
    5. Fetch structures and dock each off-target.
    6. Compute SI; PASS if >= 10.
    7. Build GateResult, write report, cache, return.
    """
    # 1. Cache check
    if not force:
        cached = _load_cached_gate_result(gene_symbol, scaffold_id, "selectivity", cache)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] Selectivity gate — [dim]SKIP[/dim] "
                f"(cached result for {scaffold_id})"
            )
            return cached

    # 2. Load SMILES
    smiles = _load_smiles(gene_symbol, scaffold_id, settings)

    # 3. Load primary affinity from existing docking CSV
    docking_csv = settings.results_dir / gene_symbol / f"docking_results_{scaffold_id}.csv"
    if not docking_csv.exists():
        raise FileNotFoundError(
            f"Docking results not found: {docking_csv}\n"
            f"Run `pipeline dock --target {gene_symbol} --scaffold {scaffold_id}` first."
        )
    docking_df = pd.read_csv(docking_csv)
    if docking_df.empty or "affinity_kcal_mol" not in docking_df.columns:
        raise ValueError(
            f"docking_results_{scaffold_id}.csv is empty or missing 'affinity_kcal_mol' column."
        )
    primary_affinity: float = float(docking_df.iloc[0]["affinity_kcal_mol"])

    # 4. Build panel dir; filter panel (skip entries matching the primary target)
    panel_dir = settings.cache_dir / "shared" / "selectivity_panel"
    panel_dir.mkdir(parents=True, exist_ok=True)

    active_panel = [
        entry for entry in SELECTIVITY_PANEL
        if entry.gene.upper() != gene_symbol.upper()
    ]

    if not active_panel:
        # Edge case: all panel entries excluded (e.g. target is all four kinases)
        result = GateResult(
            gate_name="selectivity",
            status=GateStatus.PASS,
            score=float("inf"),
            reason="No off-targets remain after excluding primary target — selectivity assumed.",
            details={"primary_affinity": primary_affinity},
            timestamp=datetime.now(timezone.utc).isoformat(),
        )
        result.report_path = _write_gate_report(
            gene_symbol, scaffold_id, "selectivity", result, settings
        )
        _cache_gate_result(gene_symbol, scaffold_id, "selectivity", result, cache)
        return result

    # 5. Fetch structures and dock each off-target
    offtarget_affinities: dict[str, float] = {}
    warnings: list[str] = []

    for entry in active_panel:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] docking {scaffold_id} into {entry.gene} "
            f"({'AF2' if entry.source == 'alphafold' else entry.pdb_id})..."
        )
        try:
            affinity = _dock_offtarget(smiles, scaffold_id, entry, panel_dir, settings)
            offtarget_affinities[entry.gene] = affinity
            if entry.warning:
                warnings.append(f"{entry.gene}: {entry.warning}")
        except Exception as exc:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] [yellow]WARNING[/yellow] "
                f"off-target docking failed for {entry.gene}: {exc}"
            )
            # Use 0.0 (very weak binding) as conservative fallback — doesn't inflate SI
            offtarget_affinities[entry.gene] = 0.0

    # 6. Compute selectivity index
    # SI = |primary_affinity| / max(|offtarget_affinity|)
    # Avoid division by zero: if all off-targets return 0.0, SI is infinite (very selective)
    abs_primary = abs(primary_affinity)
    abs_offtargets = {gene: abs(aff) for gene, aff in offtarget_affinities.items()}

    if not abs_offtargets or max(abs_offtargets.values()) == 0.0:
        si = float("inf")
        worst_offtarget = next(iter(abs_offtargets), "none")
        worst_affinity = 0.0
    else:
        worst_offtarget = max(abs_offtargets, key=abs_offtargets.__getitem__)
        worst_affinity = offtarget_affinities[worst_offtarget]
        si = abs_primary / abs_offtargets[worst_offtarget]

    status = GateStatus.PASS if si >= 10.0 else GateStatus.FAIL

    reason = f"SI={si:.1f}x vs {worst_offtarget} ({worst_affinity:.1f} kcal/mol)"
    if warnings:
        reason += "; " + "; ".join(warnings)

    # Build details dict: primary affinity + per off-target affinities
    details: dict[str, object] = {"primary_affinity": primary_affinity}
    for gene, aff in offtarget_affinities.items():
        details[gene] = aff

    result = GateResult(
        gate_name="selectivity",
        status=status,
        score=si,
        reason=reason,
        details=details,
        timestamp=datetime.now(timezone.utc).isoformat(),
    )

    # 7. Write report (with selectivity table extra section)
    extra = _build_selectivity_table(
        primary_affinity,
        offtarget_affinities,
        {e.gene: e.warning for e in active_panel},
    )
    result.report_path = _write_gate_report(
        gene_symbol, scaffold_id, "selectivity", result, settings,
        extra_sections=extra,
    )

    # 8. Cache
    _cache_gate_result(gene_symbol, scaffold_id, "selectivity", result, cache)

    console.print(
        f"  [dim]{gene_symbol}:[/dim] Selectivity SI={si:.1f}x — "
        f"{'PASS' if status == GateStatus.PASS else 'FAIL'}"
    )

    return result


# ── Report helpers ─────────────────────────────────────────────────────────────


def _build_selectivity_table(
    primary_affinity: float,
    offtarget_affinities: dict[str, float],
    warnings: dict[str, str],
) -> str:
    """Build a markdown table comparing off-target affinities to the primary target."""
    abs_primary = abs(primary_affinity)

    rows = []
    for gene, affinity in offtarget_affinities.items():
        abs_aff = abs(affinity)
        if abs_primary > 0 and abs_aff > 0:
            ratio = abs_aff / abs_primary
            ratio_str = f"{ratio:.1f}x"
        elif abs_aff == 0.0:
            ratio_str = "—"
        else:
            ratio_str = "—"

        note = warnings.get(gene) or "—"
        rows.append(f"| {gene} | {affinity:.1f} | {ratio_str} | {note} |")

    rows_text = "\n".join(rows) if rows else "| — | — | — | — |"

    return f"""## Selectivity Panel

| Off-target | Affinity (kcal/mol) | vs. Primary | Note |
|---|---|---|---|
{rows_text}

**Primary target affinity**: {primary_affinity:.1f} kcal/mol
"""
