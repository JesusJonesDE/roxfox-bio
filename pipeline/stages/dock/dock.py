from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings

_REMARK_PREFIX = "REMARK VINA RESULT:"


def _check_vina_installed() -> None:
    """Raise RuntimeError with install instructions if vina is not importable."""
    try:
        import vina  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "AutoDock Vina is not installed.\n"
            "Install with: pip install vina meeko rdkit-pypi\n"
            "Or install the docking extras: pip install 'rxpipeline[docking]'"
        )


_AD_ATOM_TYPE_MAP = {
    "C": "C", "N": "N", "O": "OA", "S": "SA", "H": "H",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "MG": "Mg", "ZN": "Zn", "FE": "Fe", "CA": "Ca",
}


def _pdb_to_pdbqt_receptor(pdb_path: Path, pdbqt_path: Path, chain_id: str = "A") -> None:
    """Convert PDB → PDBQT (single chain, protein ATOM records only, zero charge).

    Restricts to one chain to keep the receptor manageable for Vina.
    """
    from Bio.PDB import PDBParser, Select, PDBIO

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))
    first_model = next(iter(structure))

    # Find the chain with the most residues if specified chain not present
    available_chains = [c.id for c in first_model]
    target_chain = chain_id if chain_id in available_chains else available_chains[0]

    lines = []
    atom_idx = 0
    chain = first_model[target_chain]
    for residue in chain:
        if residue.id[0] != " ":
            continue  # skip HETATM / water
        for atom in residue:
            if atom.element == "H":
                continue
            atom_idx += 1
            x, y, z = atom.get_vector()
            element = (atom.element or atom.get_name()[:1]).upper().strip()
            ad_type = _AD_ATOM_TYPE_MAP.get(element, "C")
            name_raw = atom.get_name().strip()
            name_field = f" {name_raw:<3s}" if len(name_raw) < 4 else f"{name_raw:<4s}"
            res_name = residue.resname.ljust(3)
            res_seq = residue.id[1]
            lines.append(
                f"ATOM  {atom_idx:5d} {name_field} {res_name} {target_chain}{res_seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
                f"     0.000 {ad_type:<2s}\n"
            )
    pdbqt_path.write_text("".join(lines))


def _prepare_receptor(pdb_path: Path, cache_dir: Path) -> Path:
    """Convert PDB to PDBQT for AutoDock Vina. Cached."""
    pdbqt_path = cache_dir / f"{pdb_path.stem}_receptor.pdbqt"
    if pdbqt_path.exists():
        return pdbqt_path

    cache_dir.mkdir(parents=True, exist_ok=True)
    _pdb_to_pdbqt_receptor(pdb_path, pdbqt_path)
    return pdbqt_path


def _prepare_ligand(smiles: str, scaffold_id: str, cache_dir: Path) -> Path:
    """Convert SMILES → 3D conformer → PDBQT using RDKit + meeko. Cached."""
    pdbqt_path = cache_dir / f"{scaffold_id}_ligand.pdbqt"
    if pdbqt_path.exists():
        return pdbqt_path

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation
    except ImportError:
        raise RuntimeError(
            "RDKit and meeko are required for ligand preparation.\n"
            "Install with: pip install vina meeko rdkit-pypi"
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES for {scaffold_id}: {smiles!r}")

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        raise RuntimeError(f"RDKit EmbedMolecule failed for {scaffold_id}")
    AllChem.UFFOptimizeMolecule(mol)

    cache_dir.mkdir(parents=True, exist_ok=True)
    # meeko >= 0.5: prepare() returns list of RDKitMoleculeSetup; write via PDBQTWriterLegacy
    from meeko import PDBQTWriterLegacy
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise RuntimeError(f"meeko PDBQT write failed for {scaffold_id}: {err}")
    pdbqt_path.write_text(pdbqt_str)
    return pdbqt_path


def _define_box(
    binding_site_residues: list[dict],
    margin_A: float = 3.0,
) -> tuple[list[float], list[float]]:
    """Compute docking box centroid and size from Cα coordinates.

    Returns (center_xyz, box_size_xyz) where box_size is 2*(r+margin) on each axis.
    """
    xs = [r["ca_x"] for r in binding_site_residues]
    ys = [r["ca_y"] for r in binding_site_residues]
    zs = [r["ca_z"] for r in binding_site_residues]

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)

    # Max radius from centroid
    r = max(
        ((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) ** 0.5
        for x, y, z in zip(xs, ys, zs)
    )
    side = 2 * (r + margin_A)
    return [cx, cy, cz], [side, side, side]


def _run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: list[float],
    box_size: list[float],
    exhaustiveness: int,
    output_pdbqt: Path,
) -> list[dict]:
    """Run AutoDock Vina and parse REMARK VINA RESULT lines.

    Returns list of {pose_rank, affinity_kcal_mol, rmsd_lb, rmsd_ub}.
    """
    from vina import Vina

    v = Vina(sf_name="vina")
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(
        center=center,
        box_size=box_size,
    )
    v.dock(exhaustiveness=exhaustiveness, n_poses=9)
    v.write_poses(str(output_pdbqt), n_poses=9, overwrite=True)

    # Parse scores
    poses = []
    with open(output_pdbqt) as fh:
        for line in fh:
            if line.startswith(_REMARK_PREFIX):
                parts = line[len(_REMARK_PREFIX):].split()
                if len(parts) >= 3:
                    poses.append({
                        "pose_rank": len(poses) + 1,
                        "affinity_kcal_mol": float(parts[0]),
                        "rmsd_lb": float(parts[1]),
                        "rmsd_ub": float(parts[2]),
                    })
    return poses


def _run_control(
    pdb_path: Path,
    receptor_pdbqt: Path,
    center: list[float],
    box_size: list[float],
    cache_dir: Path,
) -> float:
    """Re-dock crystallographic ANP ligand from 6AC9 and return heavy-atom RMSD.

    RMSD > 2.0 Å is flagged as a warning but does not abort docking.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio import SVDSuperimposer
        import numpy as np
    except ImportError:
        raise RuntimeError("BioPython is required: pip install biopython")

    # Extract ANP (adenylyl imidodiphosphate) HETATM atoms from PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))
    ref_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname in ("ANP", "ADP", "ATP", "AMP"):
                    for atom in residue:
                        if atom.element != "H":
                            ref_coords.append(atom.get_vector().get_array())
    if not ref_coords:
        return float("nan")

    import numpy as np
    ref_arr = np.array(ref_coords)

    # Write reference ligand PDBQT (simplified: use reference coords as-is)
    ref_pdbqt = cache_dir / "control_ligand.pdbqt"
    _write_xyz_as_pdbqt(ref_arr, ref_pdbqt)

    output_pdbqt = cache_dir / "control_redock.pdbqt"
    poses = _run_vina(receptor_pdbqt, ref_pdbqt, center, box_size, 16, output_pdbqt)
    if not poses:
        return float("nan")

    # Parse top-pose coords
    top_coords = _parse_pdbqt_coords(output_pdbqt, pose_index=0)
    if top_coords is None or len(top_coords) != len(ref_arr):
        return float("nan")

    diffs = ref_arr - top_coords
    rmsd = float(np.sqrt((diffs ** 2).sum(axis=1).mean()))
    return rmsd


def _write_xyz_as_pdbqt(coords, path: Path) -> None:
    """Write bare-bones rigid PDBQT from an array of (x, y, z) heavy-atom coordinates."""
    lines = ["ROOT\n"]
    for i, (x, y, z) in enumerate(coords, 1):
        lines.append(
            f"HETATM{i:5d}  C   LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
            f"     0.000 C \n"
        )
    lines.append("ENDROOT\n")
    lines.append("TORSDOF 0\n")
    path.write_text("".join(lines))


def _parse_pdbqt_coords(pdbqt_path: Path, pose_index: int = 0):
    """Extract heavy-atom coordinates from the nth MODEL in a PDBQT file."""
    import numpy as np
    coords = []
    current_model = -1
    in_target = False
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                current_model += 1
                in_target = current_model == pose_index
                continue
            if line.startswith("ENDMDL"):
                if in_target:
                    break
                continue
            if in_target and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    element = line[76:78].strip()
                    if element != "H":
                        coords.append([x, y, z])
                except ValueError:
                    pass
    return np.array(coords) if coords else None


def _map_contacts(
    pose_pdbqt_path: Path,
    comparison_csv_path: Path,
    cutoff_A: float = 4.0,
) -> list[int]:
    """Find selectivity-candidate KLIFS positions within cutoff_A of top pose heavy atoms.

    Returns list of contacted KLIFS positions.
    """
    import numpy as np

    pose_coords = _parse_pdbqt_coords(pose_pdbqt_path, pose_index=0)
    if pose_coords is None or len(pose_coords) == 0:
        return []

    df = pd.read_csv(comparison_csv_path)
    candidates = df[df["selectivity_candidate"] == True]  # noqa: E712
    # binding_site_vrk1.csv uses ca_x/ca_y/ca_z columns
    if not {"ca_x", "ca_y", "ca_z", "klifs_position"}.issubset(df.columns):
        return []

    contacted = []
    for _, row in candidates.iterrows():
        ca = np.array([row["ca_x"], row["ca_y"], row["ca_z"]])
        dists = np.linalg.norm(pose_coords - ca, axis=1)
        if dists.min() <= cutoff_A:
            contacted.append(int(row["klifs_position"]))
    return contacted


def _extract_binding_site_coords(
    pdb_path: Path,
    binding_site_csv: Path,
) -> list[dict]:
    """Extract Cα coordinates for residues listed in binding_site_vrk1.csv.

    Falls back to ligand heavy-atom centroid if no residues can be matched.
    """
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))
    bs_df = pd.read_csv(binding_site_csv)

    # Build residue number → Cα coord mapping
    ca_map: dict[int, tuple[float, float, float]] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue  # skip HETATM
                res_num = residue.id[1]
                if "CA" in residue:
                    pos = residue["CA"].get_vector()
                    ca_map[res_num] = (pos[0], pos[1], pos[2])

    coords = []
    for _, row in bs_df.iterrows():
        res_id = row.get("residue.id")
        try:
            res_num = int(float(str(res_id).split(".")[0]))
        except (ValueError, TypeError):
            continue
        if res_num in ca_map:
            x, y, z = ca_map[res_num]
            coords.append({"ca_x": x, "ca_y": y, "ca_z": z})

    if coords:
        return coords

    # Fallback: use ANP/ligand heavy-atom centroid
    lig_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname in ("ANP", "ADP", "ATP", "AMP", "LIG"):
                    for atom in residue:
                        if atom.element != "H":
                            pos = atom.get_vector()
                            lig_coords.append({"ca_x": pos[0], "ca_y": pos[1], "ca_z": pos[2]})
    return lig_coords


def _write_report(
    gene_symbol: str,
    scaffold_id: str,
    smiles: str,
    poses: list[dict],
    control_rmsd: float,
    contacted_positions: list[int],
    comparison_df: pd.DataFrame,
    results_dir: Path,
) -> Path:
    control_str = (
        f"{control_rmsd:.2f} Å" if control_rmsd == control_rmsd else "N/A (extraction failed)"
    )
    control_warn = (
        "\n\n> **Warning**: Control RMSD > 2.0 Å — docking box may be misaligned. "
        "Results should be interpreted with caution."
        if (control_rmsd == control_rmsd and control_rmsd > 2.0)
        else ""
    )

    pose_rows = [
        f"| {p['pose_rank']} | {p['affinity_kcal_mol']:.1f} | {p['rmsd_lb']:.2f} | {p['rmsd_ub']:.2f} |"
        for p in poses[:3]
    ]

    # Selectivity candidates from comparison CSV
    cands = comparison_df[comparison_df["selectivity_candidate"] == True] if "selectivity_candidate" in comparison_df.columns else pd.DataFrame()  # noqa: E712
    contact_rows = []
    for _, r in cands.iterrows():
        kpos = int(r["klifs_position"])
        hit = "YES" if kpos in contacted_positions else "—"
        contact_rows.append(
            f"| {kpos} | {r.get('subpocket', '?')} | {r.get('vrk1_aa', '?')} "
            f"| {r.get('egfr_aa', '?')} | {hit} |"
        )

    report = f"""# {gene_symbol} Docking Report — {scaffold_id}

**Generated**: {datetime.now().strftime("%Y-%m-%d")}
**Scaffold**: {scaffold_id} | `{smiles}`
**Method**: AutoDock Vina (local), meeko ligand prep, KLIFS binding box
**Control**: ANP re-dock RMSD = {control_str}{control_warn}

---

## Top 3 Docking Poses

| Pose | Affinity (kcal/mol) | RMSD lower | RMSD upper |
|------|---------------------|-----------|-----------|
{chr(10).join(pose_rows) if pose_rows else "| — | No poses returned | — | — |"}

---

## Selectivity Candidate Contacts (≤ 4.0 Å from Cα)

Positions from spec-003 structural alignment that differ between VRK1 and EGFR.

| KLIFS Pos | Subpocket | VRK1 | EGFR | Contact? |
|-----------|-----------|------|------|----------|
{chr(10).join(contact_rows) if contact_rows else "| — | Comparison CSV not available | — | — | — |"}

---

## Interpretation

{f"Top pose affinity: {poses[0]['affinity_kcal_mol']:.1f} kcal/mol. " if poses else "No docking poses were generated. "}
{f"{len(contacted_positions)} of {len(cands)} selectivity candidates contacted." if not cands.empty else ""}
{"Strong binding predicted — proceed to analog design focusing on contacted selectivity handles." if poses and poses[0]["affinity_kcal_mol"] < -7.0 else "Moderate binding — scaffold may require optimization."}
"""
    path = results_dir / f"docking_report_{scaffold_id}.md"
    path.write_text(report)
    return path


def run_dock(
    gene_symbol: str,
    scaffold_id: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    exhaustiveness: int,
    console: Console,
) -> None:
    _check_vina_installed()

    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    cache_key = f"dock_{scaffold_id}"
    if not force:
        cached = cache.load(gene_symbol, cache_key)
        if cached is not None:
            count = len(cached) if isinstance(cached, list) else "?"
            console.print(
                f"  [dim]{gene_symbol:10}[/dim] dock {scaffold_id:12} "
                f"[yellow]SKIP[/yellow]  (cached {count} poses)"
            )
            return

    # Resolve scaffold SMILES from compounds_filtered.csv
    compounds_csv = results_dir / "compounds_filtered.csv"
    if not compounds_csv.exists():
        raise FileNotFoundError(
            f"compounds_filtered.csv not found in {results_dir}. "
            f"Run `pipeline fetch --target {gene_symbol}` first."
        )
    compounds = pd.read_csv(compounds_csv)
    row = compounds[compounds["scaffold_id"] == scaffold_id]
    if row.empty:
        # Try 'compound_id' or 'name' column fallback
        for col in ("compound_id", "name", "molecule_chembl_id"):
            if col in compounds.columns:
                row = compounds[compounds[col] == scaffold_id]
                if not row.empty:
                    break
    if row.empty:
        raise ValueError(f"Scaffold '{scaffold_id}' not found in compounds_filtered.csv")
    smiles = str(row.iloc[0].get("smiles") or row.iloc[0].get("canonical_smiles") or row.iloc[0].get("SMILES"))

    # Resolve receptor PDB from structalign outputs (check structures/ subdir first)
    pdb_dir = settings.cache_dir / gene_symbol
    pdb_files = list((pdb_dir / "structures").glob("*.pdb")) or list(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(
            f"No cached PDB found for {gene_symbol}. "
            f"Run `pipeline structalign --target {gene_symbol}` first."
        )
    pdb_path = pdb_files[0]

    dock_cache_dir = settings.cache_dir / gene_symbol / "dock"
    dock_cache_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"  [dim]{gene_symbol}:[/dim] preparing receptor {pdb_path.name}...")
    receptor_pdbqt = _prepare_receptor(pdb_path, dock_cache_dir)

    console.print(f"  [dim]{gene_symbol}:[/dim] preparing ligand {scaffold_id}...")
    ligand_pdbqt = _prepare_ligand(smiles, scaffold_id, dock_cache_dir)

    # Load binding site residues for box definition
    binding_site_csv = results_dir / "binding_site_vrk1.csv"
    if not binding_site_csv.exists():
        raise FileNotFoundError(
            f"binding_site_vrk1.csv not found. "
            f"Run `pipeline structalign --target {gene_symbol}` first."
        )
    console.print(f"  [dim]{gene_symbol}:[/dim] extracting binding site Cα coordinates...")
    binding_residues = _extract_binding_site_coords(pdb_path, binding_site_csv)
    if not binding_residues:
        raise RuntimeError(
            f"Could not extract binding site coordinates for {gene_symbol}. "
            f"Check that {pdb_path} contains the expected residues."
        )
    center, box_size = _define_box(binding_residues)

    console.print(f"  [dim]{gene_symbol}:[/dim] running ANP control docking...")
    try:
        control_rmsd = _run_control(pdb_path, receptor_pdbqt, center, box_size, dock_cache_dir)
    except Exception as e:
        console.print(f"  [dim]{gene_symbol}:[/dim] [yellow]control docking failed[/yellow]: {e}")
        control_rmsd = float("nan")
    if control_rmsd == control_rmsd and control_rmsd > 2.0:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [yellow]WARNING[/yellow] control RMSD {control_rmsd:.2f} Å > 2.0 Å"
        )

    output_pdbqt = results_dir / f"docking_poses_{scaffold_id}.pdbqt"
    console.print(
        f"  [dim]{gene_symbol}:[/dim] docking {scaffold_id} (exhaustiveness={exhaustiveness})..."
    )
    poses = _run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness, output_pdbqt)

    # Map contacts to selectivity candidates
    # Merge Cα coords (from PDB) into comparison CSV for contact mapping
    comparison_csv = results_dir / "binding_site_comparison.csv"
    comparison_df = pd.read_csv(comparison_csv) if comparison_csv.exists() else pd.DataFrame()

    contacted: list[int] = []
    if comparison_csv.exists() and binding_residues:
        # Build residue number → coords map from the PDB extraction
        bs_df = pd.read_csv(binding_site_csv)
        coord_list = _extract_binding_site_coords(pdb_path, binding_site_csv)
        # Attach coords to comparison via residue.id matching
        res_ids = bs_df["residue.id"].tolist()
        coord_map: dict[float, dict] = {}
        for rid, coord in zip(res_ids, coord_list):
            try:
                coord_map[float(rid)] = coord
            except (ValueError, TypeError):
                pass

        # Merge klifs_position → residue.id → coords into comparison_df
        if "klifs_position" in comparison_df.columns and "residue.id" in bs_df.columns:
            merged = comparison_df.merge(
                bs_df[["residue.klifs_id", "residue.id"]].rename(
                    columns={"residue.klifs_id": "klifs_position"}
                ),
                on="klifs_position",
                how="left",
            )
            merged["ca_x"] = merged["residue.id"].map(
                lambda r: coord_map.get(float(r) if r else None, {}).get("ca_x")  # type: ignore[arg-type]
            )
            merged["ca_y"] = merged["residue.id"].map(
                lambda r: coord_map.get(float(r) if r else None, {}).get("ca_y")  # type: ignore[arg-type]
            )
            merged["ca_z"] = merged["residue.id"].map(
                lambda r: coord_map.get(float(r) if r else None, {}).get("ca_z")  # type: ignore[arg-type]
            )
            tmp_csv = results_dir / "_tmp_comparison_coords.csv"
            merged.to_csv(tmp_csv, index=False)
            contacted = _map_contacts(output_pdbqt, tmp_csv)
            tmp_csv.unlink(missing_ok=True)
            comparison_df = merged

    # Save cache
    cache.save(gene_symbol, cache_key, poses, len(poses))

    # Write CSV and report
    pd.DataFrame(poses).to_csv(results_dir / f"docking_results_{scaffold_id}.csv", index=False)
    _write_report(gene_symbol, scaffold_id, smiles, poses, control_rmsd, contacted, comparison_df, results_dir)

    top_affinity = poses[0]["affinity_kcal_mol"] if poses else float("nan")
    console.print(
        f"  [dim]{gene_symbol}:[/dim] {len(poses)} poses | "
        f"top affinity: {top_affinity:.1f} kcal/mol | "
        f"control RMSD: {control_rmsd:.2f} Å | "
        f"{len(contacted)} selectivity candidates contacted"
    )
    console.print(f"  [dim]{gene_symbol}:[/dim] [green]docking_report_{scaffold_id}.md written[/green]")
