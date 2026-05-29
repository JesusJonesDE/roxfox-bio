from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings


# Module-level thin wrappers — defined here so tests can patch at this module's namespace
# without requiring vina/meeko/biopython to be installed at import time.

def _prepare_receptor(pdb_path: Path, cache_dir: Path) -> Path:
    from pipeline.stages.dock.dock import _prepare_receptor as _impl
    return _impl(pdb_path, cache_dir)


def _prepare_ligand(smiles: str, scaffold_id: str, cache_dir: Path) -> Path:
    from pipeline.stages.dock.dock import _prepare_ligand as _impl
    return _impl(smiles, scaffold_id, cache_dir)


def _run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: list[float],
    box_size: list[float],
    exhaustiveness: int,
    output_pdbqt: Path,
) -> list[dict]:
    from pipeline.stages.dock.dock import _run_vina as _impl
    return _impl(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness, output_pdbqt)


def run_screen(
    gene_symbol: str,
    library_path: Path,
    pocket: dict,
    top_n: int,
    exhaustiveness: int,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> pd.DataFrame:
    """Dock all fragments in library_path against the gene pocket; return top-N hits."""
    # 1. Locate receptor PDB
    pdb_dir = settings.cache_dir / gene_symbol
    pdb_files = list((pdb_dir / "structures").glob("*.pdb")) or list(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(
            f"No cached PDB found for {gene_symbol}. "
            f"Run `pipeline structalign --target {gene_symbol}` (or the AF2 fetch step) first."
        )
    pdb_path = pdb_files[0]

    # 2. Prepare receptor PDBQT (reuse dock cache dir)
    dock_cache_dir = settings.cache_dir / gene_symbol / "dock"
    dock_cache_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt = _prepare_receptor(pdb_path, dock_cache_dir)

    # 3. Docking box from pocket centroid
    center: list[float] = [
        pocket["centroid_x"],
        pocket["centroid_y"],
        pocket["centroid_z"],
    ]
    box_size: list[float] = [pocket["box_size_A"]] * 3  # 20 Å cube

    # 4. Read fragment library (tab-separated: SMILES<TAB>fragment_id)
    fragments: list[tuple[str, str]] = []
    with open(library_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                smiles, fid = parts[0], parts[1]
                fragments.append((smiles, fid))

    total = len(fragments)
    console.print(
        f"  [dim]{gene_symbol}:[/dim] fragment screen — {total} fragments, "
        f"top_n={top_n}, exhaustiveness={exhaustiveness}"
    )

    fragment_cache_dir = settings.cache_dir / gene_symbol / "fragment_dock"
    fragment_cache_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict] = []
    failed = 0

    # 5. Dock each fragment
    for i, (smiles, fid) in enumerate(fragments, start=1):
        cache_key = f"fragment_dock_{fid}"

        # Cache hit check
        if not force:
            cached = cache.load(gene_symbol, cache_key)
            if cached is not None:
                results.append(cached)
                # Progress reporting at 100-fragment boundaries
                if i % 100 == 0:
                    console.print(
                        f"  [dim]{gene_symbol}:[/dim] {i}/{total} fragments docked, {failed} failed"
                    )
                continue

        result: Optional[dict] = None
        try:
            ligand_pdbqt = _prepare_ligand(smiles, fid, fragment_cache_dir)
            output_pdbqt = fragment_cache_dir / f"poses_{fid}.pdbqt"
            poses = _run_vina(
                receptor_pdbqt,
                ligand_pdbqt,
                center,
                box_size,
                exhaustiveness,
                output_pdbqt,
            )
            if poses:
                result = {
                    "fragment_id": fid,
                    "smiles": smiles,
                    "affinity_kcal_mol": poses[0]["affinity_kcal_mol"],
                    "n_poses": len(poses),
                }
            else:
                result = {
                    "fragment_id": fid,
                    "smiles": smiles,
                    "affinity_kcal_mol": 0.0,
                    "n_poses": 0,
                }
        except Exception as exc:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] [yellow]fragment {fid} failed: {exc}[/yellow]"
            )
            failed += 1
            result = None

        if result is not None:
            cache.save(gene_symbol, cache_key, result, 1)
            results.append(result)

        # Progress reporting every 100 fragments
        if i % 100 == 0:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] {i}/{total} fragments docked, {failed} failed"
            )

    # 6. Sort by affinity ascending (most negative first), take top_n
    results.sort(key=lambda r: r["affinity_kcal_mol"])
    top_results = results[:top_n]

    # 7. Write fragment_hits.csv
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)
    hits_df = pd.DataFrame(
        top_results, columns=["fragment_id", "smiles", "affinity_kcal_mol", "n_poses"]
    )
    hits_df.to_csv(results_dir / "fragment_hits.csv", index=False)

    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]fragment screen complete — "
        f"{len(top_results)} hits written to fragment_hits.csv[/green]"
    )

    # 8. Return top hits DataFrame
    return hits_df
