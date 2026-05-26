from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from rich.console import Console
    from pipeline.cache import CacheManager
    from pipeline.config import Settings

from pipeline.models import Structure

METHOD_MAP = {
    "X-RAY DIFFRACTION": "X-ray",
    "ELECTRON MICROSCOPY": "Cryo-EM",
    "SOLUTION NMR": "NMR",
    "SOLID-STATE NMR": "NMR",
    "NEUTRON DIFFRACTION": "Neutron",
}


def run_structures_analysis(
    gene: str,
    settings: "Settings",
    cache: "CacheManager",
    console: "Console",
) -> list[Structure]:
    pdb_cache = cache.load(gene, "pdb")
    af_cache = cache.load(gene, "alphafold")

    structures: list[Structure] = []

    # PDB structures
    if pdb_cache:
        for entry in pdb_cache.get("structures", []):
            if "error" in entry:
                continue
            raw_method = (entry.get("method") or "Unknown").upper()
            method = METHOD_MAP.get(raw_method, raw_method.title())

            res = entry.get("resolution_angstrom")
            try:
                res = float(res) if res is not None else None
            except (ValueError, TypeError):
                res = None

            structures.append(Structure(
                structure_id=entry["pdb_id"],
                source="PDB",
                method=method,
                has_ligand=bool(entry.get("has_ligand", False)),
                resolution_angstrom=res,
                ligand_ids=entry.get("ligand_ids", []),
                chain_ids=entry.get("chain_ids", []),
                deposition_date=entry.get("deposition_date"),
                target_uniprot="",
            ))

    # Sort PDB by resolution (best first; None last)
    structures.sort(key=lambda s: (s.resolution_angstrom is None, s.resolution_angstrom or 999))

    # AlphaFold model
    if af_cache and af_cache.get("entry_id"):
        structures.append(Structure(
            structure_id=af_cache.get("entry_id", f"AF-{gene}"),
            source="AlphaFold",
            method="Predicted",
            has_ligand=False,
            mean_plddt=af_cache.get("mean_plddt"),
        ))

    # Write CSV
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "structure_id": s.structure_id,
            "source": s.source,
            "resolution_angstrom": s.resolution_angstrom,
            "method": s.method,
            "has_ligand": s.has_ligand,
            "ligand_ids": ",".join(s.ligand_ids),
            "chain_ids": ",".join(s.chain_ids),
            "mean_plddt": s.mean_plddt,
            "deposition_date": s.deposition_date,
        }
        for s in structures
    ]
    pd.DataFrame(rows).to_csv(results_dir / "structures.csv", index=False)

    pdb_count = sum(1 for s in structures if s.source == "PDB")
    af_present = any(s.source == "AlphaFold" for s in structures)
    console.print(
        f"  [dim]{gene:10}[/dim] analyze  structures     [green]OK[/green]    "
        f"({pdb_count} PDB + {'1 AlphaFold' if af_present else 'no AlphaFold'})"
    )
    return structures
