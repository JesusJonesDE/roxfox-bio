from __future__ import annotations

import statistics
from collections import defaultdict
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from rich.console import Console
    from pipeline.config import Settings

from pipeline.models import Compound, Scaffold


def _murcko_smiles(smiles: str) -> str | None:
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        scaffold = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smiles, includeChirality=False)
        return scaffold if scaffold else None
    except Exception:
        return None


def run_scaffold_analysis(
    gene: str,
    compounds: list[Compound],
    settings: "Settings",
    console: "Console",
) -> list[Compound]:
    if not compounds:
        results_dir = settings.results_dir / gene
        results_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["scaffold_id", "scaffold_smiles", "compound_count", "median_potency_nm", "best_potency_nm"]).to_csv(
            results_dir / "scaffolds.csv", index=False
        )
        console.print(f"  [dim]{gene:10}[/dim] analyze  scaffolds      [yellow]SKIP[/yellow]  (no compounds)")
        return compounds

    # Extract scaffolds
    scaffold_groups: dict[str, list[float]] = defaultdict(list)
    compound_scaffold_map: dict[str, str] = {}
    fails = 0

    for c in compounds:
        smi = _murcko_smiles(c.smiles)
        if smi is None:
            fails += 1
            continue
        scaffold_groups[smi].append(c.best_value_nm)
        compound_scaffold_map[c.molecule_chembl_id] = smi

    # Build scaffold objects with sequential IDs
    scaffolds: list[Scaffold] = []
    smi_to_id: dict[str, str] = {}
    sorted_scaffolds = sorted(scaffold_groups.items(), key=lambda x: -len(x[1]))

    for i, (smi, potencies) in enumerate(sorted_scaffolds, 1):
        sid = f"SCF-{i:03d}"
        smi_to_id[smi] = sid
        scaffolds.append(Scaffold(
            scaffold_smiles=smi,
            scaffold_id=sid,
            compound_count=len(potencies),
            median_potency_nm=round(statistics.median(potencies), 2),
            best_potency_nm=round(min(potencies), 2),
            target_gene=gene,
        ))

    # Update compound scaffold_id references
    for c in compounds:
        smi = compound_scaffold_map.get(c.molecule_chembl_id)
        if smi:
            c.scaffold_id = smi_to_id.get(smi)

    # Write scaffolds CSV
    results_dir = settings.results_dir / gene
    pd.DataFrame([
        {
            "scaffold_id": s.scaffold_id,
            "scaffold_smiles": s.scaffold_smiles,
            "compound_count": s.compound_count,
            "median_potency_nm": s.median_potency_nm,
            "best_potency_nm": s.best_potency_nm,
        }
        for s in scaffolds
    ]).to_csv(results_dir / "scaffolds.csv", index=False)

    # Re-write compounds CSV with updated scaffold_ids
    csv_path = results_dir / "compounds_filtered.csv"
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        scaffold_map_series = {c.molecule_chembl_id: c.scaffold_id for c in compounds}
        df["scaffold_id"] = df["molecule_chembl_id"].map(scaffold_map_series)
        df.to_csv(csv_path, index=False)

    console.print(
        f"  [dim]{gene:10}[/dim] analyze  scaffolds      [green]OK[/green]    "
        f"({len(compounds)} compounds → {len(scaffolds)} unique scaffolds, {fails} extraction fails)"
    )
    return compounds
