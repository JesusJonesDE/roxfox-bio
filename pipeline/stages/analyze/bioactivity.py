from __future__ import annotations

import math
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from rich.console import Console
    from pipeline.cache import CacheManager
    from pipeline.config import Settings

from pipeline.models import Compound

UNIT_TO_NM: dict[str, float] = {
    "nM": 1.0,
    "nm": 1.0,
    "uM": 1_000.0,
    "µM": 1_000.0,
    "um": 1_000.0,
    "mM": 1_000_000.0,
    "M": 1_000_000_000.0,
    "pM": 0.001,
    "pm": 0.001,
}

POTENCY_THRESHOLD_NM = 10_000.0  # 10µM


def _to_nm(value: float, units: str) -> float:
    factor = UNIT_TO_NM.get(units, 1.0)
    return value * factor


def _lipinski(smiles: str) -> dict | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rb = Descriptors.NumRotatableBonds(mol)

        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        return {
            "molecular_weight": round(mw, 2),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rb,
            "ro5_violations": violations,
            "passes_ro5": violations < 2,
        }
    except Exception:
        return None


def run_bioactivity_analysis(
    gene: str,
    settings: "Settings",
    cache: "CacheManager",
    console: "Console",
) -> list[Compound]:
    raw = cache.load(gene, "chembl")
    if not raw:
        console.print(f"  [dim]{gene:10}[/dim] analyze  bioactivity    [yellow]SKIP[/yellow]  (no ChEMBL cache)")
        return []

    activities = raw.get("activities", [])
    if not activities:
        console.print(f"  [dim]{gene:10}[/dim] analyze  bioactivity    [yellow]SKIP[/yellow]  (0 ChEMBL records)")
        _write_empty_compounds_csv(gene, settings)
        return []

    # Filter and normalise
    filtered: dict[str, dict] = {}  # molecule_chembl_id → best record per type
    skipped = 0
    for rec in activities:
        smiles = rec.get("canonical_smiles") or rec.get("molecule_canonical_smiles")
        mol_id = rec.get("molecule_chembl_id")
        std_val = rec.get("standard_value")
        std_type = rec.get("standard_type") or rec.get("activity_type", "")
        std_units = rec.get("standard_units", "nM")

        if not smiles or not mol_id or std_val is None:
            skipped += 1
            continue

        try:
            value_nm = _to_nm(float(std_val), std_units)
        except (ValueError, TypeError):
            skipped += 1
            continue

        if value_nm > POTENCY_THRESHOLD_NM:
            continue

        key = mol_id
        existing = filtered.get(key)
        if existing is None or value_nm < existing["value_nm"]:
            filtered[key] = {
                "molecule_chembl_id": mol_id,
                "smiles": smiles,
                "value_nm": value_nm,
                "assay_type": std_type,
            }

    # Compute Lipinski for each compound
    rdkit_fails = 0
    compounds: list[Compound] = []
    for mol_id, rec in filtered.items():
        props = _lipinski(rec["smiles"])
        if props is None:
            rdkit_fails += 1
            props = {
                "molecular_weight": float("nan"),
                "logp": float("nan"),
                "hbd": -1,
                "hba": -1,
                "rotatable_bonds": -1,
                "ro5_violations": -1,
                "passes_ro5": False,
            }
        compounds.append(Compound(
            molecule_chembl_id=mol_id,
            smiles=rec["smiles"],
            best_value_nm=rec["value_nm"],
            best_assay_type=rec["assay_type"],
            **props,
        ))

    # Write CSV
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / "compounds_filtered.csv"

    rows = [
        {
            "molecule_chembl_id": c.molecule_chembl_id,
            "smiles": c.smiles,
            "best_value_nm": c.best_value_nm,
            "best_assay_type": c.best_assay_type,
            "molecular_weight": c.molecular_weight,
            "logp": c.logp,
            "hbd": c.hbd,
            "hba": c.hba,
            "rotatable_bonds": c.rotatable_bonds,
            "ro5_violations": c.ro5_violations,
            "passes_ro5": c.passes_ro5,
            "scaffold_id": c.scaffold_id,
            "off_target_flags": c.off_target_flags,
            "selectivity_flag": c.selectivity_flag,
        }
        for c in compounds
    ]
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    console.print(
        f"  [dim]{gene:10}[/dim] analyze  bioactivity    [green]OK[/green]    "
        f"({len(activities)} raw → {len(filtered)} pass 10µM → {len(compounds)} compounds, "
        f"{rdkit_fails} RDKit fails)"
    )
    return compounds


def _write_empty_compounds_csv(gene: str, settings: "Settings") -> None:
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)
    columns = [
        "molecule_chembl_id", "smiles", "best_value_nm", "best_assay_type",
        "molecular_weight", "logp", "hbd", "hba", "rotatable_bonds",
        "ro5_violations", "passes_ro5", "scaffold_id", "off_target_flags", "selectivity_flag",
    ]
    pd.DataFrame(columns=columns).to_csv(results_dir / "compounds_filtered.csv", index=False)
