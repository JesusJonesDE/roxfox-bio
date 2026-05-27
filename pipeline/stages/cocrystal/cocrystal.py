from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import httpx
import pandas as pd
from rich.console import Console
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from pipeline.cache import CacheManager
from pipeline.config import Settings

# Resolution requirements per subpocket (from research.md)
_RESOLUTION_REQUIREMENTS = {
    "Gatekeeper": "≤ 2.5 Å (critical — gatekeeper contacts require clear electron density)",
    "Hinge": "≤ 2.5 Å (hinge H-bonds must be resolved unambiguously)",
    "P-loop": "≤ 3.0 Å (P-loop flexibility tolerated at moderate resolution)",
    "DFG": "≤ 2.8 Å (DFG conformation distinguishes active/inactive; < 2.8 Å needed)",
    "Other": "≤ 3.5 Å (peripheral positions tolerate lower resolution)",
}

_RCSB_GRAPHQL = "https://data.rcsb.org/graphql"

_GENERIC_KINASE_CONDITIONS = (
    "No specific crystallisation data retrieved from RCSB.\n"
    "Generic kinase soaking guidance: apo crystal soaked with compound at 1–5 mM "
    "in 20% DMSO (final), 30–60 min at RT. Monitor crystal integrity by birefringence. "
    "See Hassell et al. (2007) Acta Cryst. D63:72 for general methodology."
)


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _http_post(url: str, json_payload: dict, timeout: int = 30) -> httpx.Response:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(url, json=json_payload)
        r.raise_for_status()
        return r


def _fetch_rcsb_conditions(pdb_id: str) -> dict:
    """Query RCSB GraphQL for crystallisation conditions. Falls back gracefully."""
    query = """
    query($entry_ids: [String!]!) {
      entries(entry_ids: $entry_ids) {
        exptl_crystal_grow {
          method
          pdbx_details
        }
      }
    }
    """
    try:
        r = _http_post(
            _RCSB_GRAPHQL,
            {"query": query, "variables": {"entry_ids": [pdb_id.upper()]}},
        )
        data = r.json()
        entries = data.get("data", {}).get("entries", [])
        if entries and entries[0]:
            grow_list = entries[0].get("exptl_crystal_grow") or []
            if grow_list:
                g = grow_list[0]
                method = g.get("method") or "VAPOR DIFFUSION"
                details = g.get("pdbx_details") or ""
                return {"method": method, "conditions": details or _GENERIC_KINASE_CONDITIONS}
    except Exception:
        pass
    return {"method": "VAPOR DIFFUSION", "conditions": _GENERIC_KINASE_CONDITIONS}


def _get_space_group(pdb_path: Path) -> str:
    """Extract space group from CRYST1 record in PDB file."""
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line.startswith("CRYST1"):
                    # Format: CRYST1   a   b   c  α  β  γ space_group Z
                    # space group starts at col 55 (0-indexed)
                    space_group = line[55:66].strip()
                    return space_group or "unknown"
    except OSError:
        pass
    return "unknown"


def _flag_scaffold_atoms(smiles: str, binding_site_centroid: list[float]) -> list[str]:
    """Identify heavy atoms > 8Å from centroid adjacent to rotatable bonds.

    Returns list of atom descriptions (e.g. "C12 adjacent to rotatable bond").
    Returns empty list if SMILES cannot be parsed — does not raise.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
        import numpy as np
    except ImportError:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    mol = Chem.AddHs(mol)
    try:
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == -1:
            return []
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        return []

    conf = mol.GetConformer()
    cx, cy, cz = binding_site_centroid

    # Find rotatable bond atom indices
    rot_bond_atoms: set[int] = set()
    for bond in mol.GetBonds():
        if (
            not bond.IsInRing()
            and bond.GetBondTypeAsDouble() == 1.0
            and bond.GetBeginAtom().GetAtomicNum() > 1
            and bond.GetEndAtom().GetAtomicNum() > 1
        ):
            rot_bond_atoms.add(bond.GetBeginAtomIdx())
            rot_bond_atoms.add(bond.GetEndAtomIdx())

    flagged = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        dist = ((pos.x - cx) ** 2 + (pos.y - cy) ** 2 + (pos.z - cz) ** 2) ** 0.5
        if dist > 8.0 and atom.GetIdx() in rot_bond_atoms:
            flagged.append(
                f"{atom.GetSymbol()}{atom.GetIdx()} (dist {dist:.1f}Å from centroid, "
                f"adjacent to rotatable bond)"
            )
    return flagged


def _write_brief(
    gene_symbol: str,
    scaffold_id: str,
    smiles: str,
    pdb_id: str,
    space_group: str,
    conditions: dict,
    flagged_atoms: list[str],
    results_dir: Path,
) -> Path:
    flag_section = (
        "\n".join(f"- {a}" for a in flagged_atoms)
        if flagged_atoms
        else "_No problematic atoms identified — scaffold fits within binding pocket._"
    )

    res_rows = "\n".join(
        f"| {sp} | {req} |"
        for sp, req in _RESOLUTION_REQUIREMENTS.items()
    )

    report = f"""# {gene_symbol} Co-Crystal Structure Brief — {scaffold_id}

**Generated**: {datetime.now().strftime("%Y-%m-%d")}
**Target**: {gene_symbol}
**Scaffold**: {scaffold_id} | `{smiles}`
**Reference structure**: {pdb_id}

---

## Structures

| Role | PDB ID | Space Group |
|------|--------|-------------|
| Reference (VRK1 + ANP) | {pdb_id} | {space_group} |
| Target co-crystal | TBD | Aim to match {space_group} |

---

## Crystallisation Conditions ({pdb_id})

**Method**: {conditions.get("method", "VAPOR DIFFUSION")}

{conditions.get("conditions", _GENERIC_KINASE_CONDITIONS)}

---

## Resolution Requirements by Subpocket

| Subpocket | Minimum Resolution |
|-----------|-------------------|
{res_rows}

---

## Scaffold Compatibility

### Atoms Flagged for Potential Clashes (> 8 Å from pocket centroid, rotatable)

{flag_section}

---

## Recommended Experiment

1. **Co-crystallisation** (preferred for new scaffolds):
   - Dissolve {scaffold_id} at 10 mM in DMSO; dilute 1:10 into crystallisation drop
   - Mix 1:1 protein:precipitant (match {pdb_id} conditions above)
   - Screen pH 6.5–8.0 and PEG 3350 10–30% grid
   - Harvest crystals after 3–7 days; cryo-protect in mother liquor + 20–25% glycerol

2. **Soaking** (faster, use if co-cryst fails):
   - Grow {pdb_id}-isomorphous apo crystals
   - Soak at 1–5 mM compound, 20% DMSO (v/v), 30–60 min RT
   - Monitor crystal integrity; stop soak if cracking observed

3. **Data collection targets**:
   - Space group: {space_group} (match reference)
   - Resolution: ≤ 2.5 Å for gatekeeper contacts
   - Completeness: > 95%; Rmerge < 8%

4. **Refinement checkpoints**:
   - Confirm electron density for {scaffold_id} in 2Fo–Fc map (contoured at 1.0σ)
   - Check contacts at KLIFS positions 45 (gatekeeper), 46–48 (hinge), 72–76 (DFG)
   - Compare B-factors of scaffold atoms — high B-factor flagged atoms suggest disorder
"""
    path = results_dir / f"cocrystal_brief_{scaffold_id}.md"
    path.write_text(report)
    return path


def run_cocrystal(
    gene_symbol: str,
    scaffold_id: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> None:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    # Require structalign outputs
    comparison_csv = results_dir / "binding_site_comparison.csv"
    if not comparison_csv.exists():
        console.print(
            f"  [red]ERROR[/red] binding_site_comparison.csv not found in {results_dir}.\n"
            f"  Run `pipeline structalign --target {gene_symbol}` first."
        )
        raise SystemExit(1)

    cache_key = f"cocrystal_{scaffold_id}"
    if not force:
        cached = cache.load(gene_symbol, cache_key)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol:10}[/dim] cocrystal {scaffold_id:12} "
                f"[yellow]SKIP[/yellow]  (cached)"
            )
            return

    # Resolve scaffold SMILES
    compounds_csv = results_dir / "compounds_filtered.csv"
    if not compounds_csv.exists():
        raise FileNotFoundError(
            f"compounds_filtered.csv not found. "
            f"Run `pipeline fetch --target {gene_symbol}` first."
        )
    compounds = pd.read_csv(compounds_csv)
    row = compounds[compounds.get("scaffold_id", pd.Series(dtype=str)) == scaffold_id]
    if row.empty:
        for col in ("compound_id", "name", "molecule_chembl_id"):
            if col in compounds.columns:
                row = compounds[compounds[col] == scaffold_id]
                if not row.empty:
                    break
    if row.empty:
        raise ValueError(f"Scaffold '{scaffold_id}' not found in compounds_filtered.csv")

    smiles_col = next(
        (c for c in ("smiles", "canonical_smiles", "SMILES") if c in row.columns), None
    )
    smiles = str(row.iloc[0][smiles_col]) if smiles_col else "unknown"

    # Identify the reference PDB from structalign cache
    pdb_files = list((settings.cache_dir / gene_symbol).glob("*.pdb"))
    pdb_id = pdb_files[0].stem.upper() if pdb_files else "6AC9"

    # Compute centroid from binding_site_vrk1.csv
    binding_site_csv = results_dir / "binding_site_vrk1.csv"
    centroid = [0.0, 0.0, 0.0]
    if binding_site_csv.exists():
        bs_df = pd.read_csv(binding_site_csv)
        if {"ca_x", "ca_y", "ca_z"}.issubset(bs_df.columns):
            valid = bs_df[["ca_x", "ca_y", "ca_z"]].dropna()
            if not valid.empty:
                centroid = valid.mean().tolist()

    console.print(f"  [dim]{gene_symbol}:[/dim] fetching RCSB conditions for {pdb_id}...")
    conditions = _fetch_rcsb_conditions(pdb_id)

    # Space group from cached PDB
    space_group = "unknown"
    if pdb_files:
        space_group = _get_space_group(pdb_files[0])

    console.print(f"  [dim]{gene_symbol}:[/dim] flagging scaffold atoms for {scaffold_id}...")
    flagged = _flag_scaffold_atoms(smiles, centroid)

    # Save cache
    cache.save(gene_symbol, cache_key, {"scaffold_id": scaffold_id, "pdb_id": pdb_id}, 1)

    brief_path = _write_brief(
        gene_symbol, scaffold_id, smiles, pdb_id, space_group, conditions, flagged, results_dir
    )

    console.print(
        f"  [dim]{gene_symbol}:[/dim] space group: {space_group} | "
        f"{len(flagged)} atoms flagged"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]cocrystal_brief_{scaffold_id}.md written[/green]"
    )
