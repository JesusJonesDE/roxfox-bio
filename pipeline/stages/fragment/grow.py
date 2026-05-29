"""
Fragment pipeline — Fragment growing step.

Grows cluster representatives into drug-like molecules using BRICS
combinatorial assembly and SMARTS reaction-based transformations.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, Descriptors, rdMolDescriptors

from pipeline.cache import CacheManager
from pipeline.config import Settings


# ── SA scorer (optional) ───────────────────────────────────────────────────────

_SASCORER_PATH = (
    Path(sys.executable).parent.parent
    / "lib"
    / f"python{sys.version_info.major}.{sys.version_info.minor}"
    / "site-packages/rdkit/Contrib/SA_Score"
)

try:
    sys.path.insert(0, str(_SASCORER_PATH))
    import sascorer  # type: ignore
    HAS_SASCORER = True
except ImportError:
    HAS_SASCORER = False


# ── SMARTS reaction definitions ────────────────────────────────────────────────

REACTIONS = [
    ("amide",       "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"),
    ("n_methyl",    "[NH:1]>>[N:1]C"),
    ("n_ethyl",     "[NH:1]>>[N:1]CC"),
    ("o_methyl",    "[OH:1]>>[O:1]C"),
    ("f_subst",     "[c:1][H]>>[c:1]F"),
    ("cf3_subst",   "[c:1][H]>>[c:1]C(F)(F)F"),
    ("methyl_subst","[c:1][H]>>[c:1]C"),
    ("ethyl_subst", "[c:1][H]>>[c:1]CC"),
    ("amino_subst", "[c:1][H]>>[c:1]N"),
    ("hydroxy_subst","[c:1][H]>>[c:1]O"),
    ("sulfonamide", "[NH2:1]>>[N:1]S(=O)(=O)C"),
    ("acetamide",   "[NH2:1]>>[N:1]C(=O)C"),
    ("urea",        "[NH2:1]>>[N:1]C(=O)N"),
    ("methylamine", "[c:1][H]>>[c:1]CN"),
    ("morpholine",  "[NH:1]>>[N:1]C1CCOCC1"),
]


# ── Private helpers ────────────────────────────────────────────────────────────

def _load_building_blocks(settings: Settings) -> list[Chem.Mol]:
    """Read building_blocks.smi and return valid RDKit Mol objects."""
    bb_path = (
        settings.cache_dir
        / "shared"
        / "fragment_library"
        / "building_blocks.smi"
    )
    if not bb_path.exists():
        return []

    mols: list[Chem.Mol] = []
    with open(bb_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            smiles = parts[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mols.append(mol)
    return mols


def _grow_brics(
    fragment_mol: Chem.Mol,
    building_blocks: list[Chem.Mol],
) -> list[Chem.Mol]:
    """Grow a fragment using BRICS assembly with each building block."""
    products: list[Chem.Mol] = []
    for bb in building_blocks:
        try:
            new_mols = list(BRICS.BRICSBuild([fragment_mol, bb]))
            products.extend(new_mols)
        except Exception:
            continue

    # Deduplicate by canonical SMILES
    seen: set[str] = set()
    result: list[Chem.Mol] = []
    for m in products:
        try:
            s = Chem.MolToSmiles(m)
            if s not in seen:
                seen.add(s)
                result.append(m)
        except Exception:
            continue
    return result


def _grow_smarts(fragment_mol: Chem.Mol) -> list[Chem.Mol]:
    """Apply 15 SMARTS transformations to grow a fragment."""
    from rdkit.Chem import AddHs, RemoveHs

    # Add explicit Hs so [c:1][H] patterns can match aromatic C-H bonds
    frag_with_hs = AddHs(fragment_mol)

    products: list[Chem.Mol] = []
    for name, smarts in REACTIONS:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            # Try with explicit Hs first (needed for aromatic C-H substitution)
            try:
                products_tuple = rxn.RunReactants((frag_with_hs,))
            except Exception:
                products_tuple = rxn.RunReactants((fragment_mol,))
            for prod_tuple in products_tuple[:3]:  # limit products per reaction
                mol = prod_tuple[0]
                try:
                    # Remove explicit Hs and sanitize
                    mol = RemoveHs(mol)
                    Chem.SanitizeMol(mol)
                    products.append(mol)
                except Exception:
                    continue
        except Exception:
            continue

    # Deduplicate by canonical SMILES
    seen: set[str] = set()
    result: list[Chem.Mol] = []
    for m in products:
        try:
            s = Chem.MolToSmiles(m)
            if s not in seen:
                seen.add(s)
                result.append(m)
        except Exception:
            continue
    return result


def _filter_grown(mols: list[Chem.Mol]) -> list[tuple[Chem.Mol, dict]]:
    """Filter grown molecules by drug-likeness criteria.

    Keeps molecules with:
    - MW 300–450 Da
    - Ro5 violations == 0 (Lipinski)
    - RotatableBonds <= 8
    - SA score < 4 (if sascorer available)

    Returns list of (mol, props_dict).
    """
    if not HAS_SASCORER:
        warnings.warn(
            "sascorer not available — skipping SA score filter. "
            "Install with: pip install sa-score",
            stacklevel=2,
        )

    result: list[tuple[Chem.Mol, dict]] = []
    for mol in mols:
        try:
            mw = Descriptors.ExactMolWt(mol)
            if not (300.0 <= mw <= 450.0):
                continue

            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
            ro5_violations = int(
                (mw > 500) + (logp > 5) + (hbd > 5) + (hba > 10)
            )

            if ro5_violations != 0:
                continue

            if rotb > 8:
                continue

            sa_score: Optional[float] = None
            if HAS_SASCORER:
                try:
                    sa_score = sascorer.calculateScore(mol)
                    if sa_score >= 4.0:
                        continue
                except Exception:
                    pass

            props = {
                "molecular_weight": round(mw, 3),
                "logp": round(logp, 3),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rotb,
                "ro5_violations": ro5_violations,
                "passes_ro5": True,
                "sa_score": sa_score if sa_score is not None else float("nan"),
            }
            result.append((mol, props))
        except Exception:
            continue
    return result


# ── Public entry point ─────────────────────────────────────────────────────────

def run_grow(
    gene_symbol: str,
    clusters_df: pd.DataFrame,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> pd.DataFrame:
    """Grow cluster representatives into drug-like candidates.

    Returns a DataFrame with columns:
        candidate_id, smiles, parent_fragment_id, molecular_weight, logp,
        hbd, hba, rotatable_bonds, ro5_violations, passes_ro5, sa_score,
        grow_method
    """
    # 1. Cache check
    cache_key = "fragment_grow"
    if not force:
        cached = cache.load(gene_symbol, cache_key)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] grow [yellow]SKIP[/yellow] (cached)"
            )
            return pd.DataFrame(cached)

    # 2. Load building blocks
    building_blocks = _load_building_blocks(settings)
    console.print(
        f"  [dim]{gene_symbol}:[/dim] loaded {len(building_blocks)} building blocks"
    )

    # 3. Get cluster representatives
    if "is_representative" not in clusters_df.columns or len(clusters_df) == 0:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [yellow]no cluster representatives — "
            f"returning empty DataFrame[/yellow]"
        )
        return pd.DataFrame(
            columns=[
                "candidate_id", "smiles", "parent_fragment_id",
                "molecular_weight", "logp", "hbd", "hba",
                "rotatable_bonds", "ro5_violations", "passes_ro5",
                "sa_score", "grow_method",
            ]
        )

    reps = clusters_df[clusters_df["is_representative"] == True]

    if len(reps) == 0:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [yellow]no representatives in clusters_df — "
            f"returning empty DataFrame[/yellow]"
        )
        return pd.DataFrame(
            columns=[
                "candidate_id", "smiles", "parent_fragment_id",
                "molecular_weight", "logp", "hbd", "hba",
                "rotatable_bonds", "ro5_violations", "passes_ro5",
                "sa_score", "grow_method",
            ]
        )

    console.print(
        f"  [dim]{gene_symbol}:[/dim] growing {len(reps)} cluster representatives..."
    )

    # 4. Grow each representative
    all_candidates: list[dict] = []

    for _, rep_row in reps.iterrows():
        smiles = rep_row.get("smiles", "")
        fragment_id = rep_row.get("fragment_id", "unknown")

        if not smiles:
            continue

        frag_mol = Chem.MolFromSmiles(str(smiles))
        if frag_mol is None:
            continue

        # BRICS growing
        brics_mols = _grow_brics(frag_mol, building_blocks)
        # SMARTS growing
        smarts_mols = _grow_smarts(frag_mol)

        # Combine and deduplicate
        combined_raw = brics_mols + smarts_mols
        seen_smiles: set[str] = set()
        combined: list[tuple[Chem.Mol, str]] = []  # (mol, method)
        for mol in brics_mols:
            try:
                s = Chem.MolToSmiles(mol)
                if s not in seen_smiles:
                    seen_smiles.add(s)
                    combined.append((mol, "brics"))
            except Exception:
                continue
        for mol in smarts_mols:
            try:
                s = Chem.MolToSmiles(mol)
                if s not in seen_smiles:
                    seen_smiles.add(s)
                    combined.append((mol, "smarts"))
            except Exception:
                continue

        # Filter
        mol_list = [m for m, _ in combined]
        method_map: dict[str, str] = {}
        for mol, method in combined:
            try:
                s = Chem.MolToSmiles(mol)
                method_map[s] = method
            except Exception:
                pass

        filtered = _filter_grown(mol_list)

        # Sort by MW descending (heaviest = most elaborated), keep top 3
        filtered.sort(key=lambda x: x[1]["molecular_weight"], reverse=True)
        top3 = filtered[:3]

        for mol, props in top3:
            try:
                smi = Chem.MolToSmiles(mol)
                method = method_map.get(smi, "unknown")
                all_candidates.append(
                    {
                        "smiles": smi,
                        "parent_fragment_id": fragment_id,
                        "grow_method": method,
                        **props,
                    }
                )
            except Exception:
                continue

    # 5. Warn if fewer than 20 candidates
    if len(all_candidates) < 20:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [yellow]Warning: only {len(all_candidates)} "
            f"candidates grown (target: 20) — using all available[/yellow]"
        )

    # 6. Assign candidate IDs
    rows: list[dict] = []
    for n, cand in enumerate(all_candidates, start=1):
        rows.append(
            {
                "candidate_id": f"IGHMBP2-SCF-{n:03d}",
                "smiles": cand["smiles"],
                "parent_fragment_id": cand["parent_fragment_id"],
                "molecular_weight": cand["molecular_weight"],
                "logp": cand["logp"],
                "hbd": cand["hbd"],
                "hba": cand["hba"],
                "rotatable_bonds": cand["rotatable_bonds"],
                "ro5_violations": cand["ro5_violations"],
                "passes_ro5": cand["passes_ro5"],
                "sa_score": cand["sa_score"],
                "grow_method": cand["grow_method"],
            }
        )

    candidates_df = pd.DataFrame(
        rows,
        columns=[
            "candidate_id", "smiles", "parent_fragment_id",
            "molecular_weight", "logp", "hbd", "hba",
            "rotatable_bonds", "ro5_violations", "passes_ro5",
            "sa_score", "grow_method",
        ],
    )

    # 7. Write grown_candidates.csv
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)
    candidates_df.to_csv(results_dir / "grown_candidates.csv", index=False)
    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]grow complete — "
        f"{len(candidates_df)} candidates written to grown_candidates.csv[/green]"
    )

    # 8. Cache and return
    cache.save(
        gene_symbol, cache_key, candidates_df.to_dict(orient="records"), len(candidates_df)
    )
    return candidates_df
