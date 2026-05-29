"""
Fragment pipeline — Fragment library preparation step.

Downloads Ro3-compliant fragments from ZINC22, applies RDKit filters,
strips salts, deduplicates by Murcko scaffold, and writes a SMILES file
ready for docking. Falls back to a bundled library if downloads fail.
"""
from __future__ import annotations

import random
import warnings
from pathlib import Path
from typing import Optional

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings

# ZINC22 tranche URLs (2D, AA prefix)
_ZINC_TRANCHES = [
    "https://files.docking.org/2D/AA/AAAA.smi",
    "https://files.docking.org/2D/AA/AAAB.smi",
    "https://files.docking.org/2D/AA/AAAC.smi",
    "https://files.docking.org/2D/AA/AAAD.smi",
    "https://files.docking.org/2D/AA/AAAE.smi",
]

# Ro3 thresholds (Rule of Three, Congreve et al., 2003)
_RO3_MW = 250.0
_RO3_HBD = 3
_RO3_HBA = 3
_RO3_LOGP = 3.0
_RO3_ROTB = 3


# ── ZINC download ──────────────────────────────────────────────────────────────

@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _fetch_tranche(url: str) -> list[str]:
    """Download one ZINC tranche file and return all SMILES strings."""
    with httpx.Client(timeout=60) as client:
        r = client.get(url)
        r.raise_for_status()
    smiles_list: list[str] = []
    for line in r.text.splitlines():
        parts = line.strip().split()
        if parts:
            smiles_list.append(parts[0])
    return smiles_list


def _download_zinc_smiles(console: Console) -> Optional[list[str]]:
    """Try to fetch all 5 ZINC22 tranches. Returns None on any failure."""
    all_smiles: list[str] = []
    for url in _ZINC_TRANCHES:
        try:
            chunk = _fetch_tranche(url)
            all_smiles.extend(chunk)
            console.print(f"  [dim]library:[/dim] fetched {len(chunk)} SMILES from {url.split('/')[-1]}")
        except (httpx.HTTPStatusError, httpx.TransportError, Exception) as exc:
            console.print(f"  [yellow]ZINC download failed ({url.split('/')[-1]}): {exc} — falling back to bundled library[/yellow]")
            return None
    return all_smiles


# ── RDKit filters ──────────────────────────────────────────────────────────────

def _passes_ro3(mol) -> bool:  # type: ignore[no-untyped-def]
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mw = Descriptors.ExactMolWt(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    logp = Descriptors.MolLogP(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return mw <= _RO3_MW and hbd <= _RO3_HBD and hba <= _RO3_HBA and logp <= _RO3_LOGP and rotb <= _RO3_ROTB


def _strip_salts(mol):  # type: ignore[no-untyped-def]
    from rdkit.Chem.SaltRemover import SaltRemover
    remover = SaltRemover()
    return remover.StripMol(mol)


def _murcko_scaffold_smi(mol) -> str:  # type: ignore[no-untyped-def]
    from rdkit.Chem import MolToSmiles
    from rdkit.Chem.Scaffolds import MurckoScaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return MolToSmiles(scaffold)


def _apply_filters(smiles_list: list[str], console: Console) -> list[str]:
    """Apply Ro3 filter, salt strip, and Murcko scaffold deduplication."""
    from rdkit import Chem
    from rdkit.Chem import MolToSmiles

    seen_scaffolds: set[str] = set()
    filtered: list[str] = []
    n_invalid = 0
    n_ro3_fail = 0
    n_scaffold_dup = 0

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_invalid += 1
            continue

        mol = _strip_salts(mol)
        if mol is None:
            n_invalid += 1
            continue

        if not _passes_ro3(mol):
            n_ro3_fail += 1
            continue

        try:
            scaffold_smi = _murcko_scaffold_smi(mol)
        except Exception:
            # If scaffold computation fails, use the canonical SMILES as key
            scaffold_smi = MolToSmiles(mol)

        if scaffold_smi in seen_scaffolds:
            n_scaffold_dup += 1
            continue

        seen_scaffolds.add(scaffold_smi)
        filtered.append(MolToSmiles(mol))

    console.print(
        f"  [dim]library:[/dim] {len(filtered)} passed "
        f"(invalid={n_invalid}, Ro3-fail={n_ro3_fail}, scaffold-dup={n_scaffold_dup})"
    )
    return filtered


# ── Fallback library ───────────────────────────────────────────────────────────

def _load_fallback(settings: Settings) -> list[str]:
    fallback_path = settings.cache_dir / "shared" / "fragment_library" / "fragments_fallback.smi"
    smiles_list: list[str] = []
    if fallback_path.exists():
        for line in fallback_path.read_text().splitlines():
            parts = line.strip().split()
            if parts:
                smiles_list.append(parts[0])
    return smiles_list


# ── Writer ─────────────────────────────────────────────────────────────────────

def _write_library(lib_path: Path, smiles_list: list[str], is_fallback: bool) -> int:
    lib_path.parent.mkdir(parents=True, exist_ok=True)
    prefix = "FALLBACK" if is_fallback else "ZINC"
    fmt = "{prefix}-{i:05d}" if is_fallback else "{prefix}-{i:07d}"
    lines: list[str] = []
    for i, smi in enumerate(smiles_list, 1):
        frag_id = f"{prefix}-{i:05d}" if is_fallback else f"{prefix}-{i:07d}"
        lines.append(f"{smi}\t{frag_id}")
    lib_path.write_text("\n".join(lines) + "\n")
    return len(lines)


# ── Public entry point ─────────────────────────────────────────────────────────

def run_library(
    library_size: int,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> Path:
    """Download, filter, and write the fragment library.

    Returns the path to the written SMILES file.
    """
    lib_path = settings.cache_dir / "shared" / "fragment_library" / "fragments_ro3.smi"

    # ── 1. Cache check ────────────────────────────────────────────────────────
    if not force:
        cached = cache.load("shared", "fragment_library")
        if cached is not None and lib_path.exists():
            console.print(f"  [dim]library:[/dim] [yellow]SKIP[/yellow] (cached)")
            return lib_path

    # ── 2. Download from ZINC22 ───────────────────────────────────────────────
    is_fallback = False
    raw_smiles = _download_zinc_smiles(console)

    if raw_smiles is None:
        # Fall back to bundled library
        console.print(f"  [yellow]Warning: using bundled fallback library (ZINC download failed)[/yellow]")
        raw_smiles = _load_fallback(settings)
        is_fallback = True
        if not raw_smiles:
            raise RuntimeError(
                f"Fallback library empty or missing: "
                f"{settings.cache_dir / 'shared' / 'fragment_library' / 'fragments_fallback.smi'}"
            )

    # ── 3. Apply Ro3 filter (only for ZINC data; fallback is pre-filtered) ───
    if not is_fallback:
        filtered = _apply_filters(raw_smiles, console)
    else:
        filtered = raw_smiles

    if not filtered:
        raise RuntimeError("No fragments passed Ro3 filter")

    # ── 4. Sample to library_size ─────────────────────────────────────────────
    random.seed(42)
    if len(filtered) > library_size:
        filtered = random.sample(filtered, library_size)

    # ── 5. Write library ──────────────────────────────────────────────────────
    count = _write_library(lib_path, filtered, is_fallback)
    console.print(
        f"  [dim]library:[/dim] [green]{count} fragments written to {lib_path.name}[/green]"
        + (" [yellow](fallback)[/yellow]" if is_fallback else "")
    )

    # ── 6. Cache ──────────────────────────────────────────────────────────────
    cache.save("shared", "fragment_library", {"path": str(lib_path), "size": count}, count)

    return lib_path
