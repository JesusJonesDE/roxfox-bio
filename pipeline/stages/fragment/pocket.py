"""
Fragment pipeline — Pocket identification step.

Locates the binding pocket in the IGHMBP2 AlphaFold2 structure using
fpocket, parses the druggability descriptors, and extracts a centroid
for downstream docking box placement.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Optional

import httpx
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings

AF2_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-P38935-F1-model_v4.pdb"
AF2_PDB_FILENAME = "IGHMBP2_AF2.pdb"

# Column order in fpocket *_info.txt (whitespace-separated, after the header)
_COLUMNS = [
    "Pocket",
    "Score",
    "Drugg_score",
    "Total_SASA",
    "Polar_SASA",
    "Apolar_SASA",
    "Volume",
    "Mean_local_hydrophobic_dens",
    "Mean_acrophobicity",
    "Mean_pocketness",
    "Mean_accept_angle",
    "Mean_donor_angle",
    "Meam_water_dist",
    "Flexibility",
]

_MIN_VOLUME = 200.0
_RELAXED_VOLUME = 150.0


# ── AlphaFold2 PDB helpers ─────────────────────────────────────────────────────

def _find_pdb(gene_symbol: str, settings: Settings) -> Optional[Path]:
    """Return an existing *.pdb in the structures dir, or None."""
    struct_dir = settings.cache_dir / gene_symbol / "structures"
    if struct_dir.exists():
        pdbs = sorted(struct_dir.glob("*.pdb"))
        if pdbs:
            return pdbs[0]
    return None


def _download_af2_pdb(gene_symbol: str, settings: Settings) -> Path:
    struct_dir = settings.cache_dir / gene_symbol / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    dest = struct_dir / AF2_PDB_FILENAME
    with httpx.Client(timeout=60) as client:
        r = client.get(AF2_PDB_URL)
        r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


# ── fpocket helpers ────────────────────────────────────────────────────────────

def _run_fpocket(pdb_path: Path) -> None:
    """Run fpocket on pdb_path (cwd = pdb_path.parent)."""
    try:
        subprocess.run(
            ["fpocket", "-f", str(pdb_path)],
            check=True,
            capture_output=True,
            cwd=str(pdb_path.parent),
        )
    except FileNotFoundError:
        raise RuntimeError("fpocket not installed. Run: brew install fpocket")


def _parse_info_txt(info_path: Path) -> list[dict]:
    """Parse *_info.txt → list of row dicts keyed by _COLUMNS."""
    rows: list[dict] = []
    for line in info_path.read_text(errors="replace").splitlines():
        line = line.strip()
        # Skip blank lines, comment lines, and the header line
        if not line or line.startswith("#") or line.startswith("Pocket"):
            continue
        parts = line.split()
        if len(parts) < len(_COLUMNS):
            continue
        try:
            row = {col: parts[i] for i, col in enumerate(_COLUMNS)}
            rows.append(row)
        except (IndexError, ValueError):
            continue
    return rows


def _select_pocket(rows: list[dict], console: Console) -> dict:
    """Select the pocket with the highest Drugg_score and Volume > threshold."""
    if not rows:
        raise RuntimeError("No pockets found by fpocket")

    def _drugg(r: dict) -> float:
        try:
            return float(r["Drugg_score"])
        except (KeyError, ValueError):
            return 0.0

    def _vol(r: dict) -> float:
        try:
            return float(r["Volume"])
        except (KeyError, ValueError):
            return 0.0

    candidates = [r for r in rows if _vol(r) > _MIN_VOLUME]

    if not candidates:
        console.print(
            f"  [yellow]Warning: no pockets with volume > {_MIN_VOLUME} Å³ — "
            f"relaxing threshold to {_RELAXED_VOLUME} Å³[/yellow]"
        )
        candidates = [r for r in rows if _vol(r) > _RELAXED_VOLUME]

    if not candidates:
        raise RuntimeError("No pockets found by fpocket")

    return max(candidates, key=_drugg)


def _extract_centroid(pocket_pdb: Path) -> tuple[float, float, float]:
    """Compute the mean X/Y/Z of all ATOM/HETATM records in a pocket PDB."""
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for line in pocket_pdb.read_text(errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue
    if not xs:
        raise RuntimeError(f"No ATOM/HETATM records found in {pocket_pdb}")
    n = len(xs)
    return sum(xs) / n, sum(ys) / n, sum(zs) / n


# ── Public entry point ─────────────────────────────────────────────────────────

def run_pocket(
    gene_symbol: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> dict:
    """Identify the best druggable pocket in the AF2 structure.

    Returns a dict with keys:
        pocket_id, score, druggability_score, volume_A3,
        centroid_x, centroid_y, centroid_z, box_size_A, plddt_mean
    """
    # ── 1. Cache check ────────────────────────────────────────────────────────
    if not force:
        cached = cache.load(gene_symbol, "fragment_pocket")
        if cached is not None:
            console.print(f"  [dim]{gene_symbol}:[/dim] pocket [yellow]SKIP[/yellow] (cached)")
            return cached

    # ── 2. Locate PDB ─────────────────────────────────────────────────────────
    pdb_path = _find_pdb(gene_symbol, settings)
    if pdb_path is None:
        console.print(f"  [dim]{gene_symbol}:[/dim] AF2 PDB not found locally — downloading from EBI...")
        try:
            pdb_path = _download_af2_pdb(gene_symbol, settings)
            console.print(f"  [dim]{gene_symbol}:[/dim] downloaded {pdb_path.name}")
        except Exception as exc:
            raise RuntimeError(
                f"Failed to download AF2 PDB: {exc}. "
                f"Ensure network access or run: pipeline fetch --target {gene_symbol}"
            ) from exc

    pdb_path = _find_pdb(gene_symbol, settings)
    if pdb_path is None:
        raise RuntimeError(
            f"No PDB file found for {gene_symbol}. "
            f"Run: pipeline fetch --target {gene_symbol}"
        )

    # ── 3. Run fpocket ────────────────────────────────────────────────────────
    console.print(f"  [dim]{gene_symbol}:[/dim] running fpocket on {pdb_path.name}...")
    _run_fpocket(pdb_path)

    # ── 4. Parse output ───────────────────────────────────────────────────────
    pdb_stem = pdb_path.stem
    out_dir = pdb_path.parent / f"{pdb_stem}_out"
    if not out_dir.exists():
        raise RuntimeError(f"fpocket output directory not created: {out_dir}")

    info_path = out_dir / f"{pdb_stem}_info.txt"
    if not info_path.exists():
        raise RuntimeError(f"fpocket info file not found: {info_path}")

    rows = _parse_info_txt(info_path)
    pocket_row = _select_pocket(rows, console)
    pocket_id = int(pocket_row["Pocket"])

    # ── 5. Extract centroid ───────────────────────────────────────────────────
    pocket_pdb = out_dir / "pockets" / f"pocket{pocket_id}_atm.pdb"
    if not pocket_pdb.exists():
        raise RuntimeError(f"Pocket atom PDB not found: {pocket_pdb}")
    cx, cy, cz = _extract_centroid(pocket_pdb)

    # ── 6. Build result ───────────────────────────────────────────────────────
    result: dict = {
        "pocket_id": pocket_id,
        "score": float(pocket_row["Score"]),
        "druggability_score": float(pocket_row["Drugg_score"]),
        "volume_A3": float(pocket_row["Volume"]),
        "centroid_x": float(cx),
        "centroid_y": float(cy),
        "centroid_z": float(cz),
        "box_size_A": 20.0,
        "plddt_mean": None,
    }

    # ── 7. Write JSON ─────────────────────────────────────────────────────────
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)
    out_json = results_dir / "pocket_analysis.json"
    out_json.write_text(json.dumps(result, indent=2))
    console.print(f"  [dim]{gene_symbol}:[/dim] [green]pocket_analysis.json written[/green]")

    # ── 8. Cache ──────────────────────────────────────────────────────────────
    cache.save(gene_symbol, "fragment_pocket", result, 1)

    return result
