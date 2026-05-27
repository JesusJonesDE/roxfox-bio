"""
Druggability pocket analysis using fpocket.

Downloads X-ray PDB structures, runs fpocket, parses pocket descriptors,
and writes pocket_analysis.csv + druggability_report.md.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import httpx
import pandas as pd
from rich.console import Console

from pipeline.config import Settings

FPOCKET_BIN = Path("/tmp/fpocket_build/bin/fpocket")
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Druggability score thresholds (Schmidtke & Barril, 2010)
DRUGGABLE_SCORE = 0.5   # drug_score ≥ 0.5 → likely druggable
BORDERLINE_SCORE = 0.2  # 0.2–0.5 → borderline

# Min pocket volume (Å³) to consider as relevant
MIN_VOLUME = 200


def _download_pdb(pdb_id: str, dest: Path) -> None:
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    with httpx.Client(timeout=30) as client:
        r = client.get(url)
        r.raise_for_status()
    dest.write_bytes(r.content)


def _run_fpocket(pdb_path: Path, work_dir: Path) -> Optional[Path]:
    result = subprocess.run(
        [str(FPOCKET_BIN), "-f", str(pdb_path)],
        cwd=str(work_dir),
        capture_output=True,
        text=True,
    )
    out_dir = work_dir / f"{pdb_path.stem}_out"
    if result.returncode != 0 or not out_dir.exists():
        return None
    return out_dir


def _parse_info_file(info_path: Path) -> list[dict]:
    """Parse fpocket *_info.txt → list of pocket dicts."""
    text = info_path.read_text(errors="replace")
    pockets = []
    # Each pocket block starts with "Pocket N :"
    blocks = re.split(r"Pocket\s+(\d+)\s+:", text)
    # blocks[0] = preamble, then pairs of (pocket_num, block_text)
    for i in range(1, len(blocks), 2):
        pocket_num = int(blocks[i])
        block = blocks[i + 1]

        def _val(pattern: str) -> Optional[float]:
            m = re.search(pattern, block)
            return float(m.group(1)) if m else None

        pockets.append({
            "pocket_id": pocket_num,
            "score": _val(r"Druggability Score\s*:\s*([\d.]+)"),
            "volume_A3": _val(r"Volume\s*:\s*([\d.]+)"),
            "hydrophobicity_score": _val(r"Hydrophobicity score\s*:\s*([\d.]+)"),
            "polarity_score": _val(r"Polarity score\s*:\s*([\d.]+)"),
            "charge_score": _val(r"Charge score\s*:\s*([\d.-]+)"),
            "flexibility": _val(r"Flexibility\s*:\s*([\d.]+)"),
            "n_alpha_spheres": _val(r"Number of Alpha Spheres\s*:\s*([\d.]+)"),
            "mean_local_hydrophobic_density": _val(r"Mean local hydrophobic density\s*:\s*([\d.]+)"),
        })
    return pockets


def _druggability_label(score: Optional[float]) -> str:
    if score is None:
        return "unknown"
    if score >= DRUGGABLE_SCORE:
        return "druggable"
    if score >= BORDERLINE_SCORE:
        return "borderline"
    return "undruggable"


def run_pocket_analysis(gene: str, settings: Settings, cache: object, console: Console) -> Path:
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)

    if not FPOCKET_BIN.exists():
        raise RuntimeError(f"fpocket binary not found at {FPOCKET_BIN}. Build it first.")

    # Load PDB structure list from cache
    pdb_data = cache.load(gene, "pdb")
    if not pdb_data:
        raise FileNotFoundError(f"No PDB cache for {gene} — run `pipeline fetch --target {gene}` first")

    all_structures = pdb_data.get("structures", [])

    # Only X-ray structures (NMR ensembles confuse fpocket geometry)
    xray = [s for s in all_structures if s.get("method") == "X-RAY DIFFRACTION"]
    if not xray:
        console.print(f"  [yellow]{gene}: no X-ray structures — falling back to all structures[/yellow]")
        xray = all_structures

    console.print(f"  {gene}: {len(xray)} X-ray structure(s) → running fpocket")

    all_rows = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for struct in xray:
            pdb_id = struct["pdb_id"]
            pdb_file = tmp_path / f"{pdb_id}.pdb"
            try:
                _download_pdb(pdb_id, pdb_file)
            except Exception as exc:
                console.print(f"    [red]{pdb_id}: download failed — {exc}[/red]")
                continue

            out_dir = _run_fpocket(pdb_file, tmp_path)
            if out_dir is None:
                console.print(f"    [red]{pdb_id}: fpocket failed[/red]")
                continue

            info_files = list(out_dir.glob("*_info.txt"))
            if not info_files:
                console.print(f"    [yellow]{pdb_id}: no info file in fpocket output[/yellow]")
                continue

            pockets = _parse_info_file(info_files[0])
            for p in pockets:
                p["pdb_id"] = pdb_id
                p["resolution_A"] = struct.get("resolution_angstrom")
                p["has_known_ligand"] = struct.get("has_ligand", False)
                p["druggability"] = _druggability_label(p["score"])
            all_rows.extend(pockets)
            console.print(f"    {pdb_id}: {len(pockets)} pockets detected")

    if not all_rows:
        raise RuntimeError(f"{gene}: fpocket returned no pockets for any structure")

    df = pd.DataFrame(all_rows)
    df = df[df["volume_A3"] >= MIN_VOLUME].copy()  # drop trivially small cavities

    # Sort: best scoring first, then by volume
    df = df.sort_values(["score", "volume_A3"], ascending=[False, False]).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)

    out_csv = results_dir / "pocket_analysis.csv"
    df.to_csv(out_csv, index=False)

    _write_druggability_report(gene, df, results_dir)

    druggable = df[df["druggability"] == "druggable"]
    borderline = df[df["druggability"] == "borderline"]
    console.print(
        f"  {gene}: {len(df)} pockets (≥{MIN_VOLUME}Å³)  |  "
        f"[green]{len(druggable)} druggable[/green]  |  "
        f"[yellow]{len(borderline)} borderline[/yellow]"
    )
    return out_csv


def _write_druggability_report(gene: str, df: pd.DataFrame, results_dir: Path) -> None:
    druggable = df[df["druggability"] == "druggable"]
    borderline = df[df["druggability"] == "borderline"]

    top = df.iloc[0] if not df.empty else None

    lines = [
        f"# {gene} — Druggability Report",
        "",
        "## Summary",
        "",
        f"- Structures analysed: {df['pdb_id'].nunique()} X-ray structures",
        f"- Total pockets (volume ≥ 200 Å³): {len(df)}",
        f"- Druggable (score ≥ 0.5): **{len(druggable)}**",
        f"- Borderline (score 0.2–0.5): **{len(borderline)}**",
        f"- Undruggable: {len(df) - len(druggable) - len(borderline)}",
        "",
        "## Top Pocket",
        "",
    ]

    if top is not None:
        lines += [
            f"| Property | Value |",
            f"|---|---|",
            f"| PDB structure | {top['pdb_id']} |",
            f"| Pocket # | {int(top['pocket_id'])} |",
            f"| Druggability score | {top['score']:.3f} ({top['druggability']}) |",
            f"| Volume | {top['volume_A3']:.0f} Å³ |",
            f"| Hydrophobicity | {top['hydrophobicity_score']:.2f} |",
            f"| Alpha spheres | {int(top['n_alpha_spheres']) if top['n_alpha_spheres'] else 'N/A'} |",
            f"| Known ligand in structure | {'Yes' if top['has_known_ligand'] else 'No'} |",
            "",
        ]

    lines += [
        "## All Pockets (ranked by druggability score)",
        "",
        "| Rank | PDB | Pocket | Score | Druggability | Volume (Å³) | Hydrophobicity |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, row in df.iterrows():
        score = f"{row['score']:.3f}" if row["score"] is not None else "N/A"
        vol = f"{row['volume_A3']:.0f}" if row["volume_A3"] is not None else "N/A"
        hydro = f"{row['hydrophobicity_score']:.2f}" if row["hydrophobicity_score"] is not None else "N/A"
        lines.append(
            f"| {int(row['rank'])} | {row['pdb_id']} | {int(row['pocket_id'])} "
            f"| {score} | {row['druggability']} | {vol} | {hydro} |"
        )

    lines += [
        "",
        "## Druggability Assessment",
        "",
        "Scores use the Schmidtke & Barril (2010) scale: ≥ 0.5 = druggable, 0.2–0.5 = borderline, < 0.2 = undruggable.",
        "",
    ]

    if not druggable.empty:
        lines.append("**Druggable pockets identified** — this target has at least one geometrically suitable binding site for small-molecule drug discovery.")
    elif not borderline.empty:
        lines.append("**Borderline druggability** — binding sites exist but may require optimized compounds (e.g., higher MW, allosteric approach, or PROTAC strategy).")
    else:
        lines.append("**No druggable pockets detected** — consider alternative modalities (PROTACs, RNA targeting, protein-protein interaction disruption).")

    out_md = results_dir / "druggability_report.md"
    out_md.write_text("\n".join(lines) + "\n")
