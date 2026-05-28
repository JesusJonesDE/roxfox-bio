"""MM-GBSA validation gate.

Rescores the top docking pose using gmx_MMPBSA and returns ΔG (kcal/mol).
PASS threshold: ΔG ≤ -7.0 kcal/mol.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings
from pipeline.models import GateResult, GateStatus
from pipeline.stages.validate.validate import (
    _load_docking_pdbqt,
    _load_cached_gate_result,
    _cache_gate_result,
    _write_gate_report,
)

_GATE_NAME = "mmgbsa"
_PASS_THRESHOLD = -7.0  # kcal/mol

_MMPBSA_IN_TEMPLATE = """\
&general
  startframe=1, endframe=1, interval=1,
  verbose=2,
/
&gb
  igb=5, saltcon=0.150,
/
"""


# ── PDBQT → PDB conversion ────────────────────────────────────────────────────

def _pdbqt_to_pdb(pdbqt_path: Path, out_pdb_path: Path) -> None:
    """Convert top-pose PDBQT to PDB format.

    Reads only ATOM/HETATM lines from the first MODEL block and strips
    the last two columns (partial charge + atom type) that PDBQT appends
    beyond the standard 80-character PDB record, producing a clean PDB file.
    """
    lines_out: list[str] = []
    in_model = False

    with pdbqt_path.open() as fh:
        for line in fh:
            record = line[:6].strip()
            if record == "MODEL":
                in_model = True
                continue
            if record == "ENDMDL":
                # Only process the first MODEL block
                break
            if not in_model:
                # Lines before any MODEL record — still include ATOM/HETATM
                pass
            if record in ("ATOM", "HETATM"):
                # Standard PDB columns are 1-80; strip PDBQT extra columns
                # PDBQT appends "  <charge>  <atom_type>" after column 68
                pdb_line = line[:68].rstrip() + "\n"
                lines_out.append(pdb_line)

    out_pdb_path.write_text("".join(lines_out))


# ── gmx_MMPBSA runner ─────────────────────────────────────────────────────────

def _run_gmx_mmpbsa(
    receptor_pdb: Path,
    ligand_pdb: Path,
    work_dir: Path,
) -> float:
    """Run gmx_MMPBSA and return ΔG_bind in kcal/mol.

    Raises RuntimeError if gmx_MMPBSA is not installed or if the calculation
    fails (subprocess error or output file not found/parseable).
    """
    if shutil.which("gmx_MMPBSA") is None:
        raise RuntimeError(
            "gmx_MMPBSA is not installed or not on PATH.\n"
            "Install via: conda install -c conda-forge gmx_mmpbsa gromacs ambertools\n"
            "Then verify with: gmx_MMPBSA --version"
        )

    work_dir.mkdir(parents=True, exist_ok=True)

    mmpbsa_in = work_dir / "mmpbsa.in"
    mmpbsa_in.write_text(_MMPBSA_IN_TEMPLATE)

    results_dat = work_dir / "FINAL_RESULTS_MMGBSA.dat"
    results_csv = work_dir / "FINAL_RESULTS_MMGBSA.csv"

    cmd = [
        "gmx_MMPBSA",
        "-O",
        "-i", str(mmpbsa_in),
        "-cs", str(receptor_pdb),
        "-ci", str(ligand_pdb),
        "-cp", str(receptor_pdb),
        "-o", str(results_dat),
        "-eo", str(results_csv),
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        stderr_snippet = (exc.stderr or "")[-500:]
        raise RuntimeError(
            f"MM-GBSA calculation failed: gmx_MMPBSA exited with code {exc.returncode}.\n"
            f"stderr (last 500 chars): {stderr_snippet}"
        ) from exc

    if not results_dat.exists():
        raise RuntimeError(
            f"MM-GBSA calculation failed: FINAL_RESULTS_MMGBSA.dat not found in {work_dir}"
        )

    # Parse ΔG_bind from the results file.
    # The relevant line looks like:
    #   TOTAL                   -8.5234  0.1234  -8.6543  -8.3921
    delta_g: float | None = None
    for line in results_dat.read_text().splitlines():
        if "TOTAL" in line:
            parts = line.split()
            # Find the numeric value after "TOTAL"
            for i, part in enumerate(parts):
                if part == "TOTAL" and i + 1 < len(parts):
                    try:
                        delta_g = float(parts[i + 1])
                        break
                    except ValueError:
                        continue
            if delta_g is not None:
                break

    if delta_g is None:
        raise RuntimeError(
            f"MM-GBSA calculation failed: could not parse TOTAL ΔG from {results_dat}"
        )

    return delta_g


# ── Gate orchestrator ─────────────────────────────────────────────────────────

def run_mmgbsa_gate(
    gene_symbol: str,
    scaffold_id: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> GateResult:
    """Orchestrate the MM-GBSA gate.

    1. Return cached result unless force=True.
    2. Convert docking PDBQT → ligand PDB.
    3. Locate receptor PDB from structures cache.
    4. Run gmx_MMPBSA.
    5. Return GateResult (PASS if ΔG ≤ -7.0, FAIL otherwise, ERROR on exception).
    """
    # 1. Cache check
    if not force:
        cached = _load_cached_gate_result(gene_symbol, scaffold_id, _GATE_NAME, cache)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] MM-GBSA result loaded from cache "
                f"(ΔG = {cached.score:.2f} kcal/mol)"
            )
            return cached

    # 2. Locate docking PDBQT
    try:
        pdbqt_path = _load_docking_pdbqt(gene_symbol, scaffold_id, settings)
    except FileNotFoundError as exc:
        result = GateResult(
            gate_name=_GATE_NAME,
            status=GateStatus.ERROR,
            score=float("nan"),
            reason=str(exc),
            details={},
        )
        _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)
        _write_gate_report(gene_symbol, scaffold_id, _GATE_NAME, result, settings)
        return result

    # 3. Locate receptor PDB
    structures_dir = settings.cache_dir / gene_symbol / "structures"
    pdb_files = list(structures_dir.glob("*.pdb"))
    if not pdb_files:
        result = GateResult(
            gate_name=_GATE_NAME,
            status=GateStatus.ERROR,
            score=float("nan"),
            reason=(
                f"No receptor PDB found in {structures_dir}. "
                f"Run `pipeline fetch --target {gene_symbol}` first."
            ),
            details={},
        )
        _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)
        _write_gate_report(gene_symbol, scaffold_id, _GATE_NAME, result, settings)
        return result

    receptor_pdb = pdb_files[0]

    # 4. Create work directory and convert PDBQT → ligand PDB
    work_dir = settings.cache_dir / gene_symbol / f"mmgbsa_{scaffold_id}"
    work_dir.mkdir(parents=True, exist_ok=True)

    ligand_pdb = work_dir / f"ligand_{scaffold_id}.pdb"
    _pdbqt_to_pdb(pdbqt_path, ligand_pdb)

    # 5. Run gmx_MMPBSA
    try:
        delta_g = _run_gmx_mmpbsa(receptor_pdb, ligand_pdb, work_dir)
    except RuntimeError as exc:
        result = GateResult(
            gate_name=_GATE_NAME,
            status=GateStatus.ERROR,
            score=float("nan"),
            reason=str(exc),
            details={},
        )
        _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)
        _write_gate_report(gene_symbol, scaffold_id, _GATE_NAME, result, settings)
        return result

    # 6. Build GateResult
    passed = delta_g <= _PASS_THRESHOLD
    status = GateStatus.PASS if passed else GateStatus.FAIL
    reason = (
        f"ΔG = {delta_g:.2f} kcal/mol ({'≤' if passed else '>'} threshold {_PASS_THRESHOLD:.1f} kcal/mol)"
    )
    details = {
        "delta_G_kcal_mol": delta_g,
        "threshold_kcal_mol": _PASS_THRESHOLD,
        "receptor_pdb": str(receptor_pdb),
        "forcefield": "GB (igb=5)",
    }

    result = GateResult(
        gate_name=_GATE_NAME,
        status=status,
        score=delta_g,
        reason=reason,
        details=details,
    )

    # 7. Write report and cache
    report_path = _write_gate_report(
        gene_symbol,
        scaffold_id,
        _GATE_NAME,
        result,
        settings,
        extra_sections=(
            f"## MM-GBSA Configuration\n\n"
            f"- Force field: GB model, igb=5\n"
            f"- Salt concentration: 0.150 M\n"
            f"- Receptor: `{receptor_pdb.name}`\n"
            f"- Ligand: converted from `{pdbqt_path.name}`\n"
        ),
    )
    result.report_path = report_path
    _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)

    return result
