"""MD pose-stability validation gate.

Prepares a protein–ligand system locally with OpenMM/PDBFixer, submits a
20 ns simulation to RunPod (A100 community cloud), polls until complete,
then analyses the ligand RMSD trajectory.

PASS threshold: mean RMSD over the final 10 ns ≤ 3.0 Å
"""
from __future__ import annotations

import base64
import os
import time
from pathlib import Path

import pandas as pd

from pipeline.cache import CacheManager
from pipeline.config import Settings
from pipeline.models import GateResult, GateStatus
from pipeline.stages.validate.validate import (
    _cache_gate_result,
    _load_cached_gate_result,
    _load_docking_pdbqt,
    _write_gate_report,
)

_GATE_NAME = "md"
_RMSD_PASS_THRESHOLD_A = 3.0   # Å
_DURATION_NS = 20.0
_A100_RATE_USD_PER_HR = 1.20
_THROUGHPUT_ATOM_NS_PER_HR = 500_000  # rough A100 + 4 fs timestep estimate
_RUNPOD_DEFAULT_ENDPOINT = "md-simulation-v1"


# ── PDBQT → PDB strip ─────────────────────────────────────────────────────────

def _strip_pdbqt_to_pdb(pdbqt_path: Path, out_pdb: Path) -> None:
    """Convert the first MODEL of a PDBQT file to a plain PDB.

    PDBQT appends partial-charge and atom-type columns beyond the standard
    80-character PDB record (columns 69–80+). Strips those extra columns from
    every ATOM/HETATM line in the first MODEL block and writes a clean PDB.
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
                break  # first MODEL only
            if record in ("ATOM", "HETATM"):
                # Columns 1-68 are standard PDB; strip everything after col 68
                pdb_line = line[:68].rstrip() + "\n"
                lines_out.append(pdb_line)

    out_pdb.write_text("".join(lines_out))


# ── MD system preparation ──────────────────────────────────────────────────────

def _prepare_md_system(
    receptor_pdb: Path,
    ligand_pdbqt: Path,
    work_dir: Path,
) -> dict:
    """Prepare an OpenMM protein–ligand system for cloud MD.

    Steps:
    1. PDBFixer: add missing residues/atoms/hydrogens, fix termini.
    2. Strip ligand PDBQT → PDB.
    3. OpenMM Modeller: combine protein + ligand, add TIP3P solvent.
    4. Create system with ff14SB + tip3p forcefield and apply HMR (4 fs).
    5. Energy minimise (500 steps).
    6. Serialise system.xml + topology.pdb to work_dir.

    Returns a dict with keys: system_xml, topology_pdb, atom_count, work_dir.

    Raises RuntimeError with pip-install instructions if openmm or pdbfixer
    are not installed.
    """
    try:
        import openmm
        from openmm import app, unit, XmlSerializer
        from openmm.app import HMassRepartitioning
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM is not installed.\n"
            "Install via: conda install -c conda-forge openmm openmmforcefields pdbfixer\n"
            "Or: pip install openmm pdbfixer"
        ) from exc

    try:
        import pdbfixer
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise RuntimeError(
            "PDBFixer is not installed.\n"
            "Install via: conda install -c conda-forge pdbfixer\n"
            "Or: pip install pdbfixer"
        ) from exc

    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fix receptor with PDBFixer
    fixer = PDBFixer(filename=str(receptor_pdb))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    fixed_receptor_pdb = work_dir / "receptor_fixed.pdb"
    with fixed_receptor_pdb.open("w") as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    # 2. Convert ligand PDBQT → PDB
    ligand_pdb = work_dir / "ligand.pdb"
    _strip_pdbqt_to_pdb(ligand_pdbqt, ligand_pdb)

    # 3. Combine protein + ligand, add TIP3P solvent
    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")

    protein_pdbfile = app.PDBFile(str(fixed_receptor_pdb))
    ligand_pdbfile = app.PDBFile(str(ligand_pdb))

    modeller = app.Modeller(protein_pdbfile.topology, protein_pdbfile.positions)
    modeller.add(ligand_pdbfile.topology, ligand_pdbfile.positions)
    modeller.addSolvent(forcefield, model="tip3p", padding=1.0 * unit.nanometer)

    # 4. Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    # 5. Apply HMR (scale hydrogen masses to 4 Da → 4 fs timestep)
    HMassRepartitioning(system, modeller.topology, hydrogenMass=4 * unit.amu)

    # 6. Energy minimise
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1 / unit.picosecond,
        0.004 * unit.picoseconds,
    )
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=500)

    # 7. Serialise system.xml and topology.pdb
    system_xml_path = work_dir / "system.xml"
    topology_pdb_path = work_dir / "topology.pdb"

    system_xml_path.write_text(XmlSerializer.serialize(system))

    state = simulation.context.getState(getPositions=True)
    with topology_pdb_path.open("w") as fh:
        app.PDBFile.writeFile(modeller.topology, state.getPositions(), fh)

    atom_count = modeller.topology.getNumAtoms()

    return {
        "system_xml": str(system_xml_path),
        "topology_pdb": str(topology_pdb_path),
        "atom_count": int(atom_count),
        "work_dir": str(work_dir),
    }


# ── Cost estimation ────────────────────────────────────────────────────────────

def _estimate_runpod_cost(atom_count: int, duration_ns: float) -> float:
    """Estimate RunPod A100 community-cloud cost in USD.

    Uses a throughput estimate of 500,000 atom·ns/hr for an A100 at 4 fs
    timestep (HMR-enabled). Rate is $1.20/hr.
    """
    wall_hours = (atom_count * duration_ns) / _THROUGHPUT_ATOM_NS_PER_HR
    return round(wall_hours * _A100_RATE_USD_PER_HR, 2)


# ── RunPod submission ──────────────────────────────────────────────────────────

def _submit_runpod_job(
    system_files: dict,
    api_key: str,
    timeout_minutes: int = 90,
) -> str:
    """Submit a serverless MD job to RunPod and return the job_id.

    Reads system.xml and topology.pdb, base64-encodes them, then submits
    to the endpoint specified by the RUNPOD_ENDPOINT_ID environment variable
    (defaults to "md-simulation-v1").
    """
    try:
        import runpod
    except ImportError as exc:
        raise RuntimeError(
            "runpod SDK is not installed.\n"
            "Install via: pip install runpod"
        ) from exc

    runpod.api_key = api_key

    system_xml_b64 = base64.b64encode(
        Path(system_files["system_xml"]).read_bytes()
    ).decode()
    topology_pdb_b64 = base64.b64encode(
        Path(system_files["topology_pdb"]).read_bytes()
    ).decode()

    endpoint_id = os.environ.get("RUNPOD_ENDPOINT_ID", _RUNPOD_DEFAULT_ENDPOINT)
    endpoint = runpod.Endpoint(endpoint_id)

    run_request = endpoint.run(
        {
            "input": {
                "system_xml_b64": system_xml_b64,
                "topology_pdb_b64": topology_pdb_b64,
                "duration_ns": _DURATION_NS,
                "timestep_fs": 4,
                "timeout_minutes": timeout_minutes,
            }
        }
    )

    return run_request.id


# ── RunPod polling ─────────────────────────────────────────────────────────────

def _poll_runpod_job(
    job_id: str,
    api_key: str,
    poll_interval_s: int = 60,
) -> dict:
    """Poll RunPod until the job reaches a terminal state.

    Returns the final status dict on COMPLETED.
    Raises RuntimeError if the job FAILED or TIMED OUT with < 15 ns completed.
    """
    try:
        import runpod
    except ImportError as exc:
        raise RuntimeError(
            "runpod SDK is not installed.\n"
            "Install via: pip install runpod"
        ) from exc

    runpod.api_key = api_key

    while True:
        status = runpod.get_job_status(job_id)
        terminal_states = ("COMPLETED", "FAILED", "TIMEOUT")
        if status["status"] in terminal_states:
            if status["status"] in ("FAILED", "TIMEOUT"):
                # Check how much of the trajectory was completed
                output = status.get("output") or {}
                completed_ns = float(output.get("completed_ns", 0))
                if completed_ns < 15:
                    raise RuntimeError(
                        f"RunPod job {job_id} {status['status']} with only "
                        f"{completed_ns:.1f} ns completed (need ≥ 15 ns). "
                        f"Error: {output.get('error', 'unknown')}"
                    )
            return status
        time.sleep(poll_interval_s)


# ── RMSD analysis ──────────────────────────────────────────────────────────────

def _compute_rmsd_pass(rmsd_csv_path: Path) -> tuple[float, bool]:
    """Analyse an RMSD trajectory CSV and return (mean_rmsd_A, passed).

    Expects columns: time_ns, rmsd_A
    Analyses the final 10 ns window.
    Raises RuntimeError if fewer than 5 data points exist in the final window
    (indicates < 15 ns trajectory was completed).
    """
    df = pd.read_csv(rmsd_csv_path)
    t_max = df["time_ns"].max()
    window = df[df["time_ns"] >= t_max - 10]

    if len(window) < 5:
        raise RuntimeError(
            f"Insufficient trajectory: < 15 ns completed "
            f"(only {len(window)} data points in final 10 ns window). "
            f"Trajectory max time: {t_max:.2f} ns"
        )

    mean_rmsd = float(window["rmsd_A"].mean())
    passed = mean_rmsd <= _RMSD_PASS_THRESHOLD_A
    return mean_rmsd, passed


# ── Gate orchestrator ──────────────────────────────────────────────────────────

def run_md_gate(
    gene_symbol: str,
    scaffold_id: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console,
    md_max_cost: float = 5.0,
) -> GateResult:
    """Orchestrate the MD pose-stability gate.

    1. Return cached result unless force=True.
    2. Check RUNPOD_API_KEY env var.
    3. Prepare system locally with OpenMM/PDBFixer.
    4. Estimate cost; abort if over cap.
    5. Submit to RunPod and poll until complete.
    6. Download RMSD CSV and compute mean over final 10 ns.
    7. Return GateResult (PASS if mean RMSD ≤ 3.0 Å, FAIL otherwise, ERROR on exception).
    """
    # 1. Cache check
    if not force:
        cached = _load_cached_gate_result(gene_symbol, scaffold_id, _GATE_NAME, cache)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] MD gate — [dim]SKIP[/dim] "
                f"(cached result for {scaffold_id})"
            )
            return cached

    # 2. Check API key
    api_key = os.environ.get("RUNPOD_API_KEY")
    if not api_key:
        result = GateResult(
            gate_name=_GATE_NAME,
            status=GateStatus.ERROR,
            score=float("nan"),
            reason=(
                "RUNPOD_API_KEY environment variable not set. "
                "Get a key at runpod.io and set: export RUNPOD_API_KEY=your_key"
            ),
            details={},
        )
        _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)
        _write_gate_report(gene_symbol, scaffold_id, _GATE_NAME, result, settings)
        return result

    # 3. Locate docking PDBQT
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

    # 4. Locate receptor PDB
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

    # 5. Prepare work directory
    work_dir = settings.cache_dir / gene_symbol / f"md_{scaffold_id}"

    # 6. Prepare MD system
    console.print(
        f"  [dim]{gene_symbol}:[/dim] MD gate — preparing system "
        f"(OpenMM + PDBFixer)..."
    )
    try:
        system_files = _prepare_md_system(receptor_pdb, pdbqt_path, work_dir)
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

    atom_count = system_files["atom_count"]

    # 7. Cost estimation
    estimated_cost = _estimate_runpod_cost(atom_count, _DURATION_NS)
    if estimated_cost > md_max_cost:
        raise RuntimeError(
            f"Estimated cost ${estimated_cost:.2f} exceeds cap ${md_max_cost:.2f}. "
            f"Use --md-max-cost to increase."
        )

    # 8. Submit to RunPod
    console.print(
        f"  [dim]{gene_symbol}:[/dim] MD gate — submitting to RunPod "
        f"({atom_count:,} atoms, ~${estimated_cost:.2f} est.)..."
    )
    job_id = _submit_runpod_job(system_files, api_key)

    # 9. Poll job
    console.print(
        f"  [dim]{gene_symbol}:[/dim] MD gate — polling job {job_id}..."
    )
    try:
        job_status = _poll_runpod_job(job_id, api_key)
    except RuntimeError as exc:
        result = GateResult(
            gate_name=_GATE_NAME,
            status=GateStatus.ERROR,
            score=float("nan"),
            reason=str(exc),
            details={"job_id": job_id, "atom_count": atom_count, "cost_usd": estimated_cost},
        )
        _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)
        _write_gate_report(gene_symbol, scaffold_id, _GATE_NAME, result, settings)
        return result

    # 10. Download RMSD CSV from job output
    rmsd_dir = settings.results_dir / gene_symbol
    rmsd_dir.mkdir(parents=True, exist_ok=True)
    rmsd_csv_path = rmsd_dir / f"validate_md_{scaffold_id}_rmsd.csv"

    output = job_status.get("output") or {}
    rmsd_csv_content = output.get("rmsd_csv", "")
    if rmsd_csv_content:
        rmsd_csv_path.write_text(rmsd_csv_content)
    else:
        result = GateResult(
            gate_name=_GATE_NAME,
            status=GateStatus.ERROR,
            score=float("nan"),
            reason=(
                f"RunPod job {job_id} completed but returned no RMSD CSV data. "
                f"Job output keys: {list(output.keys())}"
            ),
            details={"job_id": job_id, "atom_count": atom_count, "cost_usd": estimated_cost},
        )
        _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)
        _write_gate_report(gene_symbol, scaffold_id, _GATE_NAME, result, settings)
        return result

    # 11. Analyse RMSD
    mean_rmsd, passed = _compute_rmsd_pass(rmsd_csv_path)

    # 12. Build GateResult
    status = GateStatus.PASS if passed else GateStatus.FAIL
    threshold_op = "≤" if passed else ">"
    reason = (
        f"Mean RMSD (final 10 ns) = {mean_rmsd:.2f} Å "
        f"({threshold_op} threshold {_RMSD_PASS_THRESHOLD_A:.1f} Å)"
    )
    details = {
        "mean_rmsd_A": mean_rmsd,
        "atom_count": atom_count,
        "cost_usd": estimated_cost,
        "duration_ns": _DURATION_NS,
        "job_id": job_id,
        "rmsd_csv": str(rmsd_csv_path),
    }

    result = GateResult(
        gate_name=_GATE_NAME,
        status=status,
        score=mean_rmsd,
        reason=reason,
        details=details,
    )

    # 13. Write report and cache
    result.report_path = _write_gate_report(
        gene_symbol,
        scaffold_id,
        _GATE_NAME,
        result,
        settings,
        extra_sections=(
            f"## MD Simulation Configuration\n\n"
            f"- Duration: {_DURATION_NS:.0f} ns\n"
            f"- Timestep: 4 fs (HMR enabled)\n"
            f"- Forcefield: ff14SB + TIP3P\n"
            f"- Solvent padding: 1.0 nm\n"
            f"- Cloud: RunPod A100 community cloud\n"
            f"- Job ID: `{job_id}`\n"
            f"- Atom count: {atom_count:,}\n"
            f"- Estimated cost: ${estimated_cost:.2f}\n\n"
            f"## RMSD Analysis\n\n"
            f"- Mean RMSD (final 10 ns): **{mean_rmsd:.2f} Å**\n"
            f"- Pass threshold: ≤ {_RMSD_PASS_THRESHOLD_A:.1f} Å\n"
            f"- RMSD CSV: `{rmsd_csv_path.name}`\n"
        ),
    )
    _cache_gate_result(gene_symbol, scaffold_id, _GATE_NAME, result, cache)

    return result
