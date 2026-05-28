"""Unit tests for the MD pose-stability validation gate."""
from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from pipeline.models import GateResult, GateStatus


# ── Helper factories ───────────────────────────────────────────────────────────

def _make_mocks():
    """Return (settings, cache, console) mocks for run_md_gate."""
    settings = MagicMock()
    settings.cache_dir = Path("/tmp/test_md_cache")
    settings.results_dir = Path("/tmp/test_md_results")

    cache = MagicMock()
    cache.load.return_value = None  # no cached result by default

    console = MagicMock()
    return settings, cache, console


def _make_rmsd_csv(time_ns_values, rmsd_values) -> Path:
    """Write a temporary RMSD CSV and return its path."""
    import tempfile
    tmp = Path(tempfile.mktemp(suffix=".csv"))
    df = pd.DataFrame({"time_ns": time_ns_values, "rmsd_A": rmsd_values})
    df.to_csv(tmp, index=False)
    return tmp


# ── Test 1: cost estimation ────────────────────────────────────────────────────

def test_estimate_runpod_cost_50k_atoms_20ns():
    """50,000 atoms × 20 ns → wall_hours = (50000*20)/500000 = 2 hr → $2.40."""
    from pipeline.stages.validate.gates.md import _estimate_runpod_cost

    cost = _estimate_runpod_cost(atom_count=50_000, duration_ns=20)
    assert cost == pytest.approx(2.40)


# ── Test 2: RMSD pass (mean 2.5 Å ≤ 3.0 Å threshold) ────────────────────────

def test_compute_rmsd_pass_mean_2_5():
    """Final 10 ns window with mean RMSD 2.5 Å should return pass=True."""
    from pipeline.stages.validate.gates.md import _compute_rmsd_pass

    # 20 data points spanning t=10..20 ns; all rmsd = 2.5
    time_ns = [10.0 + i * 0.5 for i in range(21)]  # 10.0 to 20.0
    rmsd_A = [2.5] * len(time_ns)
    csv_path = _make_rmsd_csv(time_ns, rmsd_A)

    try:
        mean_rmsd, passed = _compute_rmsd_pass(csv_path)
    finally:
        csv_path.unlink(missing_ok=True)

    assert mean_rmsd == pytest.approx(2.5)
    assert passed is True


# ── Test 3: RMSD fail (mean 3.5 Å > 3.0 Å threshold) ────────────────────────

def test_compute_rmsd_fail_mean_3_5():
    """Final 10 ns window with mean RMSD 3.5 Å should return pass=False."""
    from pipeline.stages.validate.gates.md import _compute_rmsd_pass

    time_ns = [10.0 + i * 0.5 for i in range(21)]  # 10.0 to 20.0
    rmsd_A = [3.5] * len(time_ns)
    csv_path = _make_rmsd_csv(time_ns, rmsd_A)

    try:
        mean_rmsd, passed = _compute_rmsd_pass(csv_path)
    finally:
        csv_path.unlink(missing_ok=True)

    assert mean_rmsd == pytest.approx(3.5)
    assert passed is False


# ── Test 4: missing RUNPOD_API_KEY → GateResult ERROR ────────────────────────

@patch.dict("os.environ", {}, clear=True)  # ensure RUNPOD_API_KEY is absent
@patch("pipeline.stages.validate.validate._write_gate_report")
@patch("pipeline.stages.validate.validate._cache_gate_result")
@patch("pipeline.stages.validate.validate._load_cached_gate_result")
def test_run_md_gate_missing_api_key(
    mock_load_cached,
    mock_cache_gate,
    mock_write_report,
):
    """run_md_gate without RUNPOD_API_KEY returns a GateResult with ERROR status."""
    from pipeline.stages.validate.gates.md import run_md_gate

    mock_load_cached.return_value = None
    mock_write_report.return_value = Path("/tmp/test_md_results/gene/validate_md_SCF-001.md")

    settings, cache, console = _make_mocks()

    result = run_md_gate(
        "GENE1", "SCF-001", settings, cache, force=True, console=console
    )

    assert isinstance(result, GateResult)
    assert result.status == GateStatus.ERROR
    assert "RUNPOD_API_KEY" in result.reason
    assert "runpod.io" in result.reason


# ── Test 5: cost cap exceeded → RuntimeError ─────────────────────────────────

@patch.dict("os.environ", {"RUNPOD_API_KEY": "test-key-abc"})
@patch("pipeline.stages.validate.validate._write_gate_report")
@patch("pipeline.stages.validate.validate._cache_gate_result")
@patch("pipeline.stages.validate.validate._load_cached_gate_result")
@patch("pipeline.stages.validate.gates.md._load_docking_pdbqt")
@patch("pipeline.stages.validate.gates.md._prepare_md_system")
def test_run_md_gate_cost_cap_exceeded(
    mock_prepare,
    mock_load_pdbqt,
    mock_load_cached,
    mock_cache_gate,
    mock_write_report,
    tmp_path,
):
    """When atom_count=1_000_000, estimated cost > $5 → RuntimeError raised."""
    from pipeline.stages.validate.gates.md import run_md_gate, _estimate_runpod_cost

    # Sanity-check the cost formula directly: 1M atoms × 20 ns = $48
    assert _estimate_runpod_cost(1_000_000, 20) > 5.0

    mock_load_cached.return_value = None
    mock_load_pdbqt.return_value = tmp_path / "docking_poses_SCF-002.pdbqt"
    mock_prepare.return_value = {
        "system_xml": str(tmp_path / "system.xml"),
        "topology_pdb": str(tmp_path / "topology.pdb"),
        "atom_count": 1_000_000,
        "work_dir": str(tmp_path),
    }
    mock_write_report.return_value = tmp_path / "validate_md_SCF-002.md"

    settings, cache, console = _make_mocks()
    # Point structures_dir to a real temp dir with a dummy .pdb file
    structures_dir = tmp_path / "GENE1" / "structures"
    structures_dir.mkdir(parents=True)
    (structures_dir / "receptor.pdb").write_text("ATOM fake\n")
    settings.cache_dir = tmp_path

    with pytest.raises(RuntimeError, match="Estimated cost.*exceeds cap"):
        run_md_gate(
            "GENE1", "SCF-002", settings, cache, force=True, console=console,
            md_max_cost=5.0,
        )
