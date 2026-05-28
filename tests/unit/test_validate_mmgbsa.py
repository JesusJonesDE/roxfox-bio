"""Unit tests for the MM-GBSA validation gate."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from pipeline.models import GateResult, GateStatus


# ── Helpers ───────────────────────────────────────────────────────────────────

_SAMPLE_PDBQT = """\
MODEL 1
ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00     0.000 C
HETATM    2  O1  LIG A   1       4.000   5.000   6.000  1.00  0.00    -0.500 OA
ATOM      3  N1  LIG A   1       7.000   8.000   9.000  1.00  0.00    -0.200 NA
ENDMDL
MODEL 2
ATOM      1  C1  LIG A   1      10.000  11.000  12.000  1.00  0.00     0.000 C
ENDMDL
"""


# ── Test _pdbqt_to_pdb ────────────────────────────────────────────────────────

class TestPdbqtToPdb:
    def test_output_has_only_atom_and_hetatm_lines(self, tmp_path: Path) -> None:
        """MODEL/ENDMDL/blank lines must not appear in output."""
        from pipeline.stages.validate.gates.mmgbsa import _pdbqt_to_pdb

        pdbqt = tmp_path / "pose.pdbqt"
        pdbqt.write_text(_SAMPLE_PDBQT)
        out_pdb = tmp_path / "ligand.pdb"

        _pdbqt_to_pdb(pdbqt, out_pdb)

        lines = [l for l in out_pdb.read_text().splitlines() if l.strip()]
        record_types = {l[:6].strip() for l in lines}
        assert record_types <= {"ATOM", "HETATM"}, (
            f"Unexpected record types: {record_types - {'ATOM', 'HETATM'}}"
        )

    def test_only_first_model_block_is_included(self, tmp_path: Path) -> None:
        """Lines from MODEL 2 must not appear in output."""
        from pipeline.stages.validate.gates.mmgbsa import _pdbqt_to_pdb

        pdbqt = tmp_path / "pose.pdbqt"
        pdbqt.write_text(_SAMPLE_PDBQT)
        out_pdb = tmp_path / "ligand.pdb"

        _pdbqt_to_pdb(pdbqt, out_pdb)

        content = out_pdb.read_text()
        # 3 ATOM/HETATM lines from MODEL 1; MODEL 2 has 1 extra ATOM
        lines = [l for l in content.splitlines() if l.strip()]
        assert len(lines) == 3

    def test_pdbqt_extra_columns_stripped(self, tmp_path: Path) -> None:
        """Output lines must not contain the PDBQT charge/atom-type columns."""
        from pipeline.stages.validate.gates.mmgbsa import _pdbqt_to_pdb

        pdbqt = tmp_path / "pose.pdbqt"
        pdbqt.write_text(_SAMPLE_PDBQT)
        out_pdb = tmp_path / "ligand.pdb"

        _pdbqt_to_pdb(pdbqt, out_pdb)

        for line in out_pdb.read_text().splitlines():
            if not line.strip():
                continue
            # Each output line must be no longer than 68 chars (stripped)
            assert len(line.rstrip("\n")) <= 68, (
                f"Line too long (PDBQT columns not stripped): {line!r}"
            )

    def test_empty_pdbqt_produces_empty_pdb(self, tmp_path: Path) -> None:
        from pipeline.stages.validate.gates.mmgbsa import _pdbqt_to_pdb

        pdbqt = tmp_path / "empty.pdbqt"
        pdbqt.write_text("")
        out_pdb = tmp_path / "out.pdb"

        _pdbqt_to_pdb(pdbqt, out_pdb)

        assert out_pdb.read_text().strip() == ""


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _make_settings(tmp_path: Path):
    """Return a minimal Settings-like object pointing at tmp_path."""
    settings = MagicMock()
    settings.cache_dir = tmp_path / "cache"
    settings.results_dir = tmp_path / "results"
    return settings


def _make_cache():
    cache = MagicMock()
    cache.load.return_value = None  # No cached result by default
    return cache


def _make_console():
    return MagicMock()


def _setup_gate_env(tmp_path: Path, gene: str = "VRK1", scaffold: str = "SCF-009"):
    """Create the minimum directory/file structure that the gate requires."""
    settings = _make_settings(tmp_path)

    # Receptor PDB
    structures_dir = settings.cache_dir / gene / "structures"
    structures_dir.mkdir(parents=True)
    receptor_pdb = structures_dir / "receptor.pdb"
    receptor_pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")

    # Docking PDBQT
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True)
    pdbqt = results_dir / f"docking_poses_{scaffold}.pdbqt"
    pdbqt.write_text(_SAMPLE_PDBQT)

    return settings


# ── Test run_mmgbsa_gate: PASS case ──────────────────────────────────────────

class TestRunMmgbsaGatePass:
    def test_gate_returns_pass_when_delta_g_below_threshold(self, tmp_path: Path) -> None:
        """ΔG = -8.5 → PASS (threshold -7.0)."""
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                return_value=-8.5,
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.status == GateStatus.PASS
        assert result.score == pytest.approx(-8.5)
        assert result.details["delta_G_kcal_mol"] == pytest.approx(-8.5)

    def test_gate_score_equals_delta_g(self, tmp_path: Path) -> None:
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                return_value=-9.0,
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.score == pytest.approx(-9.0)

    def test_gate_caches_result_on_pass(self, tmp_path: Path) -> None:
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                return_value=-8.5,
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._cache_gate_result"
            ) as mock_cache_save,
        ):
            run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        mock_cache_save.assert_called_once()


# ── Test run_mmgbsa_gate: FAIL case ──────────────────────────────────────────

class TestRunMmgbsaGateFail:
    def test_gate_returns_fail_when_delta_g_above_threshold(self, tmp_path: Path) -> None:
        """ΔG = -6.0 → FAIL (threshold -7.0)."""
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                return_value=-6.0,
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.status == GateStatus.FAIL
        assert result.score == pytest.approx(-6.0)

    def test_gate_returns_fail_exactly_at_threshold(self, tmp_path: Path) -> None:
        """ΔG = -7.0 → PASS (boundary: ≤ -7.0)."""
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                return_value=-7.0,
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.status == GateStatus.PASS

    def test_gate_fail_score_reflects_delta_g(self, tmp_path: Path) -> None:
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                return_value=-4.2,
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.status == GateStatus.FAIL
        assert result.score == pytest.approx(-4.2)


# ── Test run_mmgbsa_gate: ERROR case ─────────────────────────────────────────

class TestRunMmgbsaGateError:
    def test_gate_returns_error_when_gmx_mmpbsa_raises(self, tmp_path: Path) -> None:
        """RuntimeError from _run_gmx_mmpbsa → GateStatus.ERROR."""
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                side_effect=RuntimeError("gmx_MMPBSA is not installed or not on PATH."),
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.status == GateStatus.ERROR
        assert "gmx_MMPBSA" in result.reason

    def test_gate_error_score_is_nan(self, tmp_path: Path) -> None:
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate
        import math

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                side_effect=RuntimeError("calculation failed"),
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        assert result.status == GateStatus.ERROR
        assert math.isnan(result.score)

    def test_gate_error_exit_code_via_status(self, tmp_path: Path) -> None:
        """ERROR status corresponds to exit code 1 (checked by the orchestrator)."""
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _setup_gate_env(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with (
            patch(
                "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa",
                side_effect=RuntimeError("MM-GBSA calculation failed: subprocess error"),
            ),
            patch(
                "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
                return_value=tmp_path / "report.md",
            ),
        ):
            result = run_mmgbsa_gate("VRK1", "SCF-009", settings, cache, force=True, console=console)

        # The orchestrator (run_validate) maps GateStatus.ERROR → exit_code = 1
        assert result.status == GateStatus.ERROR

    def test_gate_returns_error_when_docking_pdbqt_missing(self, tmp_path: Path) -> None:
        """Missing docking PDBQT → GateStatus.ERROR with helpful message."""
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        settings = _make_settings(tmp_path)
        # Create structures dir but no results dir / pdbqt
        structures_dir = settings.cache_dir / "VRK1" / "structures"
        structures_dir.mkdir(parents=True)
        (structures_dir / "receptor.pdb").write_text("ATOM\n")
        (settings.results_dir / "VRK1").mkdir(parents=True)

        cache = _make_cache()
        console = _make_console()

        with patch(
            "pipeline.stages.validate.gates.mmgbsa._write_gate_report",
            return_value=tmp_path / "report.md",
        ):
            result = run_mmgbsa_gate(
                "VRK1", "MISSING-001", settings, cache, force=True, console=console
            )

        assert result.status == GateStatus.ERROR
        assert "MISSING-001" in result.reason or "docking" in result.reason.lower()


# ── Test cache hit ────────────────────────────────────────────────────────────

class TestRunMmgbsaGateCacheHit:
    def test_cached_result_returned_without_running_calculation(self, tmp_path: Path) -> None:
        from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate

        cached_result = GateResult(
            gate_name="mmgbsa",
            status=GateStatus.PASS,
            score=-8.5,
            reason="ΔG = -8.50 kcal/mol (≤ threshold -7.0 kcal/mol)",
            details={"delta_G_kcal_mol": -8.5},
        )

        settings = _make_settings(tmp_path)
        cache = _make_cache()
        console = _make_console()

        with patch(
            "pipeline.stages.validate.gates.mmgbsa._load_cached_gate_result",
            return_value=cached_result,
        ) as mock_load, patch(
            "pipeline.stages.validate.gates.mmgbsa._run_gmx_mmpbsa"
        ) as mock_run:
            result = run_mmgbsa_gate(
                "VRK1", "SCF-009", settings, cache, force=False, console=console
            )

        mock_load.assert_called_once()
        mock_run.assert_not_called()
        assert result.status == GateStatus.PASS
        assert result.score == pytest.approx(-8.5)
