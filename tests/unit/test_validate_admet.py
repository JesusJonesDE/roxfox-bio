"""Unit tests for the ADMET validation gate."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from pipeline.models import GateResult, GateStatus


# ── Helpers ────────────────────────────────────────────────────────────────────

def _all_pass_preds() -> dict:
    """Return a flat ADMET-AI predictions dict where all 6 properties pass.

    Field names match admet-ai >= 1.0 TDC benchmark naming.
    """
    return {
        "BBB_Martins": 0.75,          # pass: > 0.5
        "CYP1A2_Veith": 0.10,         # pass: < 0.3
        "CYP2D6_Veith": 0.15,         # pass: < 0.3
        "CYP3A4_Veith": 0.20,         # pass: < 0.3
        "Solubility_AqSolDB": -2.5,   # pass: > -4
        "HIA_Hou": 0.85,              # pass: > 0.3
    }


def _make_mocks(smiles: str = "c1ccccc1", preds: dict | None = None):
    """Return (settings, cache, console) mocks wired for run_admet_gate."""
    settings = MagicMock()
    settings.results_dir = Path("/tmp/test_admet_results")

    cache = MagicMock()
    cache.load.return_value = None  # no cached result by default

    console = MagicMock()

    return settings, cache, console


# ── Test 1: all-pass case ──────────────────────────────────────────────────────

@patch("pipeline.stages.validate.gates.admet._write_gate_report")
@patch("pipeline.stages.validate.gates.admet._cache_gate_result")
@patch("pipeline.stages.validate.gates.admet._load_smiles")
@patch("admet_ai.ADMETModel")
def test_admet_gate_all_pass(
    MockADMETModel, mock_load_smiles, mock_cache_gate, mock_write_report
):
    from pipeline.stages.validate.gates.admet import run_admet_gate

    preds = _all_pass_preds()
    MockADMETModel.return_value.predict.return_value = preds
    mock_load_smiles.return_value = "c1ccccc1"
    mock_write_report.return_value = Path("/tmp/test_admet_results/gene/validate_admet_SCF-001.md")

    settings, cache, console = _make_mocks()

    result = run_admet_gate("GENE1", "SCF-001", settings, cache, force=True, console=console)

    assert isinstance(result, GateResult)
    assert result.status == GateStatus.PASS
    assert result.score == pytest.approx(0.75)
    assert result.reason == "All 6 ADMET properties passed"
    assert len(result.details) == 6


# ── Test 2: fail case — BBB_Martini below threshold ───────────────────────────

@patch("pipeline.stages.validate.gates.admet._write_gate_report")
@patch("pipeline.stages.validate.gates.admet._cache_gate_result")
@patch("pipeline.stages.validate.gates.admet._load_smiles")
@patch("admet_ai.ADMETModel")
def test_admet_gate_fail_bbb(
    MockADMETModel, mock_load_smiles, mock_cache_gate, mock_write_report
):
    from pipeline.stages.validate.gates.admet import run_admet_gate

    preds = _all_pass_preds()
    preds["BBB_Martins"] = 0.20  # fail: threshold > 0.5
    MockADMETModel.return_value.predict.return_value = preds
    mock_load_smiles.return_value = "c1ccccc1"
    mock_write_report.return_value = Path("/tmp/validate_admet_SCF-002.md")

    settings, cache, console = _make_mocks()

    result = run_admet_gate("GENE1", "SCF-002", settings, cache, force=True, console=console)

    assert result.status == GateStatus.FAIL
    assert "BBB_Martins" in result.reason
    assert result.score == pytest.approx(0.20)


# ── Test 3: salt SMILES — largest fragment processed without error ─────────────

@patch("pipeline.stages.validate.gates.admet._write_gate_report")
@patch("pipeline.stages.validate.gates.admet._cache_gate_result")
@patch("pipeline.stages.validate.gates.admet._load_smiles")
@patch("admet_ai.ADMETModel")
def test_admet_gate_salt_smiles(
    MockADMETModel, mock_load_smiles, mock_cache_gate, mock_write_report
):
    """SMILES with counterion (Cl.) should be stripped to the largest fragment.

    _load_smiles already applies salt-stripping, so the gate sees a clean SMILES.
    This test verifies the gate runs without error when the returned SMILES is
    the stripped drug fragment.
    """
    from pipeline.stages.validate.gates.admet import run_admet_gate

    # Simulate _load_smiles returning the stripped SMILES (no 'Cl.' prefix)
    stripped = "NCCNc1ncc(C(N)=O)c(Nc2cccc(C(F)(F)F)c2)n1"
    mock_load_smiles.return_value = stripped

    preds = _all_pass_preds()
    MockADMETModel.return_value.predict.return_value = preds
    mock_write_report.return_value = Path("/tmp/validate_admet_SCF-003.md")

    settings, cache, console = _make_mocks()

    result = run_admet_gate("GENE1", "SCF-003", settings, cache, force=True, console=console)

    # No exception should be raised; result should be valid
    assert isinstance(result, GateResult)
    assert result.status in (GateStatus.PASS, GateStatus.FAIL)
    MockADMETModel.return_value.predict.assert_called_once_with(smiles=stripped)


# ── Test 4: cache hit — second call returns cached result, model not called ───

@patch("pipeline.stages.validate.gates.admet._write_gate_report")
@patch("pipeline.stages.validate.gates.admet._cache_gate_result")
@patch("pipeline.stages.validate.gates.admet._load_cached_gate_result")
@patch("pipeline.stages.validate.gates.admet._load_smiles")
@patch("admet_ai.ADMETModel")
def test_admet_gate_cache_hit(
    MockADMETModel, mock_load_smiles, mock_load_cached, mock_cache_gate, mock_write_report
):
    from pipeline.stages.validate.gates.admet import run_admet_gate

    cached_result = GateResult(
        gate_name="admet",
        status=GateStatus.PASS,
        score=0.75,
        reason="All 6 ADMET properties passed",
        details=_all_pass_preds(),
    )
    mock_load_cached.return_value = cached_result

    settings, cache, console = _make_mocks()

    # First call with force=False — should return from cache
    result = run_admet_gate("GENE1", "SCF-004", settings, cache, force=False, console=console)

    assert result.status == GateStatus.PASS
    assert result.score == pytest.approx(0.75)

    # ADMETModel should never have been instantiated or called
    MockADMETModel.assert_not_called()
    mock_load_smiles.assert_not_called()
    mock_write_report.assert_not_called()
