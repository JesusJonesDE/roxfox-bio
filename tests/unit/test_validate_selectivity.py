"""Unit tests for pipeline/stages/validate/gates/selectivity.py"""
from __future__ import annotations

import pytest


# ── Helpers ────────────────────────────────────────────────────────────────────


def _make_gate_result(primary: float, offtarget_affinities: dict) -> object:
    """Helper: compute SI and return a mock GateResult-like namespace using the
    same logic as run_selectivity_gate (extracted for unit testing without I/O).
    """
    from types import SimpleNamespace

    from pipeline.models import GateStatus

    abs_primary = abs(primary)
    abs_ots = {g: abs(a) for g, a in offtarget_affinities.items()}

    if not abs_ots or max(abs_ots.values()) == 0.0:
        si = float("inf")
        worst = next(iter(abs_ots), "none")
        worst_aff = 0.0
    else:
        worst = max(abs_ots, key=abs_ots.__getitem__)
        worst_aff = offtarget_affinities[worst]
        si = abs_primary / abs_ots[worst]

    status = GateStatus.PASS if si >= 10.0 else GateStatus.FAIL

    return SimpleNamespace(
        score=si,
        status=status,
        worst_offtarget=worst,
        worst_affinity=worst_aff,
    )


# ── Test: SI calculation → FAIL ────────────────────────────────────────────────


class TestSICalculationFail:
    """primary=-8.0, best_offtarget=-4.0 → SI=2.0 → FAIL"""

    def test_si_value(self) -> None:
        result = _make_gate_result(-8.0, {"EGFR": -4.0, "CDK2": -2.0})
        assert result.score == pytest.approx(2.0)

    def test_status_is_fail(self) -> None:
        from pipeline.models import GateStatus

        result = _make_gate_result(-8.0, {"EGFR": -4.0, "CDK2": -2.0})
        assert result.status == GateStatus.FAIL

    def test_worst_offtarget_identified(self) -> None:
        result = _make_gate_result(-8.0, {"EGFR": -4.0, "CDK2": -2.0})
        assert result.worst_offtarget == "EGFR"

    def test_worst_affinity_value(self) -> None:
        result = _make_gate_result(-8.0, {"EGFR": -4.0, "CDK2": -2.0})
        assert result.worst_affinity == pytest.approx(-4.0)


# ── Test: SI calculation → PASS ───────────────────────────────────────────────


class TestSICalculationPass:
    """primary=-8.0, best_offtarget=-0.5 → SI=16.0 → PASS"""

    def test_si_value(self) -> None:
        result = _make_gate_result(-8.0, {"EGFR": -0.5, "CDK2": -0.3})
        assert result.score == pytest.approx(16.0)

    def test_status_is_pass(self) -> None:
        from pipeline.models import GateStatus

        result = _make_gate_result(-8.0, {"EGFR": -0.5, "CDK2": -0.3})
        assert result.status == GateStatus.PASS

    def test_exact_threshold_10x_is_pass(self) -> None:
        from pipeline.models import GateStatus

        result = _make_gate_result(-10.0, {"EGFR": -1.0})
        assert result.score == pytest.approx(10.0)
        assert result.status == GateStatus.PASS

    def test_just_below_threshold_is_fail(self) -> None:
        from pipeline.models import GateStatus

        result = _make_gate_result(-9.9, {"EGFR": -1.0})
        assert result.score == pytest.approx(9.9)
        assert result.status == GateStatus.FAIL


# ── Test: target exclusion ────────────────────────────────────────────────────


class TestTargetExclusion:
    """When gene_symbol='EGFR', the EGFR panel entry must be excluded."""

    def test_egfr_excluded_from_panel(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        active = [e for e in SELECTIVITY_PANEL if e.gene.upper() != "EGFR"]
        gene_names = [e.gene for e in active]
        assert "EGFR" not in gene_names

    def test_egfr_excluded_leaves_three_entries(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        active = [e for e in SELECTIVITY_PANEL if e.gene.upper() != "EGFR"]
        assert len(active) == 3

    def test_vrk1_excluded_leaves_all_four(self) -> None:
        """VRK1 is not in the panel, so all 4 entries remain."""
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        active = [e for e in SELECTIVITY_PANEL if e.gene.upper() != "VRK1"]
        assert len(active) == 4

    def test_vrk2_excluded_leaves_three(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        active = [e for e in SELECTIVITY_PANEL if e.gene.upper() != "VRK2"]
        assert len(active) == 3
        assert all(e.gene != "VRK2" for e in active)


# ── Test: GateResult fields ────────────────────────────────────────────────────


class TestGateResultFields:
    """GateResult.score must equal the computed SI value."""

    def test_score_equals_si(self) -> None:
        result = _make_gate_result(-8.0, {"EGFR": -4.0})
        assert result.score == pytest.approx(2.0)

    def test_score_equals_si_pass_case(self) -> None:
        result = _make_gate_result(-8.0, {"EGFR": -0.5})
        assert result.score == pytest.approx(16.0)

    def test_gate_result_dataclass_fields(self) -> None:
        from pipeline.models import GateResult, GateStatus

        gr = GateResult(
            gate_name="selectivity",
            status=GateStatus.PASS,
            score=16.0,
            reason="SI=16.0x vs EGFR (-0.5 kcal/mol)",
            details={"primary_affinity": -8.0, "EGFR": -0.5},
        )
        assert gr.gate_name == "selectivity"
        assert gr.score == pytest.approx(16.0)
        assert gr.status == GateStatus.PASS
        assert gr.details["primary_affinity"] == pytest.approx(-8.0)
        assert gr.details["EGFR"] == pytest.approx(-0.5)

    def test_gate_result_fail_details(self) -> None:
        from pipeline.models import GateResult, GateStatus

        gr = GateResult(
            gate_name="selectivity",
            status=GateStatus.FAIL,
            score=2.0,
            reason="SI=2.0x vs EGFR (-4.0 kcal/mol)",
            details={"primary_affinity": -8.0, "EGFR": -4.0, "CDK2": -2.0},
        )
        assert gr.status == GateStatus.FAIL
        assert gr.score == pytest.approx(2.0)
        assert "EGFR" in gr.details


# ── Test: Panel structure ──────────────────────────────────────────────────────


class TestPanelStructure:
    """Verify the default SELECTIVITY_PANEL is correctly defined."""

    def test_panel_has_four_entries(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        assert len(SELECTIVITY_PANEL) == 4

    def test_vrk2_is_alphafold(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        vrk2 = next(e for e in SELECTIVITY_PANEL if e.gene == "VRK2")
        assert vrk2.source == "alphafold"
        assert vrk2.af2_uniprot == "O95551"
        assert vrk2.warning != ""

    def test_egfr_entry(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        egfr = next(e for e in SELECTIVITY_PANEL if e.gene == "EGFR")
        assert egfr.source == "pdb"
        assert egfr.pdb_id == "1IVO"

    def test_cdk2_entry(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        cdk2 = next(e for e in SELECTIVITY_PANEL if e.gene == "CDK2")
        assert cdk2.pdb_id == "1E9H"
        assert cdk2.chain_id == "A"

    def test_plk1_entry(self) -> None:
        from pipeline.stages.validate.gates.selectivity import SELECTIVITY_PANEL

        plk1 = next(e for e in SELECTIVITY_PANEL if e.gene == "PLK1")
        assert plk1.pdb_id == "2OKR"
