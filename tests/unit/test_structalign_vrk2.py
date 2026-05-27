"""Unit tests for VRK2 three-way comparison logic in structalign."""
import pandas as pd
import pytest

from pipeline.stages.structalign.structalign import (
    _classify_three_way,
    _build_comparison,
)


# ── _classify_three_way ────────────────────────────────────────────────────────

class TestClassifyThreeWay:
    def test_conserved_all_equal(self) -> None:
        assert _classify_three_way("M", "M", "M") == "conserved"

    def test_vrk1_specific_vrk1_differs_from_both(self) -> None:
        assert _classify_three_way("A", "G", "G") == "VRK1-specific"

    def test_pan_vrk_vrk1_equals_vrk2_differs_from_egfr(self) -> None:
        assert _classify_three_way("M", "M", "T") == "pan-VRK vs EGFR"

    def test_vrk2_specific_vrk2_differs_vrk1_equals_egfr(self) -> None:
        assert _classify_three_way("T", "A", "T") == "VRK2 vs VRK1+EGFR"

    def test_all_differ_returns_vrk1_specific(self) -> None:
        # When vrk1 ≠ vrk2 AND vrk1 ≠ egfr, it's VRK1-specific even if vrk2 ≠ egfr
        assert _classify_three_way("A", "G", "T") == "VRK1-specific"

    def test_gap_vrk2_treated_as_non_matching(self) -> None:
        # GAP in vrk2 means _eq(vrk1, vrk2) is False → vrk1≠vrk2
        # If vrk1 == egfr: "VRK2 vs VRK1+EGFR"
        assert _classify_three_way("M", "GAP", "M") == "VRK2 vs VRK1+EGFR"

    def test_gap_vrk1_treated_as_non_matching(self) -> None:
        # GAP in vrk1 → neither eq holds → VRK1-specific
        assert _classify_three_way("GAP", "M", "M") == "VRK1-specific"

    def test_low_conf_treated_as_non_matching(self) -> None:
        # LOW_CONF counts as non-residue → _eq returns False
        assert _classify_three_way("LOW_CONF", "M", "M") == "VRK1-specific"

    def test_low_conf_vrk2_with_equal_vrk1_egfr(self) -> None:
        assert _classify_three_way("M", "LOW_CONF", "M") == "VRK2 vs VRK1+EGFR"

    def test_gap_all_three(self) -> None:
        # All non-residue → no equalities → VRK1-specific
        assert _classify_three_way("GAP", "GAP", "GAP") == "VRK1-specific"

    def test_x_sentinel_treated_as_non_matching(self) -> None:
        assert _classify_three_way("X", "M", "M") == "VRK1-specific"


# ── _build_comparison backward compat ─────────────────────────────────────────

def _make_pocket_row(pos: int, letter: str, is_gap: bool = False) -> dict:
    return {
        "residue.klifs_id": pos,
        "residue.klifs_letter": letter,
        "is_gap": is_gap,
        "subpocket.name": "Hinge" if 46 <= pos <= 48 else "Other",
    }


def _make_pocket_df(positions: list[int], letter: str = "A") -> pd.DataFrame:
    return pd.DataFrame([_make_pocket_row(p, letter) for p in positions])


class TestBuildComparisonBackwardCompat:
    """When vrk2_klifs=None the comparison must have the same columns as spec-003."""

    SPEC_003_COLUMNS = {
        "klifs_position", "subpocket", "vrk1_aa", "egfr_aa",
        "identical_vrk1_egfr", "difference_type", "selectivity_candidate",
        "is_gatekeeper", "is_hinge", "notes",
    }
    VRK2_EXTRA_COLUMNS = {"vrk2_aa", "vrk1_vrk2_diff", "selectivity_class"}

    def test_two_way_run_has_no_vrk2_columns(self) -> None:
        vrk1 = _make_pocket_df([45, 46, 47], "M")
        egfr = _make_pocket_df([45, 46, 47], "T")
        df = _build_comparison(vrk1, egfr, vrk2_klifs=None)
        assert set(df.columns) == self.SPEC_003_COLUMNS

    def test_three_way_run_adds_vrk2_columns(self) -> None:
        vrk1 = _make_pocket_df([45, 46, 47], "M")
        egfr = _make_pocket_df([45, 46, 47], "T")
        vrk2 = _make_pocket_df([45, 46, 47], "M")
        df = _build_comparison(vrk1, egfr, vrk2_klifs=vrk2)
        assert self.VRK2_EXTRA_COLUMNS.issubset(set(df.columns))

    def test_two_way_gatekeeper_is_labeled(self) -> None:
        vrk1 = _make_pocket_df([45], "M")
        egfr = _make_pocket_df([45], "T")
        df = _build_comparison(vrk1, egfr, vrk2_klifs=None)
        gk = df[df["is_gatekeeper"]]
        assert len(gk) == 1
        assert "Gatekeeper" in gk.iloc[0]["notes"]

    def test_three_way_selectivity_class_coverage(self) -> None:
        vrk1 = _make_pocket_df([45, 46, 47, 48], "M")
        egfr = _make_pocket_df([45, 46, 47, 48], "T")
        vrk2 = _make_pocket_df([45, 46, 47, 48], "M")
        df = _build_comparison(vrk1, egfr, vrk2_klifs=vrk2)
        assert df["selectivity_class"].notna().all()
        # vrk1==vrk2 ("M"), vrk1!=egfr ("T") → pan-VRK vs EGFR
        assert (df["selectivity_class"] == "pan-VRK vs EGFR").all()

    def test_three_way_vrk1_specific_detected(self) -> None:
        # vrk1="A", vrk2="G", egfr="T" → all differ → VRK1-specific
        rows_v1 = [_make_pocket_row(45, "A"), _make_pocket_row(46, "A")]
        rows_v2 = [_make_pocket_row(45, "G"), _make_pocket_row(46, "G")]
        rows_eg = [_make_pocket_row(45, "T"), _make_pocket_row(46, "T")]
        df = _build_comparison(
            pd.DataFrame(rows_v1), pd.DataFrame(rows_eg), pd.DataFrame(rows_v2)
        )
        assert (df["selectivity_class"] == "VRK1-specific").all()
