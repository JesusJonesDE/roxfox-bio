"""Unit tests for DepMap dependency tier classification and lineage aggregation."""
import pandas as pd
import pytest

from pipeline.stages.depmap.depmap import (
    _classify_tier,
    _aggregate_lineages,
    _STRONGLY_DEP,
    _MODERATELY_DEP,
    _WEAKLY_DEP,
    _PAN_ESSENTIAL_FRACTION,
    _MIN_LINEAGE_LINES,
)


# ── dependency_tier boundary values ───────────────────────────────────────────

@pytest.mark.parametrize("effect, expected", [
    # Exactly at thresholds
    (_STRONGLY_DEP, "strongly_dependent"),        # -0.5 → strongly
    (_MODERATELY_DEP, "moderately_dependent"),     # -0.3 → moderately
    (_WEAKLY_DEP, "weakly_dependent"),             # -0.1 → weakly
    (0.0, "not_essential"),
    # Just below each threshold (more negative = same or stronger tier)
    (_STRONGLY_DEP - 0.001, "strongly_dependent"),
    (_MODERATELY_DEP - 0.001, "moderately_dependent"),
    (_WEAKLY_DEP - 0.001, "weakly_dependent"),
    # Just above each threshold (less negative = weaker tier)
    (_STRONGLY_DEP + 0.001, "moderately_dependent"),
    (_MODERATELY_DEP + 0.001, "weakly_dependent"),
    (_WEAKLY_DEP + 0.001, "not_essential"),
    # Extremes
    (-2.0, "strongly_dependent"),
    (1.0, "not_essential"),
])
def test_classify_tier(effect: float, expected: str) -> None:
    assert _classify_tier(effect) == expected


# ── pan_essential_flag ─────────────────────────────────────────────────────────

def _make_merged(effects: list[float], lineage: str = "Lung") -> pd.DataFrame:
    return pd.DataFrame({
        "ModelID": [f"ACH-{i:06d}" for i in range(len(effects))],
        "gene_effect": effects,
        "OncotreeLineage": lineage,
        "OncotreePrimaryDisease": "Test",
        "OncotreeSubtype": "Test",
        "CCLEName": [f"Cell{i}" for i in range(len(effects))],
    })


def test_pan_essential_flag_above_threshold() -> None:
    # 71% strongly dependent → pan_essential = True
    n = 100
    strong = 71
    effects = [-1.0] * strong + [0.0] * (n - strong)
    _, pan = _aggregate_lineages(_make_merged(effects))
    assert pan is True


def test_pan_essential_flag_exactly_at_threshold() -> None:
    # Exactly 70% → should be False (> not >=)
    n = 100
    effects = [-1.0] * 70 + [0.0] * 30
    _, pan = _aggregate_lineages(_make_merged(effects))
    assert pan is False


def test_pan_essential_flag_below_threshold() -> None:
    # 50% strongly dependent → pan_essential = False
    n = 100
    effects = [-1.0] * 50 + [0.0] * 50
    _, pan = _aggregate_lineages(_make_merged(effects))
    assert pan is False


# ── lineage filtering ──────────────────────────────────────────────────────────

def test_lineage_excluded_when_fewer_than_min_lines() -> None:
    # Create two lineages: one with MIN_LINEAGE_LINES-1 lines, one with MIN_LINEAGE_LINES
    n_small = _MIN_LINEAGE_LINES - 1
    n_large = _MIN_LINEAGE_LINES
    effects_small = [-0.6] * n_small
    effects_large = [-0.6] * n_large

    merged = pd.concat([
        _make_merged(effects_small, lineage="Small"),
        _make_merged(effects_large, lineage="Large"),
    ], ignore_index=True)

    ranked, _ = _aggregate_lineages(merged)
    assert "Small" not in ranked["OncotreeLineage"].values
    assert "Large" in ranked["OncotreeLineage"].values


def test_ranked_sorted_by_median_effect_ascending() -> None:
    merged = pd.concat([
        _make_merged([-0.2, -0.2, -0.2], lineage="Weak"),
        _make_merged([-0.8, -0.8, -0.8], lineage="Strong"),
        _make_merged([-0.5, -0.5, -0.5], lineage="Mid"),
    ], ignore_index=True)

    ranked, _ = _aggregate_lineages(merged)
    assert list(ranked["OncotreeLineage"]) == ["Strong", "Mid", "Weak"]


def test_pct_strongly_dependent_calculation() -> None:
    # 2 out of 4 lines at -0.6 = 50%
    effects = [-0.6, -0.6, 0.0, 0.0]
    merged = _make_merged(effects)
    ranked, _ = _aggregate_lineages(merged)
    assert ranked.iloc[0]["pct_strongly_dependent"] == 50.0
