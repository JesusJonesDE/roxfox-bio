"""Unit tests for the fragment output step."""
from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ── Helpers ────────────────────────────────────────────────────────────────────

_EXPECTED_COLUMNS = [
    "molecule_id",
    "smiles",
    "best_value_nm",
    "best_assay_type",
    "molecular_weight",
    "logp",
    "hbd",
    "hba",
    "rotatable_bonds",
    "ro5_violations",
    "passes_ro5",
    "scaffold_id",
    "source",
    "off_target_flags",
    "selectivity_flag",
]


def _make_candidates_df(n: int = 3) -> pd.DataFrame:
    """Return a minimal candidates DataFrame for testing."""
    rows = []
    for i in range(1, n + 1):
        rows.append(
            {
                "candidate_id": f"IGHMBP2-SCF-{i:03d}",
                "smiles": "c1ccccc1",
                "parent_fragment_id": "FRAG-001",
                "molecular_weight": 340.0 + i,
                "logp": 2.5,
                "hbd": 1,
                "hba": 2,
                "rotatable_bonds": 3,
                "ro5_violations": 0,
                "passes_ro5": True,
                "sa_score": 2.5,
                "grow_method": "smarts",
            }
        )
    return pd.DataFrame(rows)


def _make_mocks(tmp_path: Path):
    """Return wired (settings, cache, console) mocks."""
    settings = MagicMock()
    settings.results_dir = tmp_path
    cache = MagicMock()
    cache.load.return_value = None  # no cached result
    console = MagicMock()
    return settings, cache, console


# ── Test 1: Output CSV has exact column names matching VRK1-like schema ────────

class TestOutputColumns:
    def test_column_names_match_expected_schema(self, tmp_path: Path) -> None:
        """Written compounds_filtered.csv must have the exact column order."""
        from pipeline.stages.fragment.output import run_output, _OUTPUT_COLUMNS

        assert _OUTPUT_COLUMNS == _expected_columns(), (
            "Module _OUTPUT_COLUMNS must match the expected schema"
        )

    def test_written_csv_has_correct_headers(self, tmp_path: Path) -> None:
        """The CSV file on disk must have the exact expected column names."""
        from pipeline.stages.fragment.output import run_output

        settings, cache, console = _make_mocks(tmp_path / "results")
        gene_dir = tmp_path / "results" / "IGHMBP2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        candidates_df = _make_candidates_df(2)

        with patch(
            "pipeline.stages.fragment.output._run_admet_relaxed",
            return_value=(True, 0.65, {}),
        ):
            run_output("IGHMBP2", candidates_df, settings, cache, force=True, console=console)

        out_csv = gene_dir / "compounds_filtered.csv"
        assert out_csv.exists(), "compounds_filtered.csv must be written"

        with open(out_csv) as fh:
            reader = csv.reader(fh)
            header = next(reader)

        assert header == _EXPECTED_COLUMNS, (
            f"CSV header mismatch.\n  Got:      {header}\n  Expected: {_EXPECTED_COLUMNS}"
        )


def _expected_columns() -> list[str]:
    return _EXPECTED_COLUMNS


# ── Test 2: BBB threshold 0.3 used (not 0.5) ──────────────────────────────────

class TestBBBThreshold:
    def test_relaxed_thresholds_bbb_is_0_3(self) -> None:
        """_RELAXED_THRESHOLDS must use BBB threshold 0.3, not 0.5."""
        from pipeline.stages.fragment.output import _RELAXED_THRESHOLDS

        op, threshold = _RELAXED_THRESHOLDS["BBB_Martins"]
        assert op == ">", "BBB operator must be '>'"
        assert threshold == pytest.approx(0.3), (
            f"BBB threshold must be 0.3, got {threshold}"
        )
        assert threshold != pytest.approx(0.5), (
            "BBB threshold must NOT be 0.5 (the standard admet gate value)"
        )

    def test_admet_model_receives_smiles_and_uses_relaxed_threshold(
        self, tmp_path: Path
    ) -> None:
        """ADMETModel.predict must be called per candidate with relaxed BBB=0.3."""
        from pipeline.stages.fragment.output import run_output

        settings, cache, console = _make_mocks(tmp_path / "results")
        gene_dir = tmp_path / "results" / "IGHMBP2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        candidates_df = _make_candidates_df(1)
        captured_smiles: list[str] = []

        def mock_admet(smiles: str):
            captured_smiles.append(smiles)
            # BBB score 0.35 — above relaxed threshold 0.3 but below standard 0.5
            return (True, 0.35, {"BBB_Martins": 0.35})

        with patch(
            "pipeline.stages.fragment.output._run_admet_relaxed",
            side_effect=mock_admet,
        ):
            n = run_output(
                "IGHMBP2", candidates_df, settings, cache, force=True, console=console
            )

        assert n == 1, "Should write 1 row"
        assert len(captured_smiles) == 1
        assert captured_smiles[0] == "c1ccccc1"

    def test_standard_0_5_threshold_would_reject_passes_relaxed(self) -> None:
        """Verify the threshold logic: BBB=0.35 passes 0.3 but fails 0.5."""
        from pipeline.stages.fragment.output import _check_threshold

        assert _check_threshold(0.35, ">", 0.3) is True, "0.35 > 0.3 should pass"
        assert _check_threshold(0.35, ">", 0.5) is False, "0.35 > 0.5 should fail"


# ── Test 3: Empty candidates_df → empty CSV with header only, returns 0 ────────

class TestEmptyCandidates:
    def test_empty_df_writes_header_only_csv(self, tmp_path: Path) -> None:
        """Empty candidates_df should write an empty CSV with header row."""
        from pipeline.stages.fragment.output import run_output

        settings, cache, console = _make_mocks(tmp_path / "results")
        gene_dir = tmp_path / "results" / "IGHMBP2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        empty_df = pd.DataFrame(
            columns=[
                "candidate_id", "smiles", "parent_fragment_id",
                "molecular_weight", "logp", "hbd", "hba",
                "rotatable_bonds", "ro5_violations", "passes_ro5",
                "sa_score", "grow_method",
            ]
        )

        n = run_output("IGHMBP2", empty_df, settings, cache, force=True, console=console)

        assert n == 0, "Empty input should return 0 rows written"

        out_csv = gene_dir / "compounds_filtered.csv"
        assert out_csv.exists(), "CSV must be written even for empty input"

        with open(out_csv) as fh:
            content = fh.read()

        lines = [l for l in content.strip().splitlines() if l]
        assert len(lines) == 1, "Only header row expected for empty input"

        with open(out_csv) as fh:
            reader = csv.reader(fh)
            header = next(reader)

        assert header == _EXPECTED_COLUMNS

    def test_none_df_handled_gracefully(self, tmp_path: Path) -> None:
        """None candidates_df should not raise, should return 0."""
        from pipeline.stages.fragment.output import run_output

        settings, cache, console = _make_mocks(tmp_path / "results")
        gene_dir = tmp_path / "results" / "IGHMBP2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        n = run_output("IGHMBP2", None, settings, cache, force=True, console=console)
        assert n == 0


# ── Test 4: source column = "fragment_virtual_screen" ─────────────────────────

class TestSourceColumn:
    def test_source_column_value(self, tmp_path: Path) -> None:
        """All rows must have source='fragment_virtual_screen'."""
        from pipeline.stages.fragment.output import run_output

        settings, cache, console = _make_mocks(tmp_path / "results")
        gene_dir = tmp_path / "results" / "IGHMBP2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        candidates_df = _make_candidates_df(3)

        with patch(
            "pipeline.stages.fragment.output._run_admet_relaxed",
            return_value=(True, 0.60, {}),
        ):
            n = run_output(
                "IGHMBP2", candidates_df, settings, cache, force=True, console=console
            )

        assert n == 3

        out_df = pd.read_csv(gene_dir / "compounds_filtered.csv")
        assert "source" in out_df.columns, "'source' column must exist"
        assert (out_df["source"] == "fragment_virtual_screen").all(), (
            "All rows must have source='fragment_virtual_screen'"
        )

    def test_best_assay_type_column_value(self, tmp_path: Path) -> None:
        """All rows must have best_assay_type='fragment_screen_predicted'."""
        from pipeline.stages.fragment.output import run_output

        settings, cache, console = _make_mocks(tmp_path / "results")
        gene_dir = tmp_path / "results" / "IGHMBP2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        candidates_df = _make_candidates_df(2)

        with patch(
            "pipeline.stages.fragment.output._run_admet_relaxed",
            return_value=(True, 0.55, {}),
        ):
            run_output(
                "IGHMBP2", candidates_df, settings, cache, force=True, console=console
            )

        out_df = pd.read_csv(gene_dir / "compounds_filtered.csv")
        assert (out_df["best_assay_type"] == "fragment_screen_predicted").all()
