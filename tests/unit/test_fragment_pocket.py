"""Unit tests for pipeline/stages/fragment/pocket.py"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


# ── Helpers to build mock fpocket output ──────────────────────────────────────

def _make_info_txt(pockets: list[dict]) -> str:
    """Build a minimal fpocket *_info.txt string from a list of pocket dicts.

    Each dict should have keys matching _COLUMNS:
        Pocket, Score, Drugg_score, Volume, (others optional → default 0.0)
    """
    header = (
        "Pocket Score Drugg_score Total_SASA Polar_SASA Apolar_SASA Volume "
        "Mean_local_hydrophobic_dens Mean_acrophobicity Mean_pocketness "
        "Mean_accept_angle Mean_donor_angle Meam_water_dist Flexibility"
    )
    lines = [f"# fpocket info", header]
    defaults = {k: "0.0" for k in header.split()}
    for p in pockets:
        row = {**defaults, **{k: str(v) for k, v in p.items()}}
        lines.append(
            f"{row['Pocket']} {row['Score']} {row['Drugg_score']} "
            f"{row['Total_SASA']} {row['Polar_SASA']} {row['Apolar_SASA']} "
            f"{row['Volume']} {row['Mean_local_hydrophobic_dens']} "
            f"{row['Mean_acrophobicity']} {row['Mean_pocketness']} "
            f"{row['Mean_accept_angle']} {row['Mean_donor_angle']} "
            f"{row['Meam_water_dist']} {row['Flexibility']}"
        )
    return "\n".join(lines) + "\n"


def _make_pocket_pdb(atoms: list[tuple[float, float, float]]) -> str:
    """Build a minimal pocket_atm.pdb string from a list of (x, y, z) tuples."""
    lines = []
    for i, (x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"ATOM  {i:4d}  C   LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


# ── Fixtures ───────────────────────────────────────────────────────────────────

@pytest.fixture()
def settings(tmp_path):
    from pipeline.config import Settings
    s = Settings(data_dir=tmp_path)
    return s


@pytest.fixture()
def cache(settings):
    from pipeline.cache import CacheManager
    return CacheManager(settings)


@pytest.fixture()
def console():
    from rich.console import Console
    return Console(quiet=True)


# ── Test 1: Correct pocket selected from mock fpocket output ──────────────────

class TestPocketSelection:
    def test_top_druggability_selected(self, tmp_path, settings, cache, console):
        """Pocket with highest Drugg_score (and Volume > 200) is selected."""
        from pipeline.stages.fragment.pocket import (
            _parse_info_txt,
            _select_pocket,
        )

        # Pocket 2 has the highest druggability score and sufficient volume
        pockets = [
            {"Pocket": 1, "Score": 0.8, "Drugg_score": 0.6, "Volume": 210.0},
            {"Pocket": 2, "Score": 0.9, "Drugg_score": 0.85, "Volume": 250.0},
            {"Pocket": 3, "Score": 0.5, "Drugg_score": 0.3, "Volume": 300.0},
        ]
        info_txt = _make_info_txt(pockets)
        info_path = tmp_path / "test_info.txt"
        info_path.write_text(info_txt)

        rows = _parse_info_txt(info_path)
        selected = _select_pocket(rows, console)

        assert int(selected["Pocket"]) == 2
        assert float(selected["Drugg_score"]) == pytest.approx(0.85)

    def test_highest_drugg_among_volume_filtered(self, tmp_path, settings, cache, console):
        """Pocket 1 (Volume=500) beats pocket 3 (Volume=201) by drugg score;
        pocket 2 (Volume=180) is excluded by volume filter."""
        from pipeline.stages.fragment.pocket import _parse_info_txt, _select_pocket

        pockets = [
            {"Pocket": 1, "Score": 0.9, "Drugg_score": 0.9, "Volume": 500.0},
            {"Pocket": 2, "Score": 0.95, "Drugg_score": 0.95, "Volume": 180.0},  # vol < 200
            {"Pocket": 3, "Score": 0.7, "Drugg_score": 0.7, "Volume": 201.0},
        ]
        info_path = tmp_path / "info.txt"
        info_path.write_text(_make_info_txt(pockets))

        rows = _parse_info_txt(info_path)
        selected = _select_pocket(rows, console)

        assert int(selected["Pocket"]) == 1


# ── Test 2: Volume < 200 Å³ triggers threshold relaxation ────────────────────

class TestVolumeThresholdRelaxation:
    def test_relaxes_to_150_with_warning(self, tmp_path, console):
        """When all pockets have Volume < 200, threshold relaxes to 150 with warning."""
        from pipeline.stages.fragment.pocket import _parse_info_txt, _select_pocket

        # All pockets have volume between 150 and 200
        pockets = [
            {"Pocket": 1, "Score": 0.7, "Drugg_score": 0.5, "Volume": 160.0},
            {"Pocket": 2, "Score": 0.8, "Drugg_score": 0.7, "Volume": 175.0},
            {"Pocket": 3, "Score": 0.6, "Drugg_score": 0.4, "Volume": 155.0},
        ]
        info_path = tmp_path / "info.txt"
        info_path.write_text(_make_info_txt(pockets))

        rows = _parse_info_txt(info_path)

        # Capture console output to confirm warning is emitted
        from rich.console import Console
        from io import StringIO

        buf = StringIO()
        warn_console = Console(file=buf, highlight=False)

        selected = _select_pocket(rows, warn_console)

        # Pocket 2 has highest Drugg_score among the relaxed candidates
        assert int(selected["Pocket"]) == 2
        output = buf.getvalue()
        assert "relax" in output.lower() or "150" in output

    def test_no_pockets_raises(self, tmp_path, console):
        """If truly no pockets pass even the relaxed threshold, RuntimeError is raised."""
        from pipeline.stages.fragment.pocket import _parse_info_txt, _select_pocket

        # All pockets have volume < 150
        pockets = [
            {"Pocket": 1, "Score": 0.9, "Drugg_score": 0.9, "Volume": 100.0},
        ]
        info_path = tmp_path / "info.txt"
        info_path.write_text(_make_info_txt(pockets))

        rows = _parse_info_txt(info_path)
        with pytest.raises(RuntimeError, match="No pockets found"):
            _select_pocket(rows, console)


# ── Test 3: Centroid computed correctly ───────────────────────────────────────

class TestCentroidExtraction:
    def test_centroid_of_four_atoms(self, tmp_path):
        """Mean X/Y/Z of 4 atoms is computed correctly."""
        from pipeline.stages.fragment.pocket import _extract_centroid

        atoms = [
            (0.0, 0.0, 0.0),
            (4.0, 0.0, 0.0),
            (0.0, 4.0, 0.0),
            (0.0, 0.0, 4.0),
        ]
        pdb_path = tmp_path / "pocket1_atm.pdb"
        pdb_path.write_text(_make_pocket_pdb(atoms))

        cx, cy, cz = _extract_centroid(pdb_path)

        assert cx == pytest.approx(1.0)
        assert cy == pytest.approx(1.0)
        assert cz == pytest.approx(1.0)

    def test_centroid_single_atom(self, tmp_path):
        """Centroid of one atom equals its coordinates."""
        from pipeline.stages.fragment.pocket import _extract_centroid

        pdb_path = tmp_path / "pocket1_atm.pdb"
        pdb_path.write_text(_make_pocket_pdb([(7.5, 8.5, 9.5)]))

        cx, cy, cz = _extract_centroid(pdb_path)
        assert cx == pytest.approx(7.5)
        assert cy == pytest.approx(8.5)
        assert cz == pytest.approx(9.5)

    def test_centroid_symmetric_cube(self, tmp_path):
        """8 corners of a unit cube → centroid at (0.5, 0.5, 0.5)."""
        from pipeline.stages.fragment.pocket import _extract_centroid

        atoms = [
            (x, y, z)
            for x in (0.0, 1.0)
            for y in (0.0, 1.0)
            for z in (0.0, 1.0)
        ]
        pdb_path = tmp_path / "pocket_atm.pdb"
        pdb_path.write_text(_make_pocket_pdb(atoms))

        cx, cy, cz = _extract_centroid(pdb_path)
        assert cx == pytest.approx(0.5)
        assert cy == pytest.approx(0.5)
        assert cz == pytest.approx(0.5)


# ── Test 4: Integration — run_pocket with mocked subprocess ──────────────────

class TestRunPocket:
    def _setup_mock_fpocket_output(
        self,
        pdb_path: Path,
        pockets: list[dict],
        pocket_atoms: list[tuple[float, float, float]],
        best_pocket_id: int,
    ) -> None:
        """Create the fpocket output directory structure that run_pocket expects."""
        stem = pdb_path.stem
        out_dir = pdb_path.parent / f"{stem}_out"
        pockets_dir = out_dir / "pockets"
        pockets_dir.mkdir(parents=True)

        # Write _info.txt
        info_path = out_dir / f"{stem}_info.txt"
        info_path.write_text(_make_info_txt(pockets))

        # Write pocket atom PDB for the expected best pocket
        pocket_pdb = pockets_dir / f"pocket{best_pocket_id}_atm.pdb"
        pocket_pdb.write_text(_make_pocket_pdb(pocket_atoms))

    def test_run_pocket_selects_correct_pocket(self, tmp_path, settings, cache, console):
        """run_pocket returns the pocket dict with expected values."""
        from pipeline.stages.fragment.pocket import run_pocket

        gene = "IGHMBP2"
        struct_dir = settings.cache_dir / gene / "structures"
        struct_dir.mkdir(parents=True)
        pdb_path = struct_dir / "IGHMBP2_AF2.pdb"
        pdb_path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")

        pockets = [
            {"Pocket": 1, "Score": 0.7, "Drugg_score": 0.6, "Volume": 210.0},
            {"Pocket": 2, "Score": 0.9, "Drugg_score": 0.85, "Volume": 300.0},
        ]
        atoms = [(1.0, 2.0, 3.0), (3.0, 4.0, 5.0)]
        self._setup_mock_fpocket_output(pdb_path, pockets, atoms, best_pocket_id=2)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            result = run_pocket(gene, settings, cache, force=True, console=console)

        assert result["pocket_id"] == 2
        assert result["druggability_score"] == pytest.approx(0.85)
        assert result["volume_A3"] == pytest.approx(300.0)
        assert result["centroid_x"] == pytest.approx(2.0)
        assert result["centroid_y"] == pytest.approx(3.0)
        assert result["centroid_z"] == pytest.approx(4.0)
        assert result["box_size_A"] == pytest.approx(20.0)
        assert result["plddt_mean"] is None

    def test_run_pocket_writes_json(self, tmp_path, settings, cache, console):
        """run_pocket writes pocket_analysis.json to results dir."""
        from pipeline.stages.fragment.pocket import run_pocket

        gene = "IGHMBP2"
        struct_dir = settings.cache_dir / gene / "structures"
        struct_dir.mkdir(parents=True)
        pdb_path = struct_dir / "IGHMBP2_AF2.pdb"
        pdb_path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")

        pockets = [{"Pocket": 1, "Score": 0.8, "Drugg_score": 0.75, "Volume": 220.0}]
        atoms = [(10.0, 20.0, 30.0)]
        self._setup_mock_fpocket_output(pdb_path, pockets, atoms, best_pocket_id=1)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            run_pocket(gene, settings, cache, force=True, console=console)

        json_path = settings.results_dir / gene / "pocket_analysis.json"
        assert json_path.exists()
        data = json.loads(json_path.read_text())
        assert data["pocket_id"] == 1

    def test_fpocket_not_found_raises(self, tmp_path, settings, cache, console):
        """FileNotFoundError from subprocess is converted to a RuntimeError."""
        from pipeline.stages.fragment.pocket import run_pocket

        gene = "IGHMBP2"
        struct_dir = settings.cache_dir / gene / "structures"
        struct_dir.mkdir(parents=True)
        pdb_path = struct_dir / "IGHMBP2_AF2.pdb"
        pdb_path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")

        with patch("subprocess.run", side_effect=FileNotFoundError("fpocket")):
            with pytest.raises(RuntimeError, match="fpocket not installed"):
                run_pocket(gene, settings, cache, force=True, console=console)
