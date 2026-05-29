"""Unit tests for fragment screening stage (pipeline/stages/fragment/screen.py)."""
from __future__ import annotations

import io
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch, call

import pandas as pd
import pytest

from pipeline.config import Settings


# ── Helpers ────────────────────────────────────────────────────────────────────

def _make_settings(tmp_path: Path) -> Settings:
    settings = Settings(data_dir=tmp_path / "data")
    settings.cache_dir.mkdir(parents=True, exist_ok=True)
    settings.results_dir.mkdir(parents=True, exist_ok=True)
    return settings


def _make_pocket() -> dict:
    return {
        "centroid_x": 10.0,
        "centroid_y": 20.0,
        "centroid_z": 30.0,
        "box_size_A": 20.0,
    }


def _make_cache(load_return=None) -> MagicMock:
    cache = MagicMock()
    cache.load.return_value = load_return
    return cache


def _make_console() -> MagicMock:
    return MagicMock()


def _write_library(path: Path, fragments: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{smiles}\t{fid}" for smiles, fid in fragments]
    path.write_text("\n".join(lines) + "\n")


def _make_pdb(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")


# ── Test 1: Cache hit — _run_vina not called ───────────────────────────────────

class TestCacheHit:
    def test_cache_hit_skips_vina(self, tmp_path: Path) -> None:
        """Given a cached result for a fragment, _run_vina must NOT be called."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        # Create a fake receptor PDB
        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        # Library with one fragment
        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, [("c1ccccc1", "frag_001")])

        # Cache already has a result for frag_001
        cached_result = {
            "fragment_id": "frag_001",
            "smiles": "c1ccccc1",
            "affinity_kcal_mol": -5.5,
            "n_poses": 3,
        }
        cache = _make_cache(load_return=cached_result)

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor") as mock_prep_rec,
            patch("pipeline.stages.fragment.screen._prepare_ligand") as mock_prep_lig,
            patch("pipeline.stages.fragment.screen._run_vina") as mock_vina,
        ):
            mock_prep_rec.return_value = tmp_path / "receptor.pdbqt"

            from pipeline.stages.fragment.screen import run_screen
            result = run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=False,
                console=_make_console(),
            )

        mock_vina.assert_not_called()
        mock_prep_lig.assert_not_called()
        assert len(result) == 1
        assert result.iloc[0]["affinity_kcal_mol"] == pytest.approx(-5.5)

    def test_force_flag_bypasses_cache(self, tmp_path: Path) -> None:
        """With force=True, cache is not consulted and vina is called."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, [("c1ccccc1", "frag_001")])

        cache = _make_cache(load_return={"fragment_id": "frag_001", "smiles": "c1ccccc1",
                                         "affinity_kcal_mol": -5.5, "n_poses": 3})

        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", return_value=ligand_pdbqt),
            patch("pipeline.stages.fragment.screen._run_vina", return_value=[{"affinity_kcal_mol": -6.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]) as mock_vina,
        ):
            from pipeline.stages.fragment.screen import run_screen
            run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=_make_console(),
            )

        mock_vina.assert_called_once()
        # cache.load should NOT have been called when force=True
        cache.load.assert_not_called()


# ── Test 2: Fragment failure — screen continues ────────────────────────────────

class TestFragmentFailure:
    def test_failing_fragment_is_skipped(self, tmp_path: Path) -> None:
        """_prepare_ligand or _run_vina raising Exception skips that fragment; others processed."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, [
            ("c1ccccc1", "frag_001"),
            ("INVALID_SMILES_THAT_BREAKS", "frag_002"),
            ("CCO", "frag_003"),
        ])

        cache = _make_cache(load_return=None)

        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"

        def fake_prepare_ligand(smiles, fid, cache_dir):
            if fid == "frag_002":
                raise RuntimeError("meeko conversion failed")
            return ligand_pdbqt

        def fake_run_vina(rec, lig, center, box, exh, out):
            return [{"affinity_kcal_mol": -5.0 if "frag_001" in str(lig.parent) else -4.0,
                     "rmsd_lb": 0.0, "rmsd_ub": 0.0}]

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", side_effect=fake_prepare_ligand),
            patch("pipeline.stages.fragment.screen._run_vina") as mock_vina,
        ):
            mock_vina.return_value = [{"affinity_kcal_mol": -5.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]
            from pipeline.stages.fragment.screen import run_screen
            result = run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=_make_console(),
            )

        # frag_001 and frag_003 succeed; frag_002 is skipped
        assert len(result) == 2
        fragment_ids = set(result["fragment_id"].tolist())
        assert "frag_001" in fragment_ids
        assert "frag_003" in fragment_ids
        assert "frag_002" not in fragment_ids

    def test_vina_failure_skips_fragment(self, tmp_path: Path) -> None:
        """_run_vina raising an exception skips the fragment without stopping the screen."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, [("c1ccccc1", "frag_A"), ("CCN", "frag_B")])

        cache = _make_cache(load_return=None)
        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"

        call_count = {"n": 0}

        def fake_vina(rec, lig, center, box, exh, out):
            call_count["n"] += 1
            if call_count["n"] == 1:
                raise RuntimeError("Vina crashed")
            return [{"affinity_kcal_mol": -6.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", return_value=ligand_pdbqt),
            patch("pipeline.stages.fragment.screen._run_vina", side_effect=fake_vina),
        ):
            from pipeline.stages.fragment.screen import run_screen
            result = run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=_make_console(),
            )

        assert len(result) == 1
        assert result.iloc[0]["fragment_id"] == "frag_B"


# ── Test 3: Top-N selection ────────────────────────────────────────────────────

class TestTopNSelection:
    def test_returns_top_n_sorted_by_affinity(self, tmp_path: Path) -> None:
        """Given 100 docked fragments, returns the top 50 sorted by affinity (most negative first)."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        # 100 fragments with affinities from -1.0 to -100.0 (descending in library order)
        frags = [(f"CCO{i}", f"frag_{i:03d}") for i in range(100)]
        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, frags)

        cache = _make_cache(load_return=None)
        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"

        frag_affinities = {f"frag_{i:03d}": -float(i + 1) for i in range(100)}

        def fake_vina(rec, lig, center, box, exh, out):
            # out is fragment_cache_dir / f"poses_{fid}.pdbqt"
            fid = out.stem[len("poses_"):]
            aff = frag_affinities.get(fid, -1.0)
            return [{"affinity_kcal_mol": aff, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", return_value=ligand_pdbqt),
            patch("pipeline.stages.fragment.screen._run_vina", side_effect=fake_vina),
        ):
            from pipeline.stages.fragment.screen import run_screen
            result = run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=_make_console(),
            )

        assert len(result) == 50
        # Most negative affinity should be first (frag_099 → -100.0)
        assert result.iloc[0]["affinity_kcal_mol"] == pytest.approx(-100.0)
        assert result.iloc[-1]["affinity_kcal_mol"] == pytest.approx(-51.0)
        # Verify sorted order
        affs = result["affinity_kcal_mol"].tolist()
        assert affs == sorted(affs)

    def test_fewer_than_top_n_results(self, tmp_path: Path) -> None:
        """If fewer fragments succeed than top_n, returns all successful results."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, [("c1ccccc1", "frag_001"), ("CCO", "frag_002")])

        cache = _make_cache(load_return=None)
        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", return_value=ligand_pdbqt),
            patch("pipeline.stages.fragment.screen._run_vina",
                  return_value=[{"affinity_kcal_mol": -5.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]),
        ):
            from pipeline.stages.fragment.screen import run_screen
            result = run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=_make_console(),
            )

        assert len(result) == 2


# ── Test 4: Progress reporting ─────────────────────────────────────────────────

class TestProgressReporting:
    def test_progress_logged_at_100_boundary(self, tmp_path: Path) -> None:
        """Console should receive a progress message at the 100-fragment boundary."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        # 100 fragments to trigger the boundary exactly at i=100
        frags = [(f"CCO{i}", f"frag_{i:03d}") for i in range(100)]
        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, frags)

        cache = _make_cache(load_return=None)
        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"
        console = _make_console()

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", return_value=ligand_pdbqt),
            patch("pipeline.stages.fragment.screen._run_vina",
                  return_value=[{"affinity_kcal_mol": -5.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]),
        ):
            from pipeline.stages.fragment.screen import run_screen
            run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=console,
            )

        # Inspect all print calls for a progress line containing "100/100"
        all_calls = [str(c) for c in console.print.call_args_list]
        progress_calls = [c for c in all_calls if "100/100" in c]
        assert len(progress_calls) >= 1, (
            "Expected at least one progress message containing '100/100'. "
            f"Got console calls: {all_calls}"
        )

    def test_progress_message_contains_failed_count(self, tmp_path: Path) -> None:
        """Progress message at 100-fragment boundary includes the failed count."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"

        pdb_path = settings.cache_dir / gene / "structures" / "receptor.pdb"
        _make_pdb(pdb_path)

        frags = [(f"CCO{i}", f"frag_{i:03d}") for i in range(100)]
        lib_path = tmp_path / "fragments.smi"
        _write_library(lib_path, frags)

        cache = _make_cache(load_return=None)
        receptor_pdbqt = tmp_path / "rec.pdbqt"
        ligand_pdbqt = tmp_path / "lig.pdbqt"
        console = _make_console()

        call_count = {"n": 0}

        def fake_prepare_ligand(smiles, fid, cache_dir):
            call_count["n"] += 1
            # Fail frag_050 to introduce 1 failure
            if fid == "frag_050":
                raise RuntimeError("fail")
            return ligand_pdbqt

        with (
            patch("pipeline.stages.fragment.screen._prepare_receptor", return_value=receptor_pdbqt),
            patch("pipeline.stages.fragment.screen._prepare_ligand", side_effect=fake_prepare_ligand),
            patch("pipeline.stages.fragment.screen._run_vina",
                  return_value=[{"affinity_kcal_mol": -5.0, "rmsd_lb": 0.0, "rmsd_ub": 0.0}]),
        ):
            from pipeline.stages.fragment.screen import run_screen
            run_screen(
                gene_symbol=gene,
                library_path=lib_path,
                pocket=_make_pocket(),
                top_n=50,
                exhaustiveness=4,
                settings=settings,
                cache=cache,
                force=True,
                console=console,
            )

        all_calls = [str(c) for c in console.print.call_args_list]
        progress_calls = [c for c in all_calls if "100/100" in c]
        assert len(progress_calls) >= 1
        # At least one progress call should mention the failure count
        assert any("failed" in c for c in progress_calls)
