"""Unit tests for docking stage pure functions."""
import textwrap
from pathlib import Path

import pandas as pd
import pytest


# ── _check_vina_installed ──────────────────────────────────────────────────────

class TestCheckVinaInstalled:
    def test_raises_runtime_error_when_vina_missing(self, monkeypatch) -> None:
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "vina":
                raise ImportError("No module named 'vina'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)

        from pipeline.stages.dock.dock import _check_vina_installed
        with pytest.raises(RuntimeError) as exc_info:
            _check_vina_installed()

        msg = str(exc_info.value)
        assert "pip install" in msg
        assert "vina" in msg.lower()

    def test_error_message_includes_meeko(self, monkeypatch) -> None:
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "vina":
                raise ImportError
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)

        from pipeline.stages.dock.dock import _check_vina_installed
        with pytest.raises(RuntimeError) as exc_info:
            _check_vina_installed()
        assert "meeko" in str(exc_info.value)


# ── _define_box ────────────────────────────────────────────────────────────────

class TestDefineBox:
    def test_centroid_of_single_residue(self) -> None:
        from pipeline.stages.dock.dock import _define_box
        residues = [{"ca_x": 10.0, "ca_y": 20.0, "ca_z": 30.0}]
        center, size = _define_box(residues)
        assert center == pytest.approx([10.0, 20.0, 30.0])

    def test_centroid_of_two_residues(self) -> None:
        from pipeline.stages.dock.dock import _define_box
        residues = [
            {"ca_x": 0.0, "ca_y": 0.0, "ca_z": 0.0},
            {"ca_x": 4.0, "ca_y": 4.0, "ca_z": 4.0},
        ]
        center, size = _define_box(residues)
        assert center == pytest.approx([2.0, 2.0, 2.0])

    def test_box_size_equals_2_times_radius_plus_margin(self) -> None:
        from pipeline.stages.dock.dock import _define_box
        # Two residues 10 Å apart on X axis; centroid at x=5; radius = 5
        residues = [
            {"ca_x": 0.0, "ca_y": 0.0, "ca_z": 0.0},
            {"ca_x": 10.0, "ca_y": 0.0, "ca_z": 0.0},
        ]
        center, size = _define_box(residues, margin_A=3.0)
        expected_side = 2 * (5.0 + 3.0)
        assert size == pytest.approx([expected_side, expected_side, expected_side])

    def test_margin_affects_box_size(self) -> None:
        from pipeline.stages.dock.dock import _define_box
        residues = [
            {"ca_x": 0.0, "ca_y": 0.0, "ca_z": 0.0},
            {"ca_x": 10.0, "ca_y": 0.0, "ca_z": 0.0},
        ]
        _, size_3 = _define_box(residues, margin_A=3.0)
        _, size_5 = _define_box(residues, margin_A=5.0)
        assert size_5[0] > size_3[0]

    def test_symmetric_box(self) -> None:
        from pipeline.stages.dock.dock import _define_box
        residues = [{"ca_x": float(i), "ca_y": 0.0, "ca_z": 0.0} for i in range(5)]
        center, size = _define_box(residues)
        assert size[0] == size[1] == size[2]


# ── _parse_vina_remarks (via _run_vina score parsing) ─────────────────────────

def _make_pdbqt_content(*results: tuple[float, float, float]) -> str:
    """Build minimal PDBQT content with REMARK VINA RESULT lines."""
    lines = []
    for i, (aff, rmsd_lb, rmsd_ub) in enumerate(results, 1):
        lines.append(f"MODEL {i}")
        lines.append(f"REMARK VINA RESULT:    {aff:.3f}    {rmsd_lb:.3f}    {rmsd_ub:.3f}")
        lines.append("ENDMDL")
    return "\n".join(lines) + "\n"


class TestVinaRemarkParser:
    """Test the REMARK VINA RESULT parsing logic via a patched _run_vina."""

    def _parse_remarks(self, content: str) -> list[dict]:
        from pipeline.stages.dock.dock import _REMARK_PREFIX
        poses = []
        for line in content.splitlines():
            if line.startswith(_REMARK_PREFIX):
                parts = line[len(_REMARK_PREFIX):].split()
                if len(parts) >= 3:
                    poses.append({
                        "pose_rank": len(poses) + 1,
                        "affinity_kcal_mol": float(parts[0]),
                        "rmsd_lb": float(parts[1]),
                        "rmsd_ub": float(parts[2]),
                    })
        return poses

    def test_single_pose_parsed(self) -> None:
        content = _make_pdbqt_content((-8.5, 0.000, 0.000))
        poses = self._parse_remarks(content)
        assert len(poses) == 1
        assert poses[0]["affinity_kcal_mol"] == pytest.approx(-8.5)
        assert poses[0]["pose_rank"] == 1

    def test_three_poses_parsed_in_order(self) -> None:
        content = _make_pdbqt_content((-8.5, 0.0, 0.0), (-7.2, 1.5, 2.1), (-6.8, 2.0, 3.3))
        poses = self._parse_remarks(content)
        assert len(poses) == 3
        assert poses[0]["affinity_kcal_mol"] == pytest.approx(-8.5)
        assert poses[1]["affinity_kcal_mol"] == pytest.approx(-7.2)
        assert poses[2]["affinity_kcal_mol"] == pytest.approx(-6.8)

    def test_empty_content_returns_no_poses(self) -> None:
        poses = self._parse_remarks("")
        assert poses == []

    def test_non_remark_lines_ignored(self) -> None:
        content = "ATOM   1  CA  ALA A   1      0.000   0.000   0.000\n"
        poses = self._parse_remarks(content)
        assert poses == []


# ── _map_contacts ──────────────────────────────────────────────────────────────

class TestMapContacts:
    def _write_pose_pdbqt(self, path: Path, coords: list[tuple]) -> None:
        lines = ["MODEL 1\n"]
        for i, (x, y, z) in enumerate(coords, 1):
            lines.append(
                f"HETATM{i:5d}  C   LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 C\n"
            )
        lines.append("ENDMDL\n")
        path.write_text("".join(lines))

    def _write_comparison_csv(self, path: Path, rows: list[dict]) -> None:
        pd.DataFrame(rows).to_csv(path, index=False)

    def test_contact_within_cutoff_detected(self, tmp_path) -> None:
        from pipeline.stages.dock.dock import _map_contacts

        pose = tmp_path / "pose.pdbqt"
        csv = tmp_path / "comparison.csv"

        # Ligand atom at (1, 0, 0); residue Cα at (0, 0, 0) → dist = 1.0 Å < 4.0
        self._write_pose_pdbqt(pose, [(1.0, 0.0, 0.0)])
        self._write_comparison_csv(csv, [{
            "klifs_position": 45,
            "subpocket": "Gatekeeper",
            "vrk1_aa": "M",
            "egfr_aa": "T",
            "selectivity_candidate": True,
            "ca_x": 0.0, "ca_y": 0.0, "ca_z": 0.0,
        }])

        contacted = _map_contacts(pose, csv, cutoff_A=4.0)
        assert 45 in contacted

    def test_no_contact_beyond_cutoff(self, tmp_path) -> None:
        from pipeline.stages.dock.dock import _map_contacts

        pose = tmp_path / "pose.pdbqt"
        csv = tmp_path / "comparison.csv"

        # Ligand at (10, 0, 0); residue at (0, 0, 0) → dist = 10 Å > 4.0
        self._write_pose_pdbqt(pose, [(10.0, 0.0, 0.0)])
        self._write_comparison_csv(csv, [{
            "klifs_position": 45,
            "subpocket": "Gatekeeper",
            "vrk1_aa": "M",
            "egfr_aa": "T",
            "selectivity_candidate": True,
            "ca_x": 0.0, "ca_y": 0.0, "ca_z": 0.0,
        }])

        contacted = _map_contacts(pose, csv, cutoff_A=4.0)
        assert 45 not in contacted

    def test_non_candidate_positions_excluded(self, tmp_path) -> None:
        from pipeline.stages.dock.dock import _map_contacts

        pose = tmp_path / "pose.pdbqt"
        csv = tmp_path / "comparison.csv"

        self._write_pose_pdbqt(pose, [(0.5, 0.0, 0.0)])
        self._write_comparison_csv(csv, [{
            "klifs_position": 45,
            "subpocket": "Hinge",
            "vrk1_aa": "M",
            "egfr_aa": "M",
            "selectivity_candidate": False,  # not a candidate
            "ca_x": 0.0, "ca_y": 0.0, "ca_z": 0.0,
        }])

        contacted = _map_contacts(pose, csv, cutoff_A=4.0)
        assert 45 not in contacted

    def test_missing_ca_columns_returns_empty(self, tmp_path) -> None:
        from pipeline.stages.dock.dock import _map_contacts

        pose = tmp_path / "pose.pdbqt"
        csv = tmp_path / "comparison.csv"

        self._write_pose_pdbqt(pose, [(0.5, 0.0, 0.0)])
        pd.DataFrame([{"klifs_position": 45, "selectivity_candidate": True}]).to_csv(csv, index=False)

        contacted = _map_contacts(pose, csv, cutoff_A=4.0)
        assert contacted == []
