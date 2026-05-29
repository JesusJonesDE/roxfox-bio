"""Unit tests for fragment clustering stage (pipeline/stages/fragment/cluster.py)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from pipeline.config import Settings


# ── Helpers ────────────────────────────────────────────────────────────────────

def _make_settings(tmp_path: Path) -> Settings:
    settings = Settings(data_dir=tmp_path / "data")
    settings.cache_dir.mkdir(parents=True, exist_ok=True)
    settings.results_dir.mkdir(parents=True, exist_ok=True)
    return settings


def _make_cache(load_return=None) -> MagicMock:
    cache = MagicMock()
    cache.load.return_value = load_return
    return cache


def _make_console() -> MagicMock:
    return MagicMock()


def _make_hits_df(smiles_list: list[str], affinities: list[float] | None = None) -> pd.DataFrame:
    """Build a hits_df with the expected columns, sorted by affinity ascending."""
    n = len(smiles_list)
    if affinities is None:
        affinities = [-float(i + 1) for i in range(n)]  # -1.0, -2.0, …, -n (ascending)
    rows = [
        {
            "fragment_id": f"frag_{i:03d}",
            "smiles": smi,
            "affinity_kcal_mol": aff,
            "n_poses": 3,
        }
        for i, (smi, aff) in enumerate(zip(smiles_list, affinities))
    ]
    df = pd.DataFrame(rows)
    return df.sort_values("affinity_kcal_mol").reset_index(drop=True)


# ── Test 1: Identical SMILES → 1 cluster ──────────────────────────────────────

class TestIdenticalSmiles:
    def test_ten_identical_smiles_form_one_cluster(self, tmp_path: Path) -> None:
        """10 identical SMILES must all land in the same cluster."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        smiles = ["c1ccccc1"] * 10
        hits_df = _make_hits_df(smiles)

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        assert "cluster_id" in result.columns
        n_clusters = result["cluster_id"].nunique()
        assert n_clusters == 1, f"Expected 1 cluster for identical SMILES, got {n_clusters}"

    def test_single_representative_per_cluster(self, tmp_path: Path) -> None:
        """With 10 identical SMILES in one cluster, exactly 1 is_representative=True."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        smiles = ["c1ccccc1"] * 10
        hits_df = _make_hits_df(smiles)

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        assert result["is_representative"].sum() == 1


# ── Test 2: Distinct SMILES → many clusters ────────────────────────────────────

class TestDistinctSmiles:
    # Structurally distinct fragments: benzene, indole, morpholine,
    # imidazole, pyridine — all have Tanimoto << 0.6 from each other
    _DISTINCT = [
        "c1ccccc1",          # benzene
        "c1ccc2[nH]ccc2c1",  # indole
        "C1COCCN1",          # morpholine
        "c1cnc[nH]1",        # imidazole
        "c1ccncc1",          # pyridine
    ]

    def test_five_distinct_smiles_form_five_clusters(self, tmp_path: Path) -> None:
        """5 structurally distinct SMILES should each land in their own cluster."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        hits_df = _make_hits_df(self._DISTINCT)

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        n_clusters = result["cluster_id"].nunique()
        assert n_clusters == 5, (
            f"Expected 5 clusters for 5 distinct SMILES, got {n_clusters}"
        )

    def test_all_distinct_fragments_are_representatives(self, tmp_path: Path) -> None:
        """When every fragment is in its own cluster, all are representatives."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        hits_df = _make_hits_df(self._DISTINCT)

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        assert result["is_representative"].all()


# ── Test 3: Representative has best affinity ───────────────────────────────────

class TestRepresentativeAffinity:
    def test_representative_has_most_negative_affinity(self, tmp_path: Path) -> None:
        """Within each cluster, the representative must have the most negative affinity."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        # 6 benzene variants (very similar → 1 cluster) with varied affinities
        # Affinities: assigned so the most negative is NOT the first in the list
        smiles_list = ["c1ccccc1"] * 6
        # sorted ascending (most negative first) — as run_screen would deliver them
        affinities = [-9.5, -8.0, -7.5, -7.0, -6.0, -5.0]
        hits_df = _make_hits_df(smiles_list, affinities=affinities)

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        for cid in result["cluster_id"].unique():
            cluster_members = result[result["cluster_id"] == cid]
            rep = cluster_members[cluster_members["is_representative"]]
            assert len(rep) == 1, f"Cluster {cid} has {len(rep)} representatives"
            best_affinity = cluster_members["affinity_kcal_mol"].min()
            assert rep.iloc[0]["affinity_kcal_mol"] == pytest.approx(best_affinity), (
                f"Representative in cluster {cid} has affinity "
                f"{rep.iloc[0]['affinity_kcal_mol']}, but best is {best_affinity}"
            )


# ── Test 4: Single-fragment input ─────────────────────────────────────────────

class TestSingleFragment:
    def test_single_fragment_is_one_cluster_and_representative(self, tmp_path: Path) -> None:
        """A single-fragment hits_df produces 1 cluster with is_representative=True."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        hits_df = _make_hits_df(["c1ccccc1"], affinities=[-7.0])

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        assert len(result) == 1
        assert result["cluster_id"].nunique() == 1
        assert bool(result.iloc[0]["is_representative"]) is True

    def test_single_fragment_csv_written(self, tmp_path: Path) -> None:
        """Single-fragment input produces a valid fragment_clusters.csv."""
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)
        cache = _make_cache(load_return=None)

        hits_df = _make_hits_df(["CCO"], affinities=[-6.5])

        from pipeline.stages.fragment.cluster import run_cluster
        run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=True,
            console=_make_console(),
        )

        csv_path = settings.results_dir / gene / "fragment_clusters.csv"
        assert csv_path.exists(), "fragment_clusters.csv was not written"
        df = pd.read_csv(csv_path)
        assert "cluster_id" in df.columns
        assert "is_representative" in df.columns
        assert len(df) == 1


# ── Test 5: Cache hit returns cached DataFrame ─────────────────────────────────

class TestCacheHit:
    def test_cache_hit_returns_cached_dataframe(self, tmp_path: Path) -> None:
        """If cached data exists and force=False, returns cached DataFrame directly."""
        settings = _make_settings(tmp_path)
        gene = "IGHMBP2"
        (settings.results_dir / gene).mkdir(parents=True, exist_ok=True)

        cached_records = [
            {"fragment_id": "frag_001", "smiles": "c1ccccc1",
             "affinity_kcal_mol": -7.0, "cluster_id": 0, "is_representative": True}
        ]
        cache = _make_cache(load_return=cached_records)

        hits_df = _make_hits_df(["c1ccccc1"])

        from pipeline.stages.fragment.cluster import run_cluster
        result = run_cluster(
            gene_symbol=gene,
            hits_df=hits_df,
            settings=settings,
            cache=cache,
            force=False,
            console=_make_console(),
        )

        assert len(result) == 1
        assert result.iloc[0]["cluster_id"] == 0
        # Cache.save should NOT be called — we returned early
        cache.save.assert_not_called()
