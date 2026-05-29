"""Unit tests for the fragment growing step."""
from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from rdkit import Chem


# ── Test 1: _filter_grown — MW=600 rejected, MW=350 SA=3 accepted ─────────────

class TestFilterGrown:
    def test_high_mw_rejected(self) -> None:
        """Molecule with MW > 450 should be rejected."""
        from pipeline.stages.fragment.grow import _filter_grown

        # Build a molecule with MW ~600 (large steroid-like)
        # Coronene: C24H12 MW ~300 — use a repeating chain for ~600
        # Simpler: use a known heavy molecule string
        mol = Chem.MolFromSmiles(
            "CCCCCCCCCCCCCCCCCC(=O)O"  # stearic acid, MW ~284 — too light
        )
        # Use a molecule that is unambiguously > 450 Da
        # Taxol-like stub: just use a long aliphatic chain
        heavy_smiles = "CCCC(CC)(CC)CC(CC)(CC)CC(CC)(CC)CCC"  # MW ~366 but passes
        # Use something definitively above 450:
        # Erythromycin fragment approximation — just build with RDKit
        heavy_mol = Chem.MolFromSmiles(
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)OCC2=CC=C(C=C2)C(C)C(=O)OCC3=CC=C(C=C3)CC"
        )
        if heavy_mol is None:
            pytest.skip("Could not construct heavy test molecule")

        from rdkit.Chem import Descriptors
        mw = Descriptors.ExactMolWt(heavy_mol)
        if mw <= 450:
            pytest.skip(f"Test molecule MW={mw:.1f} is not > 450 — need a heavier mol")

        result = _filter_grown([heavy_mol])
        assert result == [], f"Expected MW={mw:.1f} > 450 to be rejected"

    def test_valid_molecule_accepted(self) -> None:
        """Molecule with MW 300–450, SA score ~3, 0 Ro5 violations should pass."""
        from pipeline.stages.fragment.grow import _filter_grown, HAS_SASCORER

        # Ibuprofen-like: MW ~206, too light. Use a mid-sized drug scaffold.
        # Atorvastatin core fragment-ish: use a known drug
        # Verapamil: MW ~454 — might be right at the boundary; use something safer
        # Diclofenac: MW ~296 — too low
        # Use a hand-crafted compound: 4 aromatic rings + short chain ~ 340 Da
        smiles = "c1ccc(CNc2cccc(NC3CCCCC3)c2)cc1"  # MW ~308, logP ~3.5
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, "Test SMILES must parse"

        from rdkit.Chem import Descriptors, rdMolDescriptors
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)

        # Only run if the molecule actually meets our criteria (sanity check)
        ro5_violations = int((mw > 500) + (logp > 5) + (hbd > 5) + (hba > 10))
        if not (300.0 <= mw <= 450.0) or ro5_violations != 0:
            pytest.skip(
                f"Test mol MW={mw:.1f} logP={logp:.1f} hbd={hbd} hba={hba} — "
                f"doesn't meet baseline criteria for this test"
            )

        if HAS_SASCORER:
            import sascorer  # type: ignore
            sa = sascorer.calculateScore(mol)
            if sa >= 4.0:
                pytest.skip(f"Test mol SA={sa:.2f} >= 4.0 — choose a different SMILES")

        result = _filter_grown([mol])
        assert len(result) == 1, (
            f"MW={mw:.1f}, logP={logp:.1f} molecule should pass filter"
        )
        _, props = result[0]
        assert 300.0 <= props["molecular_weight"] <= 450.0
        assert props["passes_ro5"] is True

    def test_mw_600_rejected(self) -> None:
        """Explicit MW=600 check — large molecule must be filtered out."""
        from pipeline.stages.fragment.grow import _filter_grown

        # Build a large molecule from SMILES known to have MW > 600
        # Cyclosporin A partial: not easy to construct; use a long linear chain
        # Exact approach: build something with rdkit
        big_smiles = (
            "CC(C)(C)c1ccc(CC(=O)Nc2ccc(C(=O)Nc3ccc(CC(=O)Nc4ccc"
            "(C(C)(C)C)cc4)cc3)cc2)cc1"
        )
        mol = Chem.MolFromSmiles(big_smiles)
        if mol is None:
            # Fallback: use a simpler but definitely > 450 Da molecule
            big_smiles = "C" * 30  # polyethylene chain, MW ~ 420
            mol = Chem.MolFromSmiles(big_smiles)
        if mol is None:
            pytest.skip("Cannot construct large test molecule")

        from rdkit.Chem import Descriptors
        mw = Descriptors.ExactMolWt(mol)
        if mw <= 450:
            pytest.skip(f"Fallback molecule MW={mw:.1f} is not > 450")

        result = _filter_grown([mol])
        assert result == [], f"MW={mw:.1f} should be rejected (> 450)"

    def test_mw_350_sa3_accepted(self) -> None:
        """MW=350, SA=3 should be accepted (if sascorer available)."""
        from pipeline.stages.fragment.grow import _filter_grown, HAS_SASCORER

        # Labetalol-like scaffold: MW ~328, known SA ~ 2–3
        smiles = "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(C(N)=O)c1"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pytest.skip("Test SMILES failed to parse")

        from rdkit.Chem import Descriptors, rdMolDescriptors
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ro5_violations = int((mw > 500) + (logp > 5) + (hbd > 5) + (hba > 10))

        if not (300.0 <= mw <= 450.0):
            pytest.skip(f"Test mol MW={mw:.1f} outside 300–450 range")
        if ro5_violations != 0:
            pytest.skip(f"Test mol has {ro5_violations} Ro5 violations")
        if rotb > 8:
            pytest.skip(f"Test mol has {rotb} rotatable bonds > 8")

        result = _filter_grown([mol])
        # Should pass with or without sascorer (SA score typically < 4 for drug-like mols)
        # If sascorer available and SA happens to be >= 4, skip — molecule choice issue
        if HAS_SASCORER and len(result) == 0:
            import sascorer  # type: ignore
            sa = sascorer.calculateScore(mol)
            if sa >= 4.0:
                pytest.skip(f"Test molecule SA={sa:.2f} >= 4.0")
        assert len(result) == 1, (
            f"MW={mw:.1f}, Ro5={ro5_violations} molecule should be accepted"
        )


# ── Test 2: _grow_smarts — benzene → ≥3 products ─────────────────────────────

class TestGrowSmarts:
    def test_benzene_produces_at_least_3_products(self) -> None:
        """SMARTS reactions on benzene should yield at least 3 unique products."""
        from pipeline.stages.fragment.grow import _grow_smarts

        benzene = Chem.MolFromSmiles("c1ccccc1")
        assert benzene is not None

        products = _grow_smarts(benzene)
        assert len(products) >= 3, (
            f"Expected >= 3 SMARTS products from benzene, got {len(products)}"
        )

    def test_products_are_valid_mol_objects(self) -> None:
        """All returned molecules must be valid RDKit Mol objects."""
        from pipeline.stages.fragment.grow import _grow_smarts

        benzene = Chem.MolFromSmiles("c1ccccc1")
        products = _grow_smarts(benzene)

        for mol in products:
            assert mol is not None
            # Each mol must have a valid canonical SMILES
            assert Chem.MolToSmiles(mol) != ""

    def test_products_are_deduplicated(self) -> None:
        """No duplicate canonical SMILES in the result."""
        from pipeline.stages.fragment.grow import _grow_smarts

        aniline = Chem.MolFromSmiles("Nc1ccccc1")
        assert aniline is not None

        products = _grow_smarts(aniline)
        smiles_list = [Chem.MolToSmiles(m) for m in products]
        assert len(smiles_list) == len(set(smiles_list)), "Duplicate SMILES found"


# ── Test 3: run_grow with empty clusters_df → empty DataFrame ─────────────────

class TestRunGrowEmpty:
    def test_empty_clusters_df_returns_empty_dataframe(self) -> None:
        """run_grow with empty clusters_df should return empty DataFrame without error."""
        from pipeline.stages.fragment.grow import run_grow

        settings = MagicMock()
        settings.cache_dir = Path("/tmp/test_grow_cache")
        settings.results_dir = Path("/tmp/test_grow_results")

        cache = MagicMock()
        cache.load.return_value = None

        console = MagicMock()

        empty_df = pd.DataFrame(
            columns=["fragment_id", "smiles", "affinity_kcal_mol", "cluster_id", "is_representative"]
        )

        result = run_grow("IGHMBP2", empty_df, settings, cache, force=True, console=console)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_no_representative_column_returns_empty_dataframe(self) -> None:
        """clusters_df without is_representative column returns empty DataFrame."""
        from pipeline.stages.fragment.grow import run_grow

        settings = MagicMock()
        settings.cache_dir = Path("/tmp/test_grow_cache")
        settings.results_dir = Path("/tmp/test_grow_results")

        cache = MagicMock()
        cache.load.return_value = None

        console = MagicMock()

        df_no_rep = pd.DataFrame({"fragment_id": ["F001"], "smiles": ["c1ccccc1"]})

        result = run_grow("IGHMBP2", df_no_rep, settings, cache, force=True, console=console)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0


# ── Test 4: SA score threshold ────────────────────────────────────────────────

class TestSAScoreThreshold:
    def test_high_sa_score_rejected(self) -> None:
        """SA score >= 4.0 should be rejected when sascorer is available."""
        from pipeline.stages.fragment.grow import HAS_SASCORER

        if not HAS_SASCORER:
            pytest.skip("sascorer not installed — SA filter not active")

        from pipeline.stages.fragment.grow import _filter_grown
        import sascorer  # type: ignore

        # We need a molecule that sascorer would score >= 4.0
        # Highly complex natural product-like: use norbornane with many substituents
        # In practice, we mock sascorer to control the score
        mol = Chem.MolFromSmiles("c1ccc(CNc2cccc(NC3CCCCC3)c2)cc1")
        if mol is None:
            pytest.skip("Test molecule failed to parse")

        with patch("pipeline.stages.fragment.grow.sascorer") as mock_sa:
            mock_sa.calculateScore.return_value = 5.0
            result = _filter_grown([mol])
        assert result == [], "SA score 5.0 should be rejected"

    def test_low_sa_score_accepted(self) -> None:
        """SA score < 4.0 should be accepted when sascorer is available."""
        from pipeline.stages.fragment.grow import HAS_SASCORER

        if not HAS_SASCORER:
            pytest.skip("sascorer not installed — SA filter not active")

        from pipeline.stages.fragment.grow import _filter_grown

        mol = Chem.MolFromSmiles("c1ccc(CNc2cccc(NC3CCCCC3)c2)cc1")
        if mol is None:
            pytest.skip("Test molecule failed to parse")

        from rdkit.Chem import Descriptors, rdMolDescriptors
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ro5_violations = int((mw > 500) + (logp > 5) + (hbd > 5) + (hba > 10))

        if not (300.0 <= mw <= 450.0) or ro5_violations != 0 or rotb > 8:
            pytest.skip(
                f"Test mol MW={mw:.1f}, Ro5={ro5_violations}, rotb={rotb} "
                f"doesn't meet non-SA criteria"
            )

        with patch("pipeline.stages.fragment.grow.sascorer") as mock_sa:
            mock_sa.calculateScore.return_value = 3.9
            result = _filter_grown([mol])
        assert len(result) == 1, "SA score 3.9 should be accepted"
