"""Unit tests for pipeline/stages/fragment/library.py"""
from __future__ import annotations

import textwrap
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ── Fixtures ───────────────────────────────────────────────────────────────────

@pytest.fixture()
def settings(tmp_path):
    from pipeline.config import Settings
    return Settings(data_dir=tmp_path)


@pytest.fixture()
def cache(settings):
    from pipeline.cache import CacheManager
    return CacheManager(settings)


@pytest.fixture()
def console():
    from rich.console import Console
    return Console(quiet=True)


# ── Helper: write a minimal fallback library ───────────────────────────────────

def _write_fallback(settings, smiles_list: list[str]) -> Path:
    """Write a fragments_fallback.smi to the expected location."""
    fallback_dir = settings.cache_dir / "shared" / "fragment_library"
    fallback_dir.mkdir(parents=True, exist_ok=True)
    fallback_path = fallback_dir / "fragments_fallback.smi"
    lines = [f"{smi}\tFALLBACK-{i:05d}" for i, smi in enumerate(smiles_list, 1)]
    fallback_path.write_text("\n".join(lines) + "\n")
    return fallback_path


# ── Test 1: Ro3 filter ────────────────────────────────────────────────────────

class TestRo3Filter:
    def test_mw_260_rejected(self):
        """A molecule with MW > 250 is rejected by the Ro3 filter."""
        from rdkit import Chem
        from pipeline.stages.fragment.library import _passes_ro3

        # Naphthalene (C10H8, MW ≈ 128) is fine as a base; we need a heavier mol.
        # Tryptophan: C11H12N2O2, MW ≈ 204 — passes. Use a bigger one.
        # 4-Aminobiphenyl: C12H11N, MW ≈ 169 — passes.
        # Use a MW=260-ish molecule: e.g. Stearic acid sub-fragment or just
        # construct one via SMILES that clearly fails MW test.
        # Dodecanedioic acid: HOOC-(CH2)10-COOH, MW ~ 230 — borderline.
        # Let's use: C(CCCC)(CCCC)(CCCC)C — neopentane-like, MW ~ 226, logP > 3 — may fail logP.
        # Safest: a known fragment with MW ~260 using rdkit weight.
        # We'll just compute directly.
        from rdkit.Chem import Descriptors

        # Construct a molecule with MW just over 250
        smi = "CCCCCCCCCCCCC"  # tridecane, MW ≈ 184 — too light
        # Use: CCC(CC)c1ccc(CC)cc1 — 4-sec-butyl-4'-ethylbiphenyl-like? Let's be explicit.
        # Caffeic acid phenethyl ester CAPE: MW 284 — over 250
        smi_heavy = "O=C(OCCc1ccccc1)/C=C/c1ccc(O)c(O)c1"  # CAPE, MW 284
        mol = Chem.MolFromSmiles(smi_heavy)
        assert mol is not None
        assert Descriptors.ExactMolWt(mol) > 250
        assert not _passes_ro3(mol)

    def test_mw_200_accepted(self):
        """A molecule with MW ≈ 200 and otherwise Ro3-compliant passes."""
        from rdkit import Chem
        from pipeline.stages.fragment.library import _passes_ro3
        from rdkit.Chem import Descriptors

        # Indole-3-acetic acid: C10H9NO2, MW ≈ 175 — clearly passes
        smi = "OC(=O)Cc1c[nH]c2ccccc12"
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None
        mw = Descriptors.ExactMolWt(mol)
        assert mw <= 250, f"Expected MW <= 250, got {mw:.1f}"
        assert _passes_ro3(mol)

    def test_hbd_over_3_rejected(self):
        """MW ok but HBD > 3 → rejected."""
        from rdkit import Chem
        from pipeline.stages.fragment.library import _passes_ro3
        from rdkit.Chem import rdMolDescriptors

        # Erythritol: C4H10O4, MW 122, HBD=4 → fails HBD
        smi = "OCC(O)C(O)CO"
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None
        assert rdMolDescriptors.CalcNumHBD(mol) > 3
        assert not _passes_ro3(mol)

    def test_logp_over_3_rejected(self):
        """MW ok but LogP > 3 → rejected."""
        from rdkit import Chem
        from pipeline.stages.fragment.library import _passes_ro3
        from rdkit.Chem import Descriptors

        # Phenylcyclohexane: C12H16, MW 160, LogP ~ 4.4
        smi = "C1CCC(CC1)c1ccccc1"
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None
        assert Descriptors.MolLogP(mol) > 3.0
        assert not _passes_ro3(mol)


# ── Test 2: Salt stripping ────────────────────────────────────────────────────

class TestSaltStripping:
    def test_chloride_salt_stripped(self):
        """'Cl.CCO' (ethanol hydrochloride) → 'CCO' (ethanol) after stripping."""
        from rdkit import Chem
        from rdkit.Chem import MolToSmiles
        from pipeline.stages.fragment.library import _strip_salts

        mol = Chem.MolFromSmiles("Cl.CCO")
        assert mol is not None
        stripped = _strip_salts(mol)
        assert stripped is not None
        result_smi = MolToSmiles(stripped)
        # The stripped molecule should be ethanol, not contain Cl
        assert "Cl" not in result_smi
        # Should be equivalent to CCO
        ref = Chem.MolFromSmiles("CCO")
        assert MolToSmiles(stripped) == MolToSmiles(ref)

    def test_single_component_unchanged(self):
        """A molecule without salts is returned unchanged."""
        from rdkit import Chem
        from rdkit.Chem import MolToSmiles
        from pipeline.stages.fragment.library import _strip_salts

        smi = "c1ccccc1"
        mol = Chem.MolFromSmiles(smi)
        stripped = _strip_salts(mol)
        assert stripped is not None
        assert MolToSmiles(stripped) == MolToSmiles(mol)

    def test_sodium_salt_stripped(self):
        """Sodium benzoate [Na+].[O-]C(=O)c1ccccc1 → benzoate."""
        from rdkit import Chem
        from rdkit.Chem import MolToSmiles
        from pipeline.stages.fragment.library import _strip_salts

        mol = Chem.MolFromSmiles("[Na+].[O-]C(=O)c1ccccc1")
        assert mol is not None
        stripped = _strip_salts(mol)
        # After stripping Na, should be the benzoate anion / benzoic acid form
        assert stripped is not None
        smi_out = MolToSmiles(stripped)
        assert "Na" not in smi_out


# ── Test 3: Murcko scaffold deduplication ────────────────────────────────────

class TestMurckoDedup:
    def test_same_scaffold_deduplicated(self):
        """Two benzene-substituted compounds share scaffold → only 1 retained."""
        from rdkit import Chem
        from pipeline.stages.fragment.library import _apply_filters
        from rich.console import Console

        # Both have a benzene ring as Murcko scaffold
        smiles_list = [
            "Cc1ccccc1",   # toluene (benzene scaffold)
            "CCc1ccccc1",  # ethylbenzene (benzene scaffold)
        ]
        quiet_console = Console(quiet=True)
        result = _apply_filters(smiles_list, quiet_console)
        # Both compounds pass Ro3, but share scaffold → only 1 should remain
        assert len(result) == 1

    def test_different_scaffolds_both_retained(self):
        """Pyridine and benzene have different scaffolds → both retained."""
        from pipeline.stages.fragment.library import _apply_filters
        from rich.console import Console

        smiles_list = [
            "Cc1cccnc1",   # methylpyridine (pyridine scaffold)
            "Cc1ccccc1",   # toluene (benzene scaffold)
        ]
        quiet_console = Console(quiet=True)
        result = _apply_filters(smiles_list, quiet_console)
        assert len(result) == 2

    def test_invalid_smiles_skipped(self):
        """Invalid SMILES is skipped gracefully."""
        from pipeline.stages.fragment.library import _apply_filters
        from rich.console import Console

        smiles_list = ["NOT_A_SMILES", "c1ccccc1"]
        quiet_console = Console(quiet=True)
        result = _apply_filters(smiles_list, quiet_console)
        # Only benzene passes
        assert len(result) == 1


# ── Test 4: Download failure → fallback ───────────────────────────────────────

class TestFallbackOnDownloadFailure:
    def test_fallback_used_when_download_raises(self, settings, cache, console):
        """When ZINC download fails, fallback library is loaded and warning emitted."""
        from pipeline.stages.fragment.library import run_library
        import httpx

        # Write a minimal fallback with known Ro3 compounds
        fallback_smiles = [
            "c1ccncc1",   # pyridine — passes Ro3
            "c1ccccc1",   # benzene — passes Ro3
            "CCO",        # ethanol — passes Ro3
        ]
        _write_fallback(settings, fallback_smiles)

        buf = StringIO()
        from rich.console import Console as RichConsole
        warn_console = RichConsole(file=buf, highlight=False)

        with patch(
            "pipeline.stages.fragment.library._fetch_tranche",
            side_effect=httpx.TransportError("connection failed"),
        ):
            lib_path = run_library(
                library_size=100,
                settings=settings,
                cache=cache,
                force=True,
                console=warn_console,
            )

        assert lib_path.exists()
        output = buf.getvalue()
        # Warning should mention fallback
        assert "fallback" in output.lower()

        # Result file should contain lines
        lines = [l for l in lib_path.read_text().splitlines() if l.strip()]
        assert len(lines) > 0
        # Fragment IDs should use FALLBACK prefix
        assert any("FALLBACK" in line for line in lines)

    def test_fallback_fragment_ids_use_fallback_prefix(self, settings, cache, console):
        """Fallback library entries use FALLBACK-##### IDs."""
        from pipeline.stages.fragment.library import run_library
        import httpx

        _write_fallback(settings, ["c1ccncc1", "CCO"])

        with patch(
            "pipeline.stages.fragment.library._fetch_tranche",
            side_effect=httpx.HTTPStatusError(
                "403", request=MagicMock(), response=MagicMock()
            ),
        ):
            lib_path = run_library(
                library_size=100,
                settings=settings,
                cache=cache,
                force=True,
                console=console,
            )

        content = lib_path.read_text()
        for line in content.strip().splitlines():
            parts = line.split("\t")
            assert len(parts) == 2
            assert parts[1].startswith("FALLBACK-")

    def test_cache_skip_returns_existing_path(self, settings, cache, console):
        """When cache is fresh and lib_path exists, SKIP is printed and path returned."""
        from pipeline.stages.fragment.library import run_library

        # Pre-create the library file and warm the cache
        lib_path = settings.cache_dir / "shared" / "fragment_library" / "fragments_ro3.smi"
        lib_path.parent.mkdir(parents=True, exist_ok=True)
        lib_path.write_text("c1ccccc1\tZINC-0000001\n")
        cache.save("shared", "fragment_library", {"path": str(lib_path), "size": 1}, 1)

        buf = StringIO()
        from rich.console import Console as RichConsole
        skip_console = RichConsole(file=buf, highlight=False)

        result = run_library(
            library_size=100,
            settings=settings,
            cache=cache,
            force=False,
            console=skip_console,
        )

        assert result == lib_path
        assert "skip" in buf.getvalue().lower()

    def test_library_sampled_to_requested_size(self, settings, cache, console):
        """Output contains at most library_size compounds."""
        from pipeline.stages.fragment.library import run_library

        # 5 unique Ro3 compounds in fallback (different scaffolds)
        fallback = [
            "c1ccncc1",   # pyridine
            "c1ccoc1",    # furan
            "c1ccsc1",    # thiophene
            "CCO",        # ethanol (no ring — empty Murcko → unique)
            "CC(=O)O",    # acetic acid (no ring — empty Murcko → may dedup with above)
        ]
        _write_fallback(settings, fallback)

        import httpx
        with patch(
            "pipeline.stages.fragment.library._fetch_tranche",
            side_effect=httpx.TransportError("fail"),
        ):
            lib_path = run_library(
                library_size=2,
                settings=settings,
                cache=cache,
                force=True,
                console=console,
            )

        lines = [l for l in lib_path.read_text().splitlines() if l.strip()]
        assert len(lines) <= 2
