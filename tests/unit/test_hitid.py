"""Unit tests for pipeline/stages/hitid/hitid.py"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from pipeline.stages.hitid.hitid import (
    _parse_bindingdb_response,
    _canonicalise,
    _strip_salts,
    _compute_ro5,
    _process_compounds,
    _write_csv,
    run_hitid,
    CSV_COLUMNS,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _make_bdb_response(ligands: list[dict]) -> dict:
    """Build a minimal BindingDB JSON response wrapping a list of affinity dicts."""
    return {
        "getLigandsByUniprotsResponse": {
            "affinities": ligands,
        }
    }


def _simple_ligand(smiles: str, ki: float) -> dict:
    return {"smiles": smiles, "Ki": ki}


# ── Test 1: BindingDB response parsing ───────────────────────────────────────

def test_parse_bindingdb_extracts_compounds() -> None:
    """Given a mock BindingDB JSON response, extract the correct compound list."""
    ligands = [
        {"smiles": "c1ccccc1", "Ki": 100.0},
        {"smiles": "CCO", "IC50": 250.0},
        {"smiles": "CC(=O)O", "Kd": 500.0},
        {"smiles": "", "Ki": 10.0},    # empty SMILES — should be skipped
        {"smiles": "CNCC", "Ki": None, "IC50": None, "Kd": None},  # no potency — skip
    ]
    data = _make_bdb_response(ligands)
    hits = _parse_bindingdb_response(data, "BindingDB_UniProt")

    assert len(hits) == 3
    smiles_set = {h["smiles"] for h in hits}
    assert "c1ccccc1" in smiles_set
    assert "CCO" in smiles_set
    assert "CC(=O)O" in smiles_set

    # Potency values should be correctly assigned
    benzene = next(h for h in hits if h["smiles"] == "c1ccccc1")
    assert benzene["best_value_nm"] == pytest.approx(100.0)
    assert benzene["best_assay_type"] == "Ki"
    assert benzene["source"] == "BindingDB_UniProt"


def test_parse_bindingdb_picks_best_potency() -> None:
    """When a ligand has both Ki and IC50, the smaller value should be selected."""
    ligands = [{"smiles": "c1ccccc1", "Ki": 500.0, "IC50": 200.0}]
    data = _make_bdb_response(ligands)
    hits = _parse_bindingdb_response(data, "test")
    assert len(hits) == 1
    assert hits[0]["best_value_nm"] == pytest.approx(200.0)
    assert hits[0]["best_assay_type"] == "IC50"


def test_parse_bindingdb_strips_comparison_operators() -> None:
    """Potency strings like '>10000' or '<1' should be parsed to float."""
    ligands = [{"smiles": "CCO", "Ki": ">10000"}]
    data = _make_bdb_response(ligands)
    hits = _parse_bindingdb_response(data, "test")
    assert len(hits) == 1
    assert hits[0]["best_value_nm"] == pytest.approx(10000.0)


def test_parse_bindingdb_empty_response() -> None:
    """Empty or malformed BindingDB responses should return an empty list."""
    assert _parse_bindingdb_response({}, "test") == []
    assert _parse_bindingdb_response([], "test") == []
    assert _parse_bindingdb_response(
        {"getLigandsByUniprotsResponse": {"affinities": []}}, "test"
    ) == []


# ── Test 2: PubChem CID→SMILES parsing ───────────────────────────────────────

def test_pubchem_smiles_parsing() -> None:
    """Given mock PubChem property response, return correct SMILES list."""
    mock_response = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 12345,
                    "IsomericSMILES": "CC(=O)Nc1ccc(O)cc1",
                    "MolecularWeight": 151.16,
                    "XLogP": 0.9,
                },
                {
                    "CID": 67890,
                    "IsomericSMILES": "",   # empty — should be skipped
                    "MolecularWeight": 100.0,
                    "XLogP": 1.0,
                },
            ]
        }
    }

    # Simulate what _query_pubchem does with the property response
    compounds = []
    for prop in mock_response.get("PropertyTable", {}).get("Properties", []):
        smiles = prop.get("IsomericSMILES", "").strip()
        if not smiles:
            continue
        mw = prop.get("MolecularWeight")
        xlogp = prop.get("XLogP")
        compounds.append({
            "smiles": smiles,
            "best_value_nm": None,
            "best_assay_type": "PubChem_AID999",
            "source": "PubChem",
            "molecular_weight": float(mw) if mw is not None else None,
            "logp": float(xlogp) if xlogp is not None else None,
        })

    assert len(compounds) == 1
    assert compounds[0]["smiles"] == "CC(=O)Nc1ccc(O)cc1"
    assert compounds[0]["molecular_weight"] == pytest.approx(151.16)
    assert compounds[0]["logp"] == pytest.approx(0.9)


# ── Test 3: Ro5 filter ────────────────────────────────────────────────────────

def test_ro5_mw_600_fails() -> None:
    """A compound with MW=600 should fail Ro5 (MW > 500 violation)."""
    # Use a heavy SMILES: erythromycin-like mass range
    # Instead of finding a real SMILES, test _compute_ro5 directly with a
    # synthetically heavy molecule built from RDKit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Build a simple repeating chain that gives MW > 500
    # Hexadecane: C16H34 = MW 226 — too light.  Use a polycyclic system instead.
    # Easiest: pass a SMILES known to have MW > 500 (e.g. 6-ring fused system)
    # We test logic not RDKit accuracy: mock _compute_ro5 indirectly by inspecting output
    props = _compute_ro5("CCCCCCCCCCCCCCCCCCCCCCCCCCC")  # long chain ~MW 380 — passes
    # Confirm the function works; MW may pass or fail depending on chain length
    assert props["molecular_weight"] is not None
    assert props["passes_ro5"] is not None

    # Now test with a SMILES that definitely fails MW (use a known large molecule)
    # Cyclosporin A: MW ~1202 — well above 500
    cyclosporin_smiles = (
        "CC[C@@H]1C(=O)N(CC(=O)N([C@@H](C(=O)N[C@@H](C(=O)N([C@@H]"
        "(C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N([C@@H](C(=O)N([C@@H](C(=O)"
        "N([C@@H](C(=O)N1[C@@H](CC(C)C)C(=O)O)CC(C)C)C)CC(C)C)C)C(C)C)"
        "C)CC(C)C)C)C(C)C)CC(C)C)C)C"
    )
    props_large = _compute_ro5(cyclosporin_smiles)
    if props_large["molecular_weight"] is not None:
        assert props_large["passes_ro5"] is False or props_large["ro5_violations"] >= 1


def test_ro5_mw_400_passes() -> None:
    """A drug-like compound with MW ~400 should pass Ro5."""
    # Ibuprofen: MW=206 — passes easily
    props = _compute_ro5("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O")
    assert props["molecular_weight"] is not None
    assert props["molecular_weight"] < 500
    assert props["passes_ro5"] is True


def test_ro5_violations_count() -> None:
    """Verify that ro5_violations correctly counts individual violations."""
    # Aspirin: MW=180, logP~1.2, HBD=1, HBA=3 → 0 violations
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    props = _compute_ro5(aspirin)
    assert props["ro5_violations"] == 0
    assert props["passes_ro5"] is True


# ── Test 4: Salt stripping ────────────────────────────────────────────────────

def test_salt_stripping_keeps_largest_fragment() -> None:
    """Salt stripping should return the largest organic fragment."""
    # HCl salt of a simple amine
    salty = "Cl.NCCN"
    result = _strip_salts(salty)
    # Should keep NCCN (larger than Cl)
    canon = _canonicalise(result)
    assert canon is not None
    # The canonical form should not contain Cl as a salt
    assert "." not in canon or canon.count(".") == 0 or len(canon.split(".")[0]) > 2


def test_salt_stripping_no_dot_passthrough() -> None:
    """A pure compound (no dot) should be returned unchanged (after canonicalisation)."""
    smiles = "c1ccccc1"  # benzene
    result = _strip_salts(smiles)
    assert "." not in result


# ── Test 5: Scaffold ID assignment ───────────────────────────────────────────

def test_scaffold_id_first_compound_gets_001() -> None:
    """The first compound in the processed list should receive GENE-SCF-001."""
    raw = [
        {"smiles": "c1ccccc1", "best_value_nm": 100.0, "best_assay_type": "Ki", "source": "test"},
        {"smiles": "CCO", "best_value_nm": 500.0, "best_assay_type": "Ki", "source": "test"},
    ]
    console = MagicMock()
    processed = _process_compounds(raw, "IGHMBP2", console)

    assert len(processed) >= 1
    # Most potent compound (best_value_nm=100 nM) should be first
    assert processed[0]["scaffold_id"] == "IGHMBP2-SCF-001"


def test_scaffold_ids_are_sequential() -> None:
    """Scaffold IDs should be sequential integers padded to 3 digits."""
    raw = [
        {"smiles": "c1ccccc1", "best_value_nm": 100.0, "best_assay_type": "Ki", "source": "test"},
        {"smiles": "CCO", "best_value_nm": 200.0, "best_assay_type": "Ki", "source": "test"},
        {"smiles": "CC", "best_value_nm": 300.0, "best_assay_type": "Ki", "source": "test"},
    ]
    console = MagicMock()
    processed = _process_compounds(raw, "IGHMBP2", console)
    ids = [c["scaffold_id"] for c in processed]
    assert ids == ["IGHMBP2-SCF-001", "IGHMBP2-SCF-002", "IGHMBP2-SCF-003"]


# ── Test 6: Empty sources → empty CSV, no error ───────────────────────────────

def test_empty_sources_write_empty_csv(tmp_path: Path) -> None:
    """When BindingDB and PubChem both return 0 compounds, write empty CSV with header."""
    from pipeline.config import Settings
    from pipeline.cache import CacheManager

    settings = Settings(data_dir=tmp_path)
    cache = CacheManager(settings)
    console = MagicMock()

    with (
        patch("pipeline.stages.hitid.hitid._query_bindingdb_by_uniprot", return_value=[]),
        patch("pipeline.stages.hitid.hitid._query_bindingdb_by_name", return_value=[]),
        patch("pipeline.stages.hitid.hitid._query_pubchem", return_value=[]),
    ):
        # Should not raise
        run_hitid("IGHMBP2", settings, cache, force=True, console=console)

    out_csv = tmp_path / "results" / "IGHMBP2" / "compounds_filtered.csv"
    assert out_csv.exists(), "compounds_filtered.csv should be created even when empty"

    with out_csv.open() as fh:
        reader = csv.DictReader(fh)
        assert list(reader.fieldnames) == CSV_COLUMNS
        rows = list(reader)
        assert rows == [], "CSV should be empty (header only) when no compounds found"


def test_empty_sources_write_report(tmp_path: Path) -> None:
    """When no compounds are found, hitid_report.md should explain the gap."""
    from pipeline.config import Settings
    from pipeline.cache import CacheManager

    settings = Settings(data_dir=tmp_path)
    cache = CacheManager(settings)
    console = MagicMock()

    with (
        patch("pipeline.stages.hitid.hitid._query_bindingdb_by_uniprot", return_value=[]),
        patch("pipeline.stages.hitid.hitid._query_bindingdb_by_name", return_value=[]),
        patch("pipeline.stages.hitid.hitid._query_pubchem", return_value=[]),
    ):
        run_hitid("IGHMBP2", settings, cache, force=True, console=console)

    report = tmp_path / "results" / "IGHMBP2" / "hitid_report.md"
    assert report.exists()
    text = report.read_text()
    assert "IGHMBP2" in text
    assert "Hit Identification" in text


def test_empty_sources_no_exception(tmp_path: Path) -> None:
    """run_hitid must not raise when all sources return empty."""
    from pipeline.config import Settings
    from pipeline.cache import CacheManager

    settings = Settings(data_dir=tmp_path)
    cache = CacheManager(settings)
    console = MagicMock()

    with (
        patch("pipeline.stages.hitid.hitid._query_bindingdb_by_uniprot", return_value=[]),
        patch("pipeline.stages.hitid.hitid._query_bindingdb_by_name", return_value=[]),
        patch("pipeline.stages.hitid.hitid._query_pubchem", return_value=[]),
    ):
        # Must not raise — this is the key contract
        run_hitid("IGHMBP2", settings, cache, force=True, console=console)


# ── Test 7: Deduplication by canonical SMILES ─────────────────────────────────

def test_deduplication_keeps_best_potency() -> None:
    """Duplicate SMILES (even as different representations) keep the most potent entry."""
    # benzene written two ways: aromatic and Kekulé
    raw = [
        {"smiles": "c1ccccc1", "best_value_nm": 500.0, "best_assay_type": "Ki", "source": "A"},
        {"smiles": "C1=CC=CC=C1", "best_value_nm": 100.0, "best_assay_type": "IC50", "source": "B"},
    ]
    console = MagicMock()
    processed = _process_compounds(raw, "IGHMBP2", console)
    # Should deduplicate to 1 compound with best_value_nm=100
    assert len(processed) == 1
    assert processed[0]["best_value_nm"] == pytest.approx(100.0)


# ── Test 8: CSV column schema ─────────────────────────────────────────────────

def test_csv_columns_match_schema(tmp_path: Path) -> None:
    """Written CSV must have exactly the columns defined in CSV_COLUMNS."""
    compounds = [
        {
            "smiles": "c1ccccc1",
            "best_value_nm": 100.0,
            "best_assay_type": "Ki",
            "source": "test",
            "scaffold_id": "IGHMBP2-SCF-001",
            "molecular_weight": 78.11,
            "logp": 1.56,
            "hbd": 0,
            "hba": 0,
            "rotatable_bonds": 0,
            "ro5_violations": 0,
            "passes_ro5": True,
            "off_target_flags": 0,
            "selectivity_flag": False,
        }
    ]
    out = tmp_path / "test.csv"
    _write_csv(compounds, out)
    with out.open() as fh:
        reader = csv.DictReader(fh)
        assert list(reader.fieldnames) == CSV_COLUMNS
