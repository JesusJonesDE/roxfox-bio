# Implementation Plan: VRK1 Indication Validation and Chemistry Deepening

**Branch**: `004-vrk1-indication-chem` | **Date**: 2026-05-27 | **Spec**: [spec.md](spec.md)

## Summary

Four new pipeline stages that deepen the VRK1 drug discovery case built in spec-003: a co-mutation biomarker analysis for the CNS/Brain indication (US1), a VRK1 vs VRK2 paralog pocket comparison completing the `--include-vrk2` flag already wired in spec-003 (US2), AutoDock Vina computational docking of SCF-013 into 6AC9 (US3), and a co-crystal experimental brief compiled from RCSB crystallography data (US4). All four follow the existing `pipeline <command>` pattern with caching, `--force`, and rich console output.

## Technical Context

**Language/Version**: Python 3.11+ (matches existing pipeline)

**Primary Dependencies**:
- Existing: `httpx`, `pandas`, `typer`, `rich`, `tenacity`, `biopython`
- New (US1): `scipy>=1.11` (Fisher's exact test)
- New (US3, optional): `vina>=1.2.7`, `meeko>=0.5`, `rdkit-pypi` — marked as optional extras; pipeline fails gracefully when absent

**Storage**: File-based cache in `data/cache/{gene}/` and results in `data/results/{gene}/` — identical to spec-003

**Testing**: `pytest` — unit tests for enrichment logic, three-way classification, docking box geometry, contact mapping

**Target Platform**: macOS/Linux laptop (offline-capable for all stages except initial data downloads)

**Performance Goals**:
- Biomarker: < 10 min cold start (OmicsSomaticMutations.csv ~150 MB + Fisher's tests)
- Structalign + VRK2: < 5 min additional (VRK2 KLIFS lookup + AlphaFold download if needed)
- Docking: 1–3 min per scaffold on laptop (exhaustiveness=32, local Vina)
- Co-crystal brief: < 30 sec (single RCSB GraphQL query + local PDB analysis)
- All cache-hit runs: < 5 sec

**Constraints**: Docking is an optional dependency — `pipeline biomarker`, `structalign`, and `cocrystal` must work without `vina`/`meeko` installed

**Scale/Scope**: Single-target (VRK1), four scaffold comparisons maximum for docking in v1

## Constitution Check

No project-specific constitution is defined (template is blank). No gates to evaluate.

## Project Structure

### Documentation (this feature)

```text
specs/004-vrk1-indication-chem/
├── plan.md              ← this file
├── research.md          ✓ written
├── data-model.md        ✓ written
├── quickstart.md        ✓ written
├── contracts/
│   └── cli_commands.md  ✓ written
└── tasks.md             (next: /speckit.tasks)
```

### Source Code

New files follow the established `pipeline/stages/{stage}/` pattern:

```text
pipeline/stages/
├── biomarker/
│   ├── __init__.py
│   └── biomarker.py          # US1: enrichment analysis
├── structalign/
│   └── structalign.py        # US2: extend existing file (--include-vrk2 output)
├── dock/
│   ├── __init__.py
│   └── dock.py               # US3: Vina docking + contact mapping
└── cocrystal/
    ├── __init__.py
    └── cocrystal.py          # US4: crystallography brief

pipeline/
└── cli.py                    # Add `biomarker`, `dock`, `cocrystal` commands

pyproject.toml                # Add scipy; vina/meeko/rdkit as optional extras

tests/unit/
├── test_biomarker.py         # Fisher's exact enrichment logic
├── test_structalign_vrk2.py  # Three-way classification logic
└── test_dock.py              # Box geometry, contact mapping, score parsing
```

**Structure Decision**: Extend the existing `pipeline/stages/` layout. `structalign.py` is edited in place (no new directory). Three new stage directories mirror the spec-003 pattern exactly.

## Key Design Decisions

### US1 — Biomarker
- Download `OmicsSomaticMutations.csv` from the same DepMap manifest used in spec-003; exact-match filename then prefix fallback
- Filter: drop `Variant_Classification == "Silent"`, keep `isDeleterious == True`; deduplicate to one row per (ModelID, Hugo_Symbol) before contingency table construction
- Cache the filtered, per-gene mutation table as JSON (same `CacheManager` pattern)
- Significance threshold for report header: p < 0.05 — but write full table regardless

### US2 — VRK2 Three-Way
- Extend `_build_comparison()` and `_write_report()` in `structalign.py` to accept and use the VRK2 KLIFS DataFrame already fetched by `_fetch_vrk2_structure()`
- Add `_classify_three_way(vrk1_aa, vrk2_aa, egfr_aa) -> str` pure function returning one of four selectivity class labels
- `binding_site_comparison.csv` gains three new columns when `--include-vrk2` is passed; two-way behaviour is entirely unchanged otherwise

### US3 — Docking
- Check for `vina` import at stage entry; raise `ImportError` with install instructions if missing
- Box center = centroid of binding-site residue Cα coordinates (from structalign NeighborSearch output); box size = `2 × (max_radius + 3.0)` Å isotropic
- Receptor prep: call `mk_prepare_receptor.py` via `subprocess` (simpler and more reliable than meeko's Python receptor API which has changed between versions)
- Ligand prep: RDKit `EmbedMolecule` + `UFFOptimizeMolecule` → meeko `MoleculePreparation.prepare()` → PDBQT
- Positive control: extract ANP coordinates from 6AC9 HETATM records → write reference PDBQT → re-dock → compute heavy-atom RMSD vs crystallographic coordinates
- Contact mapping: for each pose, load PDBQT via RDKit, iterate heavy atoms, measure distance to each KLIFS selectivity-candidate residue Cα (from `binding_site_comparison.csv`), flag contacts ≤ 4Å

### US4 — Co-Crystal Brief
- RCSB GraphQL query: `exptl_crystal_grow { method pdbx_details }` — confirmed working for 6AC9 (returns PEG 3350 / ammonium sulfate / HEPES pH 7.0 conditions)
- Space group from PDB CRYST1 record (BioPython `Structure.header['cryst']` or raw text parse)
- Scaffold atom flagging: identify rotatable bonds (RDKit `rdMolDescriptors.CalcNumRotatableBonds`); flag atoms adjacent to rotatable bonds that protrude > 8Å from binding site centroid as potential packing clash risks
- Resolution requirements table: hardcoded per subpocket (from research.md); not dynamically computed
- Homolog fallback: if only one VRK1 entry found (6AC9), note it explicitly and append generic kinase soaking guidance as secondary recommendation
