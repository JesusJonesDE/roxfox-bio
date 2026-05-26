# Tasks: Drug Discovery Pipeline Research Tool

**Input**: Design documents from `specs/002-pipeline-research/`

**Branch**: `002-pipeline-research`

**Stack**: Python 3.11 | Typer + Rich | httpx + tenacity | RDKit (conda-forge) | pandas | Jinja2

**Organization**: Tasks grouped by user story — each phase is independently deliverable.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project scaffolding, packaging, environment configuration

- [x] T001 Create project directory structure: `pipeline/`, `pipeline/stages/fetch/`, `pipeline/stages/analyze/`, `pipeline/stages/report/`, `tests/unit/`, `tests/integration/fixtures/`, `data/` at repo root
- [x] T002 Create `pyproject.toml` — package name `pipeline`, Python 3.11+, dependencies (typer, rich, httpx, tenacity, pandas, jinja2), entry point `pipeline = "pipeline.cli:app"`
- [x] T003 [P] Create `environment.yml` — conda env `rxpipeline`, Python 3.11, rdkit from conda-forge, pip dependencies from pyproject.toml
- [x] T004 [P] Create `.gitignore` — ignore `data/`, `__pycache__/`, `*.pyc`, `.pytest_cache/`, `dist/`, `*.egg-info/`, `.env`

**Checkpoint**: Project installs with `pip install -e .` in the conda env; `pipeline --help` prints the top-level help.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Shared models, target config, cache engine, and CLI skeleton — everything user story stages depend on

**⚠️ CRITICAL**: All user story phases build on these components

- [x] T005 Create `pipeline/models.py` — dataclasses for `Target`, `BioactivityRecord`, `Compound`, `Scaffold`, `Structure`, `CacheEntry`; field definitions per `data-model.md`
- [x] T006 [P] Create `pipeline/config.py` — `TARGETS` dict mapping gene names to `Target` instances (VRK1/Q99986/RXF-001, IGHMBP2/P38935/RXF-002, VCP/P55072/RXF-003); `Settings` dataclass (data_dir, cache_max_age_days=30)
- [x] T007 Create `pipeline/cache.py` — `CacheManager` class: `get_cache_path(target, source)`, `is_fresh(target, source, max_age_days)`, `save(target, source, data, record_count)`, `load(target, source)`, `mark_stage_complete(stage, target)`, `is_stage_complete(stage, target)`, `reset_stage(stage, target)`, read/write `data/pipeline_manifest.json`
- [x] T008 Create `pipeline/cli.py` — Typer app with five subcommands registered: `fetch`, `analyze`, `report`, `run`, `status`; shared options `--target`, `--all`, `--force`, `--max-age`, `--data-dir`; Rich console for all output; stub implementations that print "not yet implemented"

**Checkpoint**: `pipeline status` runs and prints an empty table; `pipeline fetch --target VRK1` prints "not yet implemented".

---

## Phase 3: User Story 1 — Run Full Pipeline with Idempotency (Priority: P1) 🎯 MVP

**Goal**: `pipeline run --target VRK1` executes all stages end-to-end; re-runs skip cached stages; `--force` overrides selectively.

**Independent Test**: Run `pipeline run --target VRK1` from scratch — confirm `data/cache/VRK1/` contains JSON files for all 5 sources, `data/results/VRK1/` contains all 4 output files, and a completion summary is printed. Re-run immediately — confirm all stages report SKIP in < 30 seconds.

### Implementation

- [x] T009 [US1] Implement `pipeline/stages/fetch/chembl.py` — `fetch_target_id(uniprot_id)` (GET target endpoint, extract `target_chembl_id`), `fetch_bioactivity(chembl_target_id)` (paginated GET activity endpoint, `activity_type__in=IC50,Ki,Kd`, `assay_type__in=B,F`, limit=1000, collect all pages, return list of raw dicts); httpx client with tenacity retry on 429/503
- [x] T010 [P] [US1] Implement `pipeline/stages/fetch/open_targets.py` — `fetch_genetic_evidence(ensembl_id)` using GraphQL POST to Open Targets API; query fields: target symbol, associated diseases with scores, tractability summary; resolve Ensembl ID from UniProt via ChEMBL target lookup or hardcode mapping
- [x] T011 [P] [US1] Implement `pipeline/stages/fetch/pdb.py` — `fetch_structures(uniprot_id)` using RCSB Search API POST query filtered to UniProt accession; for each PDB ID retrieve metadata via `https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}`; extract resolution, method, ligand IDs, chain IDs, deposition date
- [x] T012 [P] [US1] Implement `pipeline/stages/fetch/alphafold.py` — `fetch_model(uniprot_id)` GET `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`; extract model accession, mean pLDDT, download URL; return None gracefully if no model exists
- [x] T013 [P] [US1] Implement `pipeline/stages/fetch/clinical_trials.py` — `fetch_trials(gene_name)` GET ClinicalTrials.gov v2 API `https://clinicaltrials.gov/api/v2/studies` with `query.term={gene_name}`; extract study title, phase, status, sponsor, intervention names
- [x] T014 [US1] Implement `fetch` subcommand body in `pipeline/cli.py` — for each requested target: call all 5 fetch functions, check cache freshness before each call (skip if fresh), save results via `CacheManager`, print per-source status line (SKIP/OK, record count, elapsed); mark fetch stage complete in manifest
- [x] T015 [US1] Implement `run` subcommand body in `pipeline/cli.py` — call `fetch`, `analyze`, `report` stages in sequence for each target; pass `--force` through to each stage; print per-target section header; exit non-zero if any stage fails
- [x] T016 [US1] Implement `status` subcommand body in `pipeline/cli.py` — read manifest, print Rich table with columns: Target, Fetch (date/STALE/—), Analyze (date/STALE/—), Report (date/STALE/—), Status (✓/⚠/✗)

**Checkpoint**: US1 complete. `pipeline run --target VRK1` works end-to-end. Re-run is fully cached. `pipeline status` shows correct state.

---

## Phase 4: User Story 2 — Compound Data and Scaffold Analysis (Priority: P2)

**Goal**: `pipeline analyze --target VRK1` produces `compounds_filtered.csv` and `scaffolds.csv` with correct potency filter, Lipinski properties, and Murcko scaffold clusters.

**Independent Test**: After a fetch run, run `pipeline analyze --target VRK1`. Open `data/results/VRK1/compounds_filtered.csv` — confirm rows exist with IC50 ≤ 10µM, MW/logP/HBD/HBA/ro5_violations/passes_ro5 columns populated. Open `data/results/VRK1/scaffolds.csv` — confirm multiple scaffolds with compound counts and potency values.

### Implementation

- [x] T017 [US2] Implement `pipeline/stages/analyze/bioactivity.py` — `filter_and_enrich(raw_records)`: (1) normalise units to nM, (2) filter to ≤ 10,000 nM, (3) deduplicate by compound ID keeping best (lowest) value per assay type, (4) calculate Lipinski properties via RDKit for each unique SMILES (MW, logP, HBD, HBA, rotatable_bonds, ro5_violations, passes_ro5), (5) handle RDKit sanitization failures gracefully (skip + log), (6) return list of `Compound` objects; write `data/results/{target}/compounds_filtered.csv`
- [x] T018 [US2] Implement `pipeline/stages/analyze/scaffolds.py` — `cluster_by_scaffold(compounds)`: extract Murcko scaffold SMILES per compound via RDKit `MurckoScaffoldSmilesFromSmiles`, group by scaffold, compute cluster size + median/best potency per cluster, assign `SCF-{N}` IDs, return list of `Scaffold` objects; write `data/results/{target}/scaffolds.csv`; update `scaffold_id` field in compounds CSV
- [x] T019 [US2] Implement `pipeline/stages/analyze/selectivity.py` — `profile_selectivity(compounds)`: for each compound in the filtered list, query ChEMBL activity endpoint for that `molecule_chembl_id` across all targets; count distinct human targets with activity ≤ 1,000 nM; set `off_target_flags` and `selectivity_flag` (flag if > 3 unrelated targets); update compounds_filtered.csv with these columns
- [x] T020 [US2] Implement `analyze` subcommand body in `pipeline/cli.py` — check fetch stage is complete (error if not), call bioactivity → scaffolds → selectivity in sequence, print per-step status lines, mark analyze stage complete in manifest

**Checkpoint**: US2 complete. Compounds CSV and scaffolds CSV are correct and non-empty for VRK1. Selectivity flags are set.

---

## Phase 5: User Story 3 — Structural Data Inventory (Priority: P2)

**Goal**: `pipeline analyze --target VRK1` also produces `structures.csv` listing PDB entries and AlphaFold model with quality metadata.

**Independent Test**: After fetch run, run `pipeline analyze --target IGHMBP2`. Open `data/results/IGHMBP2/structures.csv` — confirm rows with PDB IDs, resolution values, method, has_ligand flag. Confirm AlphaFold row present if no/few PDB entries. For a target with no PDB structures, confirm AlphaFold row still present and CSV is non-empty.

### Implementation

- [x] T021 [US3] Implement `pipeline/stages/analyze/structures.py` — `build_inventory(pdb_cache, alphafold_cache)`: parse PDB cache into `Structure` objects (structure_id, source=PDB, resolution_angstrom, method, has_ligand, ligand_ids, chain_ids, deposition_date); parse AlphaFold cache into one `Structure` object (source=AlphaFold, mean_plddt); sort PDB entries by resolution ascending; write `data/results/{target}/structures.csv`
- [x] T022 [US3] Wire structural inventory into `analyze` subcommand in `pipeline/cli.py` — call `build_inventory` after existing analyze steps, print status line (PDB count + AlphaFold presence), update manifest

**Checkpoint**: US3 complete. `structures.csv` exists with correct data for all three targets.

---

## Phase 6: User Story 4 — Research Dossier (Priority: P3)

**Goal**: `pipeline report --target VRK1` produces `dossier.md` with all 8 required sections, fully populated from analysis results, readable without any other files.

**Independent Test**: After a full run for VCP, open `data/results/VCP/dossier.md`. Run `grep "^##" data/results/VCP/dossier.md` — confirm all 8 section headers present. Read each section — confirm data is present or a "no data" explanation is given; confirm no raw JSON, no code, and no placeholder text remains.

### Implementation

- [x] T023 [US4] Create `pipeline/stages/report/dossier_template.md` — Jinja2 template with all 8 sections: `# {{target.gene_name}} — Research Dossier`, `## Overview`, `## Genetic Evidence`, `## Bioactivity Summary`, `## Scaffold Highlights`, `## Structural Data`, `## Selectivity Profile`, `## Competitive Landscape`, `## Data Gaps & Limitations`; each section uses template variables to render data or a "no data available" message
- [x] T024 [US4] Implement `pipeline/stages/report/dossier.py` — `generate_dossier(target, results_dir, cache_dir)`: load all CSVs from results, load Open Targets + ClinicalTrials cache, build template context dict (compound counts, top scaffolds list, best structure, genetic scores, trial list, competitive drugs, data gap flags), render template via Jinja2, write `data/results/{target}/dossier.md`
- [x] T025 [US4] Implement `report` subcommand body in `pipeline/cli.py` — check analyze stage complete (error if not), call `generate_dossier`, print output path, mark report stage complete in manifest

**Checkpoint**: US4 complete. All three targets have complete, readable dossiers after `pipeline run --all`.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Tests, documentation, validation against quickstart.md scenarios

- [x] T026 Create `tests/unit/test_bioactivity.py` — unit tests: potency filter at boundary (10000 nM passes, 10001 nM fails), unit normalisation (µM → nM), Lipinski pass/fail at rule boundaries, dedup keeps best value, invalid SMILES handled gracefully
- [x] T027 [P] Create `tests/unit/test_scaffolds.py` — unit tests: Murcko extraction from known SMILES, clustering groups identical scaffolds, scaffold ID assignment is stable, sanitization failure skips compound without crash
- [x] T028 [P] Create `tests/unit/test_cache.py` — unit tests: freshness check with mock timestamps, manifest write/read round-trip, stage completion tracking, --force resets correct stages, corrupted manifest handled
- [x] T029 Save sample API responses as fixtures in `tests/integration/fixtures/` — `chembl_vrk1.json` (first page of bioactivity), `open_targets_vrk1.json`, `pdb_vrk1.json` (first 5 structures), `alphafold_vrk1.json`, `clinical_trials_vrk1.json`
- [x] T030 Create `tests/integration/test_fetch_clients.py` — test each fetch client against fixture files using httpx mock transport; verify parsed output matches expected `BioactivityRecord`/`Structure` counts
- [x] T031 Create `README.md` at project root — setup instructions (conda env + pip install), quickstart (first run command), data directory structure, how to add a new target
- [x] T032 Run full quickstart.md validation — execute all 7 scenarios in `quickstart.md`; confirm each expected result matches; fix any failures

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1
- **Phase 3 (US1)**: Depends on Phase 2 — **MVP stops here** (fetch works end-to-end)
- **Phase 4 (US2)**: Depends on Phase 2 — can work in parallel with Phase 3
- **Phase 5 (US3)**: Depends on Phase 2 — can work in parallel with Phases 3+4
- **Phase 6 (US4)**: Depends on Phases 4+5 (needs results from both analyze stages)
- **Phase 7 (Polish)**: Depends on all prior phases complete

### User Story Dependencies

- **US1**: Independent after Foundation — delivers end-to-end pipeline skeleton + caching
- **US2**: Independent after Foundation — delivers compound + scaffold analysis
- **US3**: Independent after Foundation — delivers structural inventory
- **US4**: Depends on US2+US3 being complete (dossier reads all result files)

### Parallel Opportunities

```
# Phase 1 — run together:
T003 (environment.yml) + T004 (.gitignore)

# Phase 2 — T006 can run with T005+T007+T008:
T006 (config.py) can be written while T005 (models.py) is in progress

# Phase 3 — fetch clients run together after T009 starts:
T010 (open_targets) + T011 (pdb) + T012 (alphafold) + T013 (clinical_trials)

# Phase 7 — tests run together:
T026 (test_bioactivity) + T027 (test_scaffolds) + T028 (test_cache)
```

---

## Implementation Strategy

### MVP (Phases 1–3 only)

1. Complete Phase 1: scaffolding
2. Complete Phase 2: models, cache, CLI skeleton
3. Complete Phase 3: all 5 fetch clients + run/status commands
4. **STOP and VALIDATE**: `pipeline run --target VRK1` succeeds; re-run is cached; `status` shows correct state
5. All raw data is in `data/cache/VRK1/` — manually inspect JSON files

**MVP delivers**: an idempotent data fetcher with caching — usable for manual inspection immediately.

### Full Research Tool (All Phases)

Add Phase 4 (compound analysis) → Phase 5 (structural analysis) → Phase 6 (dossier) → Phase 7 (tests + docs).

**Total tasks**: 32 | **MVP tasks**: 16 | **Parallel opportunities**: 3 batches

---

## Notes

- No TDD requested — tests are in Phase 7 (Polish), written after implementation
- RDKit must be installed via conda-forge; document this clearly in error messages if import fails
- `data/` directory is fully gitignored — never commit API responses or results
- Each fetch client is an isolated module; failures in one do not block others
- The `pipeline run` command is the primary user-facing interface; subcommands are for granular control
