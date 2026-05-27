# Tasks: VRK1 Indication Validation and Chemistry Deepening

**Input**: Design documents from `/specs/004-vrk1-indication-chem/`

**Branch**: `004-vrk1-indication-chem`

**Format**: `[ID] [P?] [Story?] Description with file path`
- **[P]**: Parallelisable (different files, no blocking dependencies)
- **[Story]**: US1 = Biomarker, US2 = VRK2 Three-Way, US3 = Docking, US4 = Co-Crystal

---

## Phase 1: Setup

**Purpose**: Add new dependencies and prepare package skeleton.

- [ ] T001 Add `scipy>=1.11` to `[project.dependencies]` and `[project.optional-dependencies] docking = ["vina>=1.2.7", "meeko>=0.5", "rdkit-pypi"]` in pyproject.toml

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Package init files required before any stage code can be imported.

**⚠️ CRITICAL**: All user story work depends on these being present.

- [ ] T002 [P] Create empty `pipeline/stages/biomarker/__init__.py`
- [ ] T003 [P] Create empty `pipeline/stages/dock/__init__.py`
- [ ] T004 [P] Create empty `pipeline/stages/cocrystal/__init__.py`

**Checkpoint**: Package structure ready — all four user story phases can begin.

---

## Phase 3: User Story 1 — CNS/Brain Biomarker Analysis (Priority: P1) 🎯 MVP

**Goal**: `pipeline biomarker --target VRK1 --lineage "CNS/Brain"` downloads `OmicsSomaticMutations.csv`, runs Fisher's exact test per gene against the cached DepMap dependency data, and writes a ranked co-mutation report.

**Independent Test**: Run command cold; verify `biomarker_report.md` and `biomarker_results.csv` appear in `data/results/VRK1/`; re-run and verify `SKIP (cached)` in < 5 s.

- [ ] T005 [P] [US1] Implement `_download_mutations(manifest: dict, cache_dir: Path) -> pd.DataFrame` — exact-match `OmicsSomaticMutations.csv` from DepMap manifest, filter `Variant_Classification != "Silent"` and `isDeleterious == True`, deduplicate to one row per (ModelID, Hugo_Symbol), cache result in `pipeline/stages/biomarker/biomarker.py`
- [ ] T006 [P] [US1] Implement `_compute_enrichment(dep_records: list, mutations: pd.DataFrame, lineage: str, min_lines: int) -> pd.DataFrame` — split cell lines into strongly-dependent (gene_effect ≤ −0.5) vs non-dependent, build 2×2 contingency per gene, call `scipy.stats.fisher_exact(alternative="greater")`, return BiomarkerResult DataFrame sorted by p_value ascending, in `pipeline/stages/biomarker/biomarker.py`
- [ ] T007 [US1] Implement `_write_report()`, `update_research_report()`, and `run_biomarker(gene_symbol, lineage, settings, cache, force, min_lines, console)` — orchestrate download → enrich → write `biomarker_results.csv` + `biomarker_report.md` (Markdown table: gene, odds_ratio, p_value, n_dep_with_mut, n_nondep_with_mut, direction) + inject into `research_report.md`; follow spec-003 SKIP/force/console pattern, in `pipeline/stages/biomarker/biomarker.py`
- [ ] T008 [US1] Add `biomarker` Typer command to `pipeline/cli.py` with options: `--target/-t TEXT`, `--lineage TEXT`, `--all`, `--force`, `--min-lines INT` (default 3), `--data-dir PATH`; import and call `run_biomarker()`
- [ ] T009 [P] [US1] Write `tests/unit/test_biomarker.py` — test `_compute_enrichment()` with synthetic dep_records and mutation DataFrame: (a) gene enriched in dependent lines → odds_ratio > 1, direction = "enriched_in_dependent"; (b) exactly at min-lines threshold → included; (c) below threshold → excluded; (d) no enriched genes → report states "no significant biomarker found"

**Checkpoint**: US1 fully functional and cached. `biomarker_report.md` readable by a non-technical stakeholder.

---

## Phase 4: User Story 2 — VRK2 Three-Way Comparison (Priority: P2)

**Goal**: `pipeline structalign --target VRK1 --include-vrk2 --force` extends the existing two-way comparison with a `selectivity_class` column and a VRK2 section in the report. Running without the flag is identical to spec-003.

**Independent Test**: Run with `--include-vrk2`; verify `binding_site_comparison.csv` has `selectivity_class` column with 100% coverage across 85 positions; run without flag and diff CSV — must be identical to spec-003 output.

- [X] T010 [P] [US2] Add `_classify_three_way(vrk1_aa: str, vrk2_aa: str, egfr_aa: str) -> str` pure function to `pipeline/stages/structalign/structalign.py` — returns "VRK1-specific" (vrk1≠vrk2 AND vrk1≠egfr), "pan-VRK vs EGFR" (vrk1==vrk2 AND vrk1≠egfr), "VRK2 vs VRK1+EGFR" (vrk1≠vrk2 AND vrk1==egfr), or "conserved" (all equal); treat GAP/LOW_CONF as non-matching
- [X] T011 [US2] Extend `_build_comparison()` in `pipeline/stages/structalign/structalign.py` to populate three new columns when `vrk2_klifs` is not None: `vrk2_aa`, `vrk1_vrk2_diff` (using existing `_classify_difference()`), `selectivity_class` (using `_classify_three_way()`); when vrk2_klifs is None columns are omitted so two-way output is unchanged
- [X] T012 [US2] Extend `_write_report()` in `pipeline/stages/structalign/structalign.py` to add a "## VRK2 Three-Way Comparison" section after the existing selectivity hypothesis when `vrk2_source` is provided — show: VRK2 source label, count of VRK1-specific / pan-VRK / conserved positions, gatekeeper three-way comparison row, and a table of VRK1-specific selectivity candidates (positions where vrk1≠vrk2 AND selectivity_candidate==True)
- [X] T013 [P] [US2] Write `tests/unit/test_structalign_vrk2.py` — test `_classify_three_way()` for all four class labels including GAP and LOW_CONF inputs; test that running `_build_comparison()` with `vrk2_klifs=None` produces identical column set to spec-003 (no vrk2 columns present)

**Checkpoint**: US2 complete. Three-way report visible; two-way regression confirmed.

---

## Phase 5: User Story 3 — SCF-013 Docking (Priority: P3)

**Goal**: `pipeline dock --target VRK1 --scaffold SCF-013` runs AutoDock Vina locally, validates with ANP positive control, and maps SCF-013 poses onto the 10 selectivity candidates from spec-003.

**Independent Test**: (a) With Vina absent → clean error + install instructions. (b) With Vina present → `docking_report.md` exists, control RMSD reported, top-3 poses listed, selectivity contacts mapped. Cache hit < 5 s.

- [X] T014 [P] [US3] Implement `_check_vina_installed()` (import test → raises `RuntimeError` with `pip install vina meeko rdkit-pypi` message if missing), `_prepare_receptor(pdb_path, cache_dir) -> Path` (call `mk_prepare_receptor.py` via `subprocess.run`, cache PDBQT), and `_prepare_ligand(smiles, scaffold_id, cache_dir) -> Path` (RDKit `EmbedMolecule` + `UFFOptimizeMolecule` → meeko `MoleculePreparation.prepare()` → PDBQT, cache result) in `pipeline/stages/dock/dock.py`
- [X] T015 [P] [US3] Implement `_define_box(binding_site_residues: list) -> tuple[list, list]` — compute centroid of Cα coordinates, max radius, return `(center_xyz, [2*(r+3)]*3)` as Vina box args; and `_run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness, output_pdbqt) -> list[dict]` — run Vina, parse `REMARK VINA RESULT:` lines into list of `{pose_rank, affinity_kcal_mol, rmsd_lb, rmsd_ub}` dicts, in `pipeline/stages/dock/dock.py`
- [X] T016 [US3] Implement `_run_control(pdb_path, receptor_pdbqt, center, box_size, cache_dir) -> float` — extract ANP HETATM coordinates from 6AC9 PDB, write reference PDBQT, re-dock with exhaustiveness=16, compute heavy-atom RMSD of top pose vs crystallographic coordinates using BioPython; and `_map_contacts(pose_pdbqt_path, comparison_csv_path, cutoff_A=4.0) -> list[int]` — load pose heavy atoms, measure distance to each selectivity-candidate residue Cα, return list of contacted KLIFS positions, in `pipeline/stages/dock/dock.py`
- [X] T017 [US3] Implement `_write_report()` and `run_dock(gene_symbol, scaffold_id, settings, cache, force, exhaustiveness, console)` — orchestrate: check Vina → resolve scaffold SMILES from `compounds_filtered.csv` → prepare receptor/ligand → define box from `binding_site_vrk1.csv` Cα coords → run control → run SCF-013 → map contacts → write `docking_results_{scaffold}.csv` + `docking_poses_{scaffold}.pdbqt` + `docking_report.md`; follow spec-003 SKIP/force/console pattern, in `pipeline/stages/dock/dock.py`
- [X] T018 [US3] Add `dock` Typer command to `pipeline/cli.py` with options: `--target/-t TEXT`, `--scaffold TEXT`, `--all-scaffolds`, `--force`, `--exhaustiveness INT` (default 32), `--data-dir PATH`; catch `RuntimeError` from `_check_vina_installed()` and exit with code 1
- [X] T019 [P] [US3] Write `tests/unit/test_dock.py` — test `_define_box()` with known coordinates (centroid and box size correct); test `_run_vina()` score parser with hardcoded PDBQT REMARK lines; test missing-Vina error message content; test `_map_contacts()` with synthetic atom positions and known KLIFS positions

**Checkpoint**: US3 complete. Docking report maps SCF-013 contacts to VRK1 selectivity handles.

---

## Phase 6: User Story 4 — Co-Crystal Structure Brief (Priority: P4)

**Goal**: `pipeline cocrystal --target VRK1 --scaffold SCF-013` queries RCSB for 6AC9 crystallisation conditions and generates a structured experimental brief.

**Independent Test**: Run command; verify `cocrystal_brief.md` cites PEG 3350/ammonium sulfate/HEPES pH 7.0 conditions and states minimum resolution requirements per subpocket.

- [X] T020 [P] [US4] Implement `_fetch_rcsb_conditions(pdb_id: str) -> dict` — POST to `https://data.rcsb.org/graphql` with query `{ entries(entry_ids: ["PDB_ID"]) { exptl_crystal_grow { method pdbx_details } } }`, return `{method, conditions}` dict; fall back to generic kinase soaking guidance if request fails or returns empty; and `_get_space_group(pdb_path: Path) -> str` — parse CRYST1 record line from PDB file, in `pipeline/stages/cocrystal/cocrystal.py`
- [X] T021 [P] [US4] Implement `_flag_scaffold_atoms(smiles: str, binding_site_centroid: list[float]) -> list[str]` — generate RDKit 3D conformer, align to centroid, identify heavy atoms > 8Å from centroid that are adjacent to rotatable bonds (`rdMolDescriptors.CalcNumRotatableBonds`), return list of atom descriptions; return empty list if SMILES cannot be parsed (do not raise), in `pipeline/stages/cocrystal/cocrystal.py`
- [X] T022 [US4] Implement `_write_brief()` and `run_cocrystal(gene_symbol, scaffold_id, settings, cache, force, console)` — check `binding_site_comparison.csv` exists (exit code 1 with actionable message if not), resolve scaffold SMILES, fetch RCSB conditions, parse space group, flag scaffold atoms, compute centroid from `binding_site_vrk1.csv`, write `cocrystal_brief.md` with sections: Structures, Crystallisation Conditions, Resolution Requirements (hardcoded per subpocket from research.md), Scaffold Compatibility, Recommended Experiment; follow spec-003 SKIP/force pattern, in `pipeline/stages/cocrystal/cocrystal.py`
- [X] T023 [US4] Add `cocrystal` Typer command to `pipeline/cli.py` with options: `--target/-t TEXT`, `--scaffold TEXT`, `--force`, `--data-dir PATH`

**Checkpoint**: US4 complete. `cocrystal_brief.md` is a self-contained experimental planning document.

---

## Phase 7: Polish & Validation

- [X] T024 Run quickstart.md end-to-end validation — cold start for all four commands in sequence (`biomarker` → `structalign --include-vrk2` → `dock` → `cocrystal`), verify all output files match quickstart.md checklists, confirm cache-hit runs complete in < 5 s
- [X] T025 [P] Verify structalign regression — run `pipeline structalign --target VRK1 --force` (no `--include-vrk2`), diff `binding_site_comparison.csv` column set against spec-003 output (must have no `vrk2_aa` or `selectivity_class` columns)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (T001)**: No dependencies — start immediately
- **Foundational (T002–T004)**: Depends on T001 — blocks all story phases
- **US1 (T005–T009)**: Depends on Foundational — no dependencies on US2/3/4
- **US2 (T010–T013)**: Depends on Foundational — no dependencies on US1/3/4; edits existing structalign.py
- **US3 (T014–T019)**: Depends on Foundational — reads spec-003 outputs (binding_site_vrk1.csv, compounds_filtered.csv); no dependency on US1/2/4
- **US4 (T020–T023)**: Depends on Foundational — reads spec-003 structalign outputs; no dependency on US1/2/3
- **Polish (T024–T025)**: Depends on all US phases complete

### Within Each User Story

- US1: T005 and T006 are parallel → T007 depends on both → T008 depends on T007
- US2: T010 is parallel → T011 depends on T010 → T012 depends on T011
- US3: T014 and T015 are parallel → T016 depends on both → T017 depends on T016 → T018 depends on T017
- US4: T020 and T021 are parallel → T022 depends on both → T023 depends on T022

### Parallel Opportunities

All four user story phases can execute simultaneously after Foundational completes. Within each story, tasks marked [P] run in parallel.

---

## Parallel Example: US3 (Docking)

```
# Launch in parallel (different function groups, same file):
Task T014: _check_vina_installed, _prepare_receptor, _prepare_ligand
Task T015: _define_box, _run_vina

# Sequential after T014 + T015 complete:
Task T016: _run_control, _map_contacts
Task T017: run_dock + _write_report
Task T018: CLI command
```

---

## Implementation Strategy

### MVP (US1 only — 5 tasks)

1. T001 Setup → T002 Foundational → T005+T006 parallel → T007 → T008
2. **Validate**: `pipeline biomarker --target VRK1 --lineage "CNS/Brain"` produces ranked biomarker table
3. Delivers the patient stratification evidence needed for investor conversations

### Full Delivery Order

P1 (Biomarker) → P2 (VRK2) → P3 (Docking) → P4 (Co-Crystal) → Polish

Each story is independently releasable after its checkpoint.

---

## Notes

- US2 edits `pipeline/stages/structalign/structalign.py` in place — no new directory
- US3 requires optional pip install (`vina meeko rdkit-pypi`); pipeline must not break without them
- `compounds_filtered.csv` (produced by spec-003 `pipeline fetch`) is the scaffold SMILES source for US3 and US4
- `binding_site_vrk1.csv` (produced by spec-003 `pipeline structalign`) is the binding site residue source for US3 and US4
- Control RMSD > 2.0 Å is a warning, not a hard failure — docking report is still written but flagged as unvalidated
