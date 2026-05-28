# Tasks: Computational Validation Gates

**Input**: Design documents from `/specs/005-validation-gates/`

**Prerequisites**: plan.md ✓, spec.md ✓, research.md ✓, data-model.md ✓, contracts/ ✓

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Install new dependencies and create the validate stage skeleton that all gates share.

- [ ] T001 Install admet-ai: `pip install admet-ai` and verify import in pipeline env
- [ ] T002 Install gmx_mmpbsa + gromacs + ambertools: `conda install -c conda-forge gmx_mmpbsa gromacs ambertools` and verify `gmx_MMPBSA --version`
- [ ] T003 [P] Install openmm + openmmforcefields + openff-toolkit + pdbfixer: `conda install -c conda-forge openmm openmmforcefields openff-toolkit pdbfixer` and verify imports
- [ ] T004 [P] Install runpod: `pip install runpod` and verify import
- [ ] T005 Create `pipeline/stages/validate/__init__.py` (empty)
- [ ] T006 Create `pipeline/stages/validate/gates/__init__.py` (empty)
- [ ] T007 Add `GateStatus` enum and `GateResult` dataclass to `pipeline/models.py` per data-model.md (fields: gate_name, status, score, reason, details, report_path, duration_s, timestamp)
- [ ] T008 Add `ValidationResult` dataclass to `pipeline/models.py` (fields: gene_symbol, scaffold_id, smiles, gates, overall_pass, handoff_ready, created_at)

**Checkpoint**: Dependencies installed, dataclasses defined — gate implementations can begin

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Shared utilities needed by all four gate modules before any gate can be implemented.

- [ ] T009 Implement SMILES salt-stripping helper in `pipeline/stages/validate/validate.py`: strips counterions from multi-fragment SMILES using RDKit, returns largest fragment
- [ ] T010 Implement `_load_smiles(gene_symbol, scaffold_id, settings)` in `pipeline/stages/validate/validate.py`: reads SMILES from `compounds_filtered.csv`, applies salt stripping
- [ ] T011 Implement `_load_docking_pdbqt(gene_symbol, scaffold_id, settings)` in `pipeline/stages/validate/validate.py`: locates `docking_poses_{scaffold}.pdbqt`, raises clear error if not found (directs user to `pipeline dock`)
- [ ] T012 Implement `_cache_gate_result(gene_symbol, scaffold_id, gate_name, result, cache)` and `_load_cached_gate_result(...)` in `pipeline/stages/validate/validate.py` using existing CacheManager with key `validate_{gate_name}_{scaffold_id}`
- [ ] T013 Implement `_write_gate_report(gene_symbol, scaffold_id, gate_name, result, settings)` in `pipeline/stages/validate/validate.py`: writes markdown report to `data/results/{gene}/validate_{gate}_{scaffold}.md`

**Checkpoint**: Shared utilities complete — all four gate modules can now be implemented independently

---

## Phase 3: User Story 1 — ADMET Gate (Priority: P1) 🎯 MVP

**Goal**: `pipeline validate --gate admet` predicts 6 ADMET properties from SMILES and returns PASS/FAIL within 60 seconds.

**Independent Test**: `pipeline validate --target VRK1 --scaffold SCF-009 --gate admet` completes in < 60 s and writes `validate_admet_SCF-009.md` with BBB, CYP, solubility, and bioavailability scores.

### Implementation

- [ ] T014 [US1] Implement `pipeline/stages/validate/gates/admet.py`: `run_admet_gate(gene_symbol, scaffold_id, settings, cache, force, console) -> GateResult`
  - Import `admet_ai.ADMETModel`, instantiate model
  - Call `_load_smiles()` to get SMILES
  - Run `model.predict(smiles=smiles)` → dict of scores
  - Apply thresholds: BBB_Martini > 0.5, CYP1A2/2D6/3A4_Inhibitor < 0.3, Solubility logS > −4, HIA_Hou > 0.3
  - Build `GateResult` with status, score (BBB score as primary), reason listing failing properties, details dict
  - Call `_cache_gate_result()` and `_write_gate_report()`
  - Return `GateResult`
- [ ] T015 [US1] Write ADMET gate markdown report template in `_write_gate_report()`: table of all 6 properties with score, threshold, and PASS/FAIL per property; overall gate decision header
- [ ] T016 [US1] Add unit tests in `tests/unit/test_validate_admet.py`: test pass case (all thresholds met), fail case (BBB fails), salt SMILES handling, cache hit returns same result

**Checkpoint**: `pipeline validate --gate admet` works end-to-end on SCF-009, report written

---

## Phase 4: User Story 2 — MM-GBSA Gate (Priority: P2)

**Goal**: `pipeline validate --gate mmgbsa` rescores the top docking pose with MM-GBSA and returns ΔG (kcal/mol) with PASS/FAIL against −7.0 kcal/mol threshold.

**Independent Test**: `pipeline validate --target VRK1 --scaffold SCF-009 --gate mmgbsa` returns a ΔG value and writes `validate_mmgbsa_SCF-009.md` in under 10 minutes.

### Implementation

- [ ] T017 [US2] Implement `_pdbqt_to_pdb(pdbqt_path, out_pdb_path)` helper in `pipeline/stages/validate/gates/mmgbsa.py`: converts top-pose PDBQT to PDB using RDKit/meeko for gmx_MMPBSA input
- [ ] T018 [US2] Implement `_run_gmx_mmpbsa(receptor_pdb, ligand_pdb, work_dir)` in `pipeline/stages/validate/gates/mmgbsa.py`: prepares topology (GAFF2 + ff19SB), runs `gmx_MMPBSA` subprocess, parses FINAL_RESULTS_MMGBSA.dat for ΔG_bind value; raises RuntimeError (→ ERROR status) on convergence failure
- [ ] T019 [US2] Implement `run_mmgbsa_gate(gene_symbol, scaffold_id, settings, cache, force, console) -> GateResult` in `pipeline/stages/validate/gates/mmgbsa.py`
  - Check docking PDBQT exists (via `_load_docking_pdbqt()`); error if not
  - Create temp work dir under `data/cache/{gene}/mmgbsa_{scaffold}/`
  - Call `_pdbqt_to_pdb()` and `_run_gmx_mmpbsa()`
  - Build GateResult: status PASS if ΔG ≤ −7.0, FAIL otherwise, ERROR on exception
  - Cache and write report
- [ ] T020 [US2] Write MM-GBSA report: ΔG value, threshold, PASS/FAIL, energy components (if available), note on forcefield used
- [ ] T021 [US2] Add unit tests in `tests/unit/test_validate_mmgbsa.py`: test GateResult construction for pass/fail/error cases; mock subprocess call

**Checkpoint**: `pipeline validate --gate mmgbsa` works on SCF-009, ΔG returned

---

## Phase 5: User Story 3 — Selectivity Docking Panel Gate (Priority: P2)

**Goal**: `pipeline validate --gate selectivity` docks the scaffold into VRK2, EGFR, CDK2, PLK1 and computes selectivity index vs. the primary target.

**Independent Test**: `pipeline validate --target VRK1 --scaffold SCF-009 --gate selectivity` produces a selectivity table and SI value and writes `validate_selectivity_SCF-009.md` within 2 hours.

### Implementation

- [ ] T022 [US3] Define `SELECTIVITY_PANEL` constant in `pipeline/stages/validate/gates/selectivity.py`: list of `OffTargetEntry` dataclasses for VRK2 (AF2, O95551), EGFR (1M17, AQ4, A), CDK2 (1E9H, ATP, A), PLK1 (2OKR, ADP, A)
- [ ] T023 [US3] Implement `_fetch_offtarget_structure(entry, panel_dir)` in `pipeline/stages/validate/gates/selectivity.py`: downloads PDB from RCSB or AF2 from EBI if not cached in `data/cache/shared/selectivity_panel/`; adds warning flag for AF2 entries
- [ ] T024 [US3] Implement `_dock_offtarget(smiles, scaffold_id, entry, panel_dir, settings)` in `pipeline/stages/validate/gates/selectivity.py`: reuses `_prepare_receptor()`, `_prepare_ligand()`, `_extract_ligand_centroid()`, `_define_box()`, `_run_vina()` from `pipeline/stages/dock/dock.py`; returns top pose affinity
- [ ] T025 [US3] Implement `run_selectivity_gate(gene_symbol, scaffold_id, settings, cache, force, console) -> GateResult` in `pipeline/stages/validate/gates/selectivity.py`
  - Skip off-targets that match the primary target (e.g. exclude EGFR panel entry when target is EGFR)
  - Fetch structures for all panel entries
  - Dock scaffold into each off-target (sequentially; each ~20–30 min)
  - Compute SI = |primary_dg| / max(|offtarget_dg|)
  - Build GateResult: PASS if SI ≥ 10, FAIL otherwise; details dict includes per-off-target affinity and SI
  - Cache and write report
- [ ] T026 [US3] Write selectivity report: table of off-target affinities vs. primary target, SI calculation, PASS/FAIL, AF2 warning if VRK2 used
- [ ] T027 [US3] Add unit tests in `tests/unit/test_validate_selectivity.py`: SI calculation logic, primary-target exclusion logic, GateResult for pass/fail cases

**Checkpoint**: `pipeline validate --gate selectivity` works on SCF-009, SI value and selectivity table produced

---

## Phase 6: User Story 4 — MD Pose Stability Gate (Priority: P3)

**Goal**: `pipeline validate --gate md` prepares the system locally, submits to RunPod A100, and returns RMSD-based PASS/FAIL within 2 hours at cost ≤ $5.

**Independent Test**: `pipeline validate --target VRK1 --scaffold SCF-009 --gate md` (with `RUNPOD_API_KEY` set) submits a 20 ns job, polls until complete, downloads RMSD CSV, and writes `validate_md_SCF-009.md`.

### Implementation

- [ ] T028 [US4] Implement `_prepare_md_system(receptor_pdb, ligand_pdbqt, work_dir)` in `pipeline/stages/validate/gates/md.py`
  - PDB → PDBFixer (cap termini, add missing residues/atoms, add hydrogens)
  - PDBQT → PDB (using RDKit), parametrise with GAFF2 via `openff-toolkit`
  - OpenMM Modeller: combine protein + ligand, solvate TIP3P 10 Å shell
  - Apply HMR (4 fs timestep) via `openmm.app.HMassRepartitioning`
  - Energy minimise (500 steps max)
  - Serialise system to OpenMM XML: `system.xml`, `topology.pdb`
  - Return dict with paths to XML files and atom count
- [ ] T029 [US4] Implement `_estimate_runpod_cost(atom_count, duration_ns)` in `pipeline/stages/validate/gates/md.py`: estimates cost from community A100 rate (~$1.20/hr) × estimated wall time (atom_count / 500000 ns/day × duration_ns); returns float USD
- [ ] T030 [US4] Implement `_submit_runpod_job(system_files, md_config, api_key)` in `pipeline/stages/validate/gates/md.py`: uses `runpod` SDK to submit serverless job with system XML, 20 ns config, 90-minute timeout; returns job_id
- [ ] T031 [US4] Implement `_poll_runpod_job(job_id, api_key, poll_interval_s=60)` in `pipeline/stages/validate/gates/md.py`: polls job status until COMPLETED/FAILED/TIMEOUT; downloads RMSD CSV and trajectory summary on completion
- [ ] T032 [US4] Implement `_compute_rmsd_pass(rmsd_csv_path)` in `pipeline/stages/validate/gates/md.py`: reads RMSD CSV (time_ns, rmsd_A columns), computes mean RMSD over final 10 ns; returns (mean_rmsd, pass_bool); handles timeout case where < 15 ns available → ERROR
- [ ] T033 [US4] Implement `run_md_gate(gene_symbol, scaffold_id, settings, cache, force, console, md_max_cost) -> GateResult` in `pipeline/stages/validate/gates/md.py`
  - Check `RUNPOD_API_KEY` env var; raise ERROR with helpful message if missing
  - Call `_prepare_md_system()` (~5 min local)
  - Call `_estimate_runpod_cost()`; raise ERROR if > md_max_cost
  - Call `_submit_runpod_job()` and `_poll_runpod_job()`
  - Call `_compute_rmsd_pass()`; build GateResult
  - Save RMSD CSV to `data/results/{gene}/validate_md_{scaffold}_rmsd.csv`
  - Cache and write report
- [ ] T034 [US4] Write MD report: cost estimate, actual cost, atom count, simulation length, mean RMSD, RMSD plot (ASCII sparkline), PASS/FAIL
- [ ] T035 [US4] Add unit tests in `tests/unit/test_validate_md.py`: cost estimation logic, RMSD pass/fail threshold, missing API key error, timeout handling (< 15 ns → ERROR)

**Checkpoint**: MD gate runs end-to-end with real RunPod API key on SCF-009

---

## Phase 7: User Story 5 — Dashboard and Wet-Lab Handoff Report (Priority: P2)

**Goal**: `pipeline validate --dashboard` shows a Rich table of all scaffolds × gates and auto-generates wet-lab handoff reports for all-pass scaffolds.

**Independent Test**: After running gates for SCF-009 and SCF-156, `pipeline validate --target VRK1 --dashboard` prints a 2-row grid and writes `validation_dashboard.md`.

### Implementation

- [ ] T036 [US5] Implement `_load_all_validation_results(gene_symbol, settings, cache)` in `pipeline/stages/validate/validate.py`: scans CacheManager for all `validate_*_{scaffold}` keys for a target; returns dict of scaffold_id → ValidationResult
- [ ] T037 [US5] Implement `_render_dashboard(gene_symbol, results, console)` in `pipeline/stages/validate/validate.py`: builds Rich Table with scaffolds as rows and ADMET/MM-GBSA/Selectivity/MD/Handoff as columns; colour-codes PASS (green), FAIL (red), ERROR (yellow), NOT_RUN (dim)
- [ ] T038 [US5] Implement `_write_dashboard_md(gene_symbol, results, settings)` in `pipeline/stages/validate/validate.py`: writes markdown table + links to individual gate reports to `data/results/{gene}/validation_dashboard.md`; writes machine-readable `validation_dashboard.json`
- [ ] T039 [US5] Implement `_write_wetlab_handoff(gene_symbol, scaffold_id, result, settings)` in `pipeline/stages/validate/validate.py`: generates `wetlab_handoff_{scaffold}.md` with ADMET summary table, MM-GBSA ΔG, selectivity index, MD RMSD, links to all gate reports and existing docking report; only called when `result.handoff_ready == True`
- [ ] T040 [US5] Implement `run_dashboard(gene_symbol, settings, cache, console)` in `pipeline/stages/validate/validate.py`: calls `_load_all_validation_results()`, renders dashboard, writes markdown, triggers handoff report for any `handoff_ready` scaffold

**Checkpoint**: Dashboard renders correctly with mixed PASS/FAIL/NOT_RUN states; handoff report generated for passing scaffolds

---

## Phase 8: CLI Integration & Orchestrator

**Purpose**: Wire all gates into the `pipeline validate` CLI command with correct sequencing and `--all-scaffolds` support.

- [ ] T041 Implement `run_validate(gene_symbol, scaffold_id, gate, settings, cache, force, console, md_max_cost)` orchestrator in `pipeline/stages/validate/validate.py`
  - If `gate` specified: run only that gate
  - If no `gate`: run ADMET → MM-GBSA → Selectivity → MD in sequence; skip MD if any of first 3 = FAIL; skip subsequent gates if any gate = ERROR
  - Persist ValidationResult after each gate completes
  - Trigger handoff report automatically if all gates PASS
- [ ] T042 Add `validate` command to `pipeline/cli.py` following existing command pattern:
  - Options: `--target`/`-t`, `--all`, `--scaffold`, `--all-scaffolds`, `--top-n`, `--gate`, `--dashboard`, `--force`, `--md-max-cost` (default 5.0), `--data-dir`
  - `--dashboard` mode: call `run_dashboard()` and exit
  - Normal mode: resolve scaffold list via `_resolve_scaffolds()` (existing helper), call `run_validate()` per scaffold
  - Exit code 0 on PASS/FAIL; exit code 1 on ERROR
- [ ] T043 Add `Bash(pipeline validate *)` to `.claude/settings.json` permissions allow list

**Checkpoint**: `pipeline validate --target VRK1 --scaffold SCF-009` runs full gate sequence from CLI

---

## Phase 9: Polish & Cross-Cutting Concerns

- [ ] T044 [P] Add `pipeline validate` to `pipeline status` command output: show gate completion counts per target alongside existing fetch/analyze/report stages
- [ ] T045 [P] Update `pipeline rank` command to include ADMET BBB score and MM-GBSA ΔG columns alongside Vina affinity in `scaffold_ranking.md`
- [ ] T046 [P] Update `quickstart.md` with verified working command examples after end-to-end test on SCF-009
- [ ] T047 Run end-to-end validation: `pipeline validate --target VRK1 --scaffold SCF-009` through all gates; confirm output files exist and dashboard renders correctly
- [ ] T048 Run end-to-end validation for SCF-156 in parallel with SCF-009 to confirm dashboard shows 2-row grid correctly

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — install deps and create skeleton first
- **Foundational (Phase 2)**: Depends on Phase 1 — shared utilities needed by all gates
- **US1 ADMET (Phase 3)**: Depends on Phase 2 — can start after foundational complete; no dependency on other gates
- **US2 MM-GBSA (Phase 4)**: Depends on Phase 2 — independent of ADMET
- **US3 Selectivity (Phase 5)**: Depends on Phase 2 — independent of ADMET and MM-GBSA
- **US4 MD (Phase 6)**: Depends on Phase 2 — independent implementation; requires RunPod API key to test end-to-end
- **US5 Dashboard (Phase 7)**: Depends on at least one gate producing results (Phase 3 minimum)
- **CLI Integration (Phase 8)**: Depends on all gate implementations (Phases 3–7)
- **Polish (Phase 9)**: Depends on Phase 8

### User Story Dependencies

- **US1, US2, US3, US4** can all be implemented in parallel after Phase 2
- **US5** can be started after US1 (needs at least one gate to produce a result to test)
- **CLI (Phase 8)** must wait for all gate implementations

### Parallel Opportunities

- T001–T004 (dependency installs) can all run in parallel
- T007–T008 (dataclasses) can run in parallel with T005–T006 (init files)
- T014–T016 (ADMET), T017–T021 (MM-GBSA), T022–T027 (Selectivity), T028–T035 (MD) can all be implemented in parallel after Phase 2
- T044–T046 (polish) can run in parallel

---

## Implementation Strategy

### MVP First (US1 — ADMET gate only)

1. Complete Phase 1 (Setup) + Phase 2 (Foundational)
2. Complete Phase 3 (US1 — ADMET gate)
3. Add US5 dashboard skeleton (T036–T040) — works with ADMET results only
4. Add minimal CLI (T042) for `--gate admet` only
5. **STOP and VALIDATE**: `pipeline validate --target VRK1 --scaffold SCF-009 --gate admet` works end-to-end

This gives a working gate in ~2 days with zero cloud dependencies.

### Incremental Delivery

1. MVP: ADMET gate → validates drug-likeness for all scaffolds in minutes
2. Add MM-GBSA → replaces Vina scores with physics-based ranking
3. Add Selectivity → identifies VRK1-selective scaffolds
4. Add MD → confirms pose stability for handoff-ready candidates
5. Each gate adds value independently; dashboard improves as more gates produce results

---

## Notes

- [P] tasks can run in parallel (different files, no cross-dependencies)
- Gate modules import lazily — if a dependency is not installed, the gate raises a clear RuntimeError at call time rather than at import
- All gate functions follow the same signature: `run_X_gate(gene_symbol, scaffold_id, settings, cache, force, console, **kwargs) -> GateResult`
- SMILES salt-stripping (T009) must be applied before any gate receives the SMILES — enforced in `_load_smiles()`
- The RunPod job script (run on the cloud worker) is a separate artefact; for now, assume a pre-built RunPod template image with GROMACS + OpenMM + MDTraj; T030 should document the template ID
