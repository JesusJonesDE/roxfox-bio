# Tasks: IGHMBP2 Fragment-Based Virtual Screening

**Input**: Design documents from `/specs/006-ighmbp2-fragment-screen/`

**Prerequisites**: plan.md ✓, spec.md ✓, research.md ✓, data-model.md ✓, contracts/ ✓

---

## Phase 1: Setup

**Purpose**: Install new dependency and create the fragment stage skeleton.

- [X] T001 Install `sa-score` package: `pip install sa-score` and verify `import sascorer`
- [X] T002 Create `pipeline/stages/fragment/__init__.py` (empty)
- [X] T003 Commit bundled fallback fragment library to `data/cache/shared/fragment_library/fragments_fallback.smi` — 500 Ro3-compliant SMILES from ZINC (generate from a fixed seed using RDKit random walk or use a known public set)
- [X] T004 Create `data/cache/shared/fragment_library/building_blocks.smi` — 20 drug-like BRICS building blocks for fragment growing (common amine, acid, heterocycle building blocks)

**Checkpoint**: Dependency installed, skeleton created, fallback data committed

---

## Phase 2: Foundational

**Purpose**: Shared orchestrator and state tracking used by all steps.

- [X] T005 Implement `FragmentState` dataclass and `_load_state(gene, settings)` / `_save_state(state, settings)` in `pipeline/stages/fragment/fragment.py` — reads/writes `data/cache/{gene}/fragment_state.json` tracking which steps are complete
- [X] T006 Implement `run_fragment(gene_symbol, step, settings, cache, force, top_n, exhaustiveness, library_size, console)` orchestrator in `pipeline/stages/fragment/fragment.py` — calls step functions in order (pocket → library → dock → cluster → grow → admet), skips completed steps unless force, stops on failure with clear message
- [X] T007 Implement `_write_fragment_report(gene_symbol, state, pocket, hits, clusters, candidates, settings)` in `pipeline/stages/fragment/fragment.py` — writes `data/results/{gene}/fragment_screen_report.md` with pocket summary, top-10 hits table, cluster summary, grown candidates, ADMET results, next-steps section

**Checkpoint**: Orchestrator wires all steps; `--step` and `--force` flags work

---

## Phase 3: User Story 1 — Pocket Identification (Priority: P1) 🎯

**Goal**: `pipeline fragment --target IGHMBP2 --step pocket` identifies the binding pocket in the AlphaFold2 structure in under 60 seconds.

**Independent Test**: Produces `data/results/IGHMBP2/pocket_analysis.json` with volume > 200 Å³ and valid centroid coordinates.

- [X] T008 [US1] Implement `run_pocket(gene_symbol, settings, cache, force, console) -> dict` in `pipeline/stages/fragment/pocket.py`:
  - Locate AF2 PDB: `data/cache/{gene}/alphafold_*.json` → extract PDB URL → check if PDB already cached in `data/cache/{gene}/structures/`; if not, download from EBI (`https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb`) to that path
  - Run fpocket via subprocess: `fpocket -f {pdb_path}`; parse output dir `{stem}_out/{stem}_info.txt` with pandas (whitespace-separated)
  - Select top pocket: sort by Druggability_Score descending, filter Volume > 200 Å³; if none found, relax to > 150 Å³ with warning
  - Extract centroid: parse `{stem}_out/pockets/pocket{N}_atm.pdb`, compute mean of all heavy-atom X/Y/Z coordinates
  - Write `data/results/{gene}/pocket_analysis.json`: `{pocket_id, score, druggability_score, volume_A3, centroid_x, centroid_y, centroid_z, box_size_A: 20.0, plddt_mean}`
  - Cache result with key `"fragment_pocket"` via CacheManager
  - Return pocket dict

- [X] T009 [US1] Add unit tests in `tests/unit/test_fragment_pocket.py`:
  - Mock fpocket subprocess returning a sample `_info.txt`; verify correct pocket selected
  - Pocket with volume < 200 Å³ triggers threshold relaxation warning
  - Correct centroid computed from mock pocket PDB

**Checkpoint**: `pipeline fragment --target IGHMBP2 --step pocket` writes `pocket_analysis.json`

---

## Phase 4: User Story 2 — Fragment Library (Priority: P1)

**Goal**: `pipeline fragment --target IGHMBP2 --step library` downloads and prepares ≥5,000 Ro3-filtered fragments in under 10 minutes.

**Independent Test**: `data/cache/shared/fragment_library/fragments_ro3.smi` exists with ≥5,000 SMILES, all passing MW ≤ 250, HBD ≤ 3, HBA ≤ 3, logP ≤ 3, RotB ≤ 3.

- [X] T010 [US2] Implement `run_library(library_size, settings, cache, force, console) -> Path` in `pipeline/stages/fragment/library.py`:
  - Check if `data/cache/shared/fragment_library/fragments_ro3.smi` exists and is cached; if so, return path (skip download)
  - Download: fetch 5 ZINC22 tranche SMILES files from `https://files.docking.org/2D/{AA}/{AAAA}.smi` using httpx + tenacity; parse line by line (format: `SMILES zinc_id`)
  - Apply Ro3 filter using RDKit: MW ≤ 250, HBD ≤ 3, HBA ≤ 3, logP ≤ 3, RotB ≤ 3; skip invalid SMILES
  - Strip salts: RDKit `SaltRemover`, keep largest fragment
  - Murcko scaffold deduplication: `MurckoScaffold.GetScaffoldForMol`; keep first representative per scaffold
  - Random sample to `library_size` (seed=42 for reproducibility)
  - Write to `data/cache/shared/fragment_library/fragments_ro3.smi` (tab-separated: `SMILES\tfragment_id`)
  - If ZINC download fails: load `fragments_fallback.smi` with warning; log that fallback was used in state
  - Return path to SMILES file

- [X] T011 [US2] Add unit tests in `tests/unit/test_fragment_library.py`:
  - Ro3 filter: MW=260 fails, MW=200 passes
  - Salt stripping: `Cl.CCO` → `CCO`
  - Murcko dedup: two compounds with same scaffold → 1 retained
  - Fallback triggered when download mock raises exception

**Checkpoint**: Fragment library file ready; Ro3 compliance verified

---

## Phase 5: User Story 3 — Fragment Docking Screen (Priority: P1)

**Goal**: `pipeline fragment --target IGHMBP2 --step dock` docks all library fragments against the IGHMBP2 pocket, returning top-50 hits.

**Independent Test**: `data/results/IGHMBP2/fragment_hits.csv` exists with ≥50 rows, affinity column populated, cache-based resumption works on re-run.

- [X] T012 [US3] Implement `run_screen(gene_symbol, library_path, pocket, top_n, exhaustiveness, settings, cache, force, console) -> pd.DataFrame` in `pipeline/stages/fragment/screen.py`:
  - Reuse from `pipeline/stages/dock/dock.py`: `_prepare_receptor`, `_prepare_ligand`, `_extract_ligand_centroid` (not needed — use pocket centroid), `_define_box`, `_run_vina`
  - Prepare receptor PDBQT once (reuse dock cache `data/cache/{gene}/dock/{pdb_stem}_receptor.pdbqt`)
  - Box: use pocket centroid as center, box_size=20.0 Å (from pocket dict)
  - For each fragment in library_path:
    - Check cache key `fragment_dock_{fragment_id}`; if cached, load and skip
    - Call `_prepare_ligand(smiles, fragment_id, dock_cache)` — catch RuntimeError (meeko failure) → log and skip
    - Call `_run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness, output_pdbqt)` → catch Exception → log and skip
    - Cache result `{fragment_id, smiles, affinity_kcal_mol, n_poses}`
  - Report progress every 100 fragments: `{N}/{total} fragments docked, {failed} failed`
  - Sort all results by affinity, take top_n; write `data/results/{gene}/fragment_hits.csv`
  - Return DataFrame of top hits

- [X] T013 [US3] Add unit tests in `tests/unit/test_fragment_screen.py`:
  - Cache hit returns stored result without calling vina
  - meeko failure on one fragment → fragment skipped, screen continues
  - Top-N selection: given 100 docked fragments, returns top 50 by affinity
  - Progress reporting outputs correct counts

**Checkpoint**: Fragment screen runs, cache-based resumption works, top-50 CSV written

---

## Phase 6: User Story 4 — Fragment Clustering (Priority: P2)

**Goal**: `pipeline fragment --target IGHMBP2 --step cluster` groups top-50 hits into ≥3 distinct clusters.

**Independent Test**: `data/results/IGHMBP2/fragment_clusters.csv` exists; each row has cluster_id; ≥3 distinct cluster_ids present.

- [X] T014 [US4] Implement `run_cluster(gene_symbol, hits_df, settings, cache, force, console) -> pd.DataFrame` in `pipeline/stages/fragment/cluster.py`:
  - Compute Morgan ECFP4 fingerprints: `AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)`
  - Compute pairwise Tanimoto distance matrix via `DataStructs.BulkTanimotoSimilarity`
  - Run Butina clustering: `Butina.ClusterData(dists, n, cutoff=0.4, isDistData=True)`
  - Assign cluster_id and is_representative (best affinity per cluster)
  - Write `data/results/{gene}/fragment_clusters.csv`: fragment_id, smiles, affinity_kcal_mol, cluster_id, is_representative
  - Cache result with key `"fragment_cluster"`
  - Return annotated DataFrame

- [X] T015 [US4] Add unit tests in `tests/unit/test_fragment_cluster.py`:
  - 10 identical SMILES → 1 cluster; 10 very different SMILES → 10 clusters
  - Representative = fragment with best (most negative) affinity in cluster
  - Handles single-fragment input (1 cluster, is_representative=True)

**Checkpoint**: `fragment_clusters.csv` written; cluster diversity verified

---

## Phase 7: User Story 5 — Fragment Growing (Priority: P2)

**Goal**: `pipeline fragment --target IGHMBP2 --step grow` grows cluster representatives into ≥20 drug-like candidates (MW 300–450).

**Independent Test**: `data/results/IGHMBP2/grown_candidates.csv` exists with ≥20 rows; all rows have MW 300–450, passes_ro5=True, sa_score < 4.

- [X] T016 [US5] Implement `_load_building_blocks(settings) -> list[Mol]` in `pipeline/stages/fragment/grow.py`: reads `data/cache/shared/fragment_library/building_blocks.smi`, returns list of RDKit Mol objects
- [X] T017 [US5] Implement `_grow_brics(fragment_mol, building_blocks) -> list[Mol]` in `pipeline/stages/fragment/grow.py`:
  - Use `BRICS.BRICSBuild([fragment_mol, bb])` for each building block
  - Collect all products, deduplicate by canonical SMILES
  - Return list of RDKit Mol objects

- [X] T018 [US5] Implement `_grow_smarts(fragment_mol) -> list[Mol]` in `pipeline/stages/fragment/grow.py`:
  - Apply 15 SMARTS transformations (amide, sulfonamide, N-methyl, N-ethyl, F-substitution, CF3, ring closure 5/6, ether, thioether, methylenation, Boc-deprotection, Cbz-deprotection, hydroxyl-to-OMe, NH2-to-NHCH3, introduction of pyridine ring)
  - Run each reaction via `AllChem.ReactionFromSmarts(smarts).RunReactants((mol,))`
  - Collect all products, deduplicate

- [X] T019 [US5] Implement `_filter_grown(mols) -> list[tuple[Mol, dict]]` in `pipeline/stages/fragment/grow.py`:
  - Apply filters: MW 300–450, Ro5 (0 violations), RotB ≤ 8
  - Compute SA score: `sascorer.calculateScore(mol)` → keep SA < 4
  - Return list of (Mol, {mw, logp, hbd, hba, rotb, sa_score})

- [X] T020 [US5] Implement `run_grow(gene_symbol, clusters_df, settings, cache, force, console) -> pd.DataFrame` in `pipeline/stages/fragment/grow.py`:
  - For each cluster representative: run `_grow_brics` + `_grow_smarts`, combine, deduplicate, apply `_filter_grown`
  - Keep top 3 by MW-adjusted score (heavier = more complete) per representative
  - Assign candidate IDs: `IGHMBP2-SCF-{N:03d}`
  - Write `data/results/{gene}/grown_candidates.csv`: candidate_id, smiles, parent_fragment_id, molecular_weight, logp, hbd, hba, rotatable_bonds, ro5_violations, passes_ro5, sa_score, grow_method
  - Cache result with key `"fragment_grow"`
  - If <20 candidates produced: warn; use top representatives directly as candidates

- [X] T021 [US5] Add unit tests in `tests/unit/test_fragment_grow.py`:
  - `_filter_grown`: MW=600 rejected, MW=350 + SA=3 accepted
  - `_grow_smarts`: benzene fragment → at least 3 valid products
  - `run_grow` with 0 cluster reps → returns empty DataFrame without error
  - SA score threshold: SA=5 rejected, SA=3.9 accepted

**Checkpoint**: `grown_candidates.csv` with ≥20 MW 300–450 candidates

---

## Phase 8: User Story 6 — ADMET + Final Output (Priority: P2)

**Goal**: `pipeline fragment --target IGHMBP2 --step admet` produces `compounds_filtered.csv` in VRK1 schema from grown candidates.

**Independent Test**: `compounds_filtered.csv` exists with correct column names matching VRK1 schema; `pipeline dock --target IGHMBP2 --all-scaffolds` runs without schema errors.

- [X] T022 [US6] Implement `run_output(gene_symbol, candidates_df, settings, cache, force, console)` in `pipeline/stages/fragment/output.py`:
  - For each grown candidate: call `run_admet_gate` from `pipeline/stages/validate/gates/admet.py` with relaxed BBB threshold 0.3 (pass a custom threshold dict override)
  - Build `compounds_filtered.csv` in VRK1 schema: `molecule_id=candidate_id, smiles, best_value_nm=None, best_assay_type="fragment_screen_predicted", molecular_weight, logp, hbd, hba, rotatable_bonds, ro5_violations, passes_ro5, scaffold_id=candidate_id, source="fragment_virtual_screen", off_target_flags=0, selectivity_flag=False`
  - Sort by composite rank: ADMET BBB score descending (higher BBB = better CNS)
  - Write `data/results/{gene}/compounds_filtered.csv`
  - Cache result with key `"fragment_output"`

- [X] T023 [US6] Implement `_write_fragment_report()` final version in `pipeline/stages/fragment/fragment.py` — fill in full report with: pocket details table, top-10 hits table, cluster summary table, top-10 grown candidates table (with SA score and ADMET BBB), next steps section pointing to `pipeline dock --target IGHMBP2`

- [X] T024 [US6] Add unit tests in `tests/unit/test_fragment_output.py`:
  - Output CSV has exact column names matching VRK1 schema (test against VRK1 CSV header)
  - BBB threshold override: 0.3 used instead of default 0.5
  - Empty candidates_df → empty CSV with header only, no error

**Checkpoint**: `compounds_filtered.csv` written; `pipeline dock --target IGHMBP2 --all-scaffolds` resolves compound IDs correctly

---

## Phase 9: User Story 7 — CLI Command (Priority: P1)

**Goal**: `pipeline fragment --target IGHMBP2` runs the complete pipeline end-to-end.

**Independent Test**: Single command runs all 6 steps in sequence; `--step` flag runs only the specified step; `--force` clears cache and re-runs.

- [X] T025 [US7] Add `fragment` command to `pipeline/cli.py` following the existing command pattern:
  ```python
  @app.command()
  def fragment(
      target: Optional[str] = typer.Option(None, "--target", "-t"),
      all_targets: bool = typer.Option(False, "--all"),
      step: Optional[str] = typer.Option(None, "--step"),
      force: bool = typer.Option(False, "--force"),
      top_n: int = typer.Option(50, "--top-n"),
      exhaustiveness: int = typer.Option(4, "--exhaustiveness"),
      library_size: int = typer.Option(8000, "--library-size"),
      data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
  ):
  ```
  Import and call `run_fragment()` from `pipeline.stages.fragment.fragment`

- [X] T026 [US7] Add `"Bash(pipeline fragment *)"` to `.claude/settings.local.json` permissions allow list

**Checkpoint**: `pipeline fragment --help` shows all options; `--step pocket` runs only pocket step

---

## Phase 10: Polish & Cross-Cutting Concerns

- [ ] T027 [P] Run end-to-end test: `pipeline fragment --target IGHMBP2 --step pocket` → verify `pocket_analysis.json` created with real IGHMBP2 AF2 structure
- [ ] T028 [P] Run end-to-end test: `pipeline fragment --target IGHMBP2 --step library --library-size 100` → verify `fragments_ro3.smi` created with 100 compounds
- [ ] T029 [P] Update `pipeline status` to include fragment screen progress (n_fragments_docked, n_candidates) per target
- [ ] T030 Run full pipeline integration test on 100-fragment subset: `pipeline fragment --target IGHMBP2 --library-size 100 --exhaustiveness 4` → all steps complete, `compounds_filtered.csv` produced

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies
- **Foundational (Phase 2)**: Depends on Phase 1
- **US1 Pocket (Phase 3)**: Depends on Phase 2; also requires AF2 structure in cache (from `pipeline fetch`)
- **US2 Library (Phase 4)**: Depends on Phase 2; independent of US1 (can run in parallel)
- **US3 Screen (Phase 5)**: Depends on US1 (pocket centroid) AND US2 (library)
- **US4 Cluster (Phase 6)**: Depends on US3
- **US5 Grow (Phase 7)**: Depends on US4
- **US6 Output (Phase 8)**: Depends on US5
- **US7 CLI (Phase 9)**: Depends on all gate implementations being complete (US1–US6)
- **Polish (Phase 10)**: Depends on Phase 9

### Parallel Opportunities

- T008 (pocket.py) and T010 (library.py) can be implemented in parallel (different files)
- T009, T011, T013, T015, T021, T024 (unit tests) can be written in parallel with their corresponding implementations

---

## Implementation Strategy

### MVP First (Phases 1–3 + US7 CLI)

1. Complete Phase 1 (install sa-score)
2. Complete Phase 2 (orchestrator skeleton)
3. Complete Phase 3 (pocket step) + Phase 4 (library step)
4. Complete Phase 5 (dock screen) — the core scientific value
5. Add minimal CLI (T025)
6. **STOP and VALIDATE**: Run `pipeline fragment --target IGHMBP2 --library-size 100 --step pocket`, `--step library`, `--step dock` to confirm end-to-end flow

This MVP delivers fragment hits without growing, useful standalone as a prioritisation tool.

### Full delivery

Add clustering (Phase 6), growing (Phase 7), ADMET output (Phase 8) in sequence.

---

## Notes

- [P] tasks can run in parallel
- Fragment docking is the longest step (~8–12h for 8k fragments) — all other steps are minutes
- The fallback fragment library (T003) must be committed to the repo so the pipeline can run offline
- Building blocks file (T004) must contain chemically valid, commercially available fragments — use Enamine REAL building blocks subset
- sascorer.calculateScore returns float 1.0–10.0; lower = easier to synthesise
- BBB threshold relaxed to 0.3 for IGHMBP2 (neurological disease context) — document clearly in report
