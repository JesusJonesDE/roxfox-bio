# Tasks: VRK1 Oncology Evidence — DepMap + Structural Alignment

**Input**: Design documents from `specs/003-vrk1-oncology-evidence/`

**Branch**: `003-vrk1-oncology-evidence`

**Stack**: Python 3.11 | Typer + Rich | pandas + httpx + tenacity | biopython≥1.83 | opencadd≥0.8

**Organization**: Tasks grouped by user story — each phase is independently deliverable.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Add new dependencies and create new stage package directories

- [x] T001 Add `biopython>=1.83` and `opencadd>=0.8` to `[project].dependencies` in `pyproject.toml`
- [x] T002 [P] Create `pipeline/stages/depmap/__init__.py` (empty module init)
- [x] T003 [P] Create `pipeline/stages/structalign/__init__.py` (empty module init)

**Checkpoint**: `pip install -e .` completes without errors; `python -c "import Bio; import opencadd"` succeeds in the rxpipeline conda env.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Add shared cache path convention used by both new commands

**⚠️ CRITICAL**: Required before structalign can cache the EGFR reference structure

- [x] T004 Add `shared_structures_dir` property to `Settings` class in `pipeline/config.py` returning `self.data_dir / "cache" / "shared" / "structures"`; this is the cache location for reference structures shared across targets (e.g. EGFR 1M17)

**Checkpoint**: `from pipeline.config import Settings; s = Settings(); print(s.shared_structures_dir)` prints the expected path.

---

## Phase 3: User Story 1 — Cancer Dependency Landscape (Priority: P1) 🎯 MVP

**Goal**: `pipeline depmap --target VRK1` produces a ranked cancer lineage dependency report from DepMap CRISPR Chronos data.

**Independent Test**: Running `pipeline depmap --target VRK1` from a clean state downloads data and writes `data/results/VRK1/depmap_lineage_summary.csv` (≥15 lineage rows) and `data/results/VRK1/depmap_report.md`. Re-running within the cache window prints a SKIP line and completes in < 5 seconds.

### Implementation

- [x] T005 [US1] Implement DepMap manifest fetch in `pipeline/stages/depmap/depmap.py`: GET `https://depmap.org/portal/api/download/files` → parse returned CSV → locate latest `CRISPRGeneEffect.csv` and `Model.csv` rows → return their signed GCS download URLs; use `httpx` with `tenacity` retry (3 attempts, exponential backoff) on 429/5xx
- [x] T006 [P] [US1] Implement `CRISPRGeneEffect.csv` download in `pipeline/stages/depmap/depmap.py`: use `pd.read_csv(url, index_col=0, usecols=lambda c: c == "ModelID" or c.startswith(f"{gene_symbol} ("))` to load only the target gene column; raise `ValueError` with actionable message if target column is absent
- [x] T007 [P] [US1] Implement `Model.csv` download in `pipeline/stages/depmap/depmap.py`: use `pd.read_csv(url, usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "CCLEName"])`; coerce null `OncotreeLineage` to `"Unknown"`
- [x] T008 [US1] Implement merge + `dependency_tier` classification in `pipeline/stages/depmap/depmap.py`: left-join gene effects onto Model on `ModelID`; drop rows where `gene_effect` is null; assign `dependency_tier` enum using thresholds ≤−0.5=`strongly_dependent`, −0.5 to −0.3=`moderately_dependent`, −0.3 to −0.1=`weakly_dependent`, >−0.1=`not_essential` (depends on T006, T007)
- [x] T009 [US1] Implement `LineageSummary` aggregation in `pipeline/stages/depmap/depmap.py`: groupby `OncotreeLineage` → compute `n_lines`, `median_effect`, `mean_effect`, `pct_strongly_dependent`, `n_strongly_dependent`, `dependency_tier` (from median); compute scalar `pan_essential_flag` as `True` if pct strongly dependent across **all** screened lines > 70%; exclude lineages with `n_lines < 3` from ranked output (depends on T008)
- [x] T010 [US1] Implement cache read/write in `pipeline/stages/depmap/depmap.py`: on cache miss — save merged per-cell-line records to `data/cache/<GENE>/depmap_<TIMESTAMP>.json` via `CacheManager`; on cache hit — load from JSON and skip download; add SKIP console line per contract format (depends on T005)
- [x] T011 [US1] Implement output in `pipeline/stages/depmap/depmap.py`: write `data/results/<GENE>/depmap_lineage_summary.csv` (ranked by `median_effect` ascending, columns per `LineageSummary` entity in `data-model.md`); write `data/results/<GENE>/depmap_report.md` (lineage table + pan-essential finding + plain-language top-3 interpretation paragraph); print console header block per contract format (depends on T009, T010)
- [x] T012 [US1] Add `depmap` subcommand to `pipeline/cli.py`: flags `--target/-t`, `--all`, `--force`, `--max-age` (default 30), `--data-dir`; call `run_depmap(target, settings, cache, force)`; exit 0 on success, exit 1 with Rich error on failure (depends on T011)

**Checkpoint**: US1 complete. `pipeline depmap --target VRK1` produces both output files; second run prints SKIP and completes in < 5 seconds.

---

## Phase 4: User Story 2 — Structural Selectivity Explanation (Priority: P1)

**Goal**: `pipeline structalign --target VRK1` produces a residue-level VRK1 vs. EGFR binding site comparison with gatekeeper identification and a selectivity hypothesis.

**Independent Test**: Running `pipeline structalign --target VRK1` writes `data/results/VRK1/binding_site_comparison.csv` (> 30 rows) and `data/results/VRK1/structural_selectivity_report.md`. The report contains "Met131" (VRK1 gatekeeper) and "Thr790" (EGFR gatekeeper). Re-running uses cached PDB files and completes without re-downloading.

### Implementation

- [x] T013 [P] [US2] Implement VRK1 best-structure selection in `pipeline/stages/structalign/structalign.py`: scan `data/cache/VRK1/structures/*.pdb`; for each file parse `REMARK   2 RESOLUTION` line for resolution (Å) and scan `HETATM` records for non-water ligands (`res.id[0].startswith("H_")`); return path + PDB ID of highest-resolution ligand-bound structure (prefer 6AC9); raise informative error if no structures found
- [x] T014 [P] [US2] Implement EGFR reference structure fetch in `pipeline/stages/structalign/structalign.py`: if `settings.shared_structures_dir / "1M17.pdb"` is absent, create the directory and stream-download `https://files.rcsb.org/download/1M17.pdb` with `httpx`; always return the local path (depends on T004)
- [x] T015 [P] [US2] Implement KLIFS pocket residue mapping in `pipeline/stages/structalign/structalign.py` using `opencadd.databases.klifs`: for each kinase/PDB ID call `KLIFS_SESSION.structures.by_structure_pdb_id(pdb_id)` → get `structure_ID` → call `KLIFS_SESSION.pocket.by_structure_klifs_id(structure_id)` → extract 85-position table with `klifs_position`, `amino_acid`, `pdb_residue_id`, `subpocket`; tag `is_gatekeeper=True` for position 45, `is_hinge=True` for positions 46–48 (depends on T013, T014)
- [x] T016 [P] [US2] Implement BioPython NeighborSearch binding site extraction in `pipeline/stages/structalign/structalign.py`: load PDB with `Bio.PDB.PDBParser`; for each structure collect ligand heavy atoms (HETATM, `res.id[0].startswith("H_")`, exclude water `"W"`); select only conformer A atoms (`atom.get_altloc() in ("", "A")`); run `Bio.PDB.NeighborSearch(all_atoms).search(ligand_center, cutoff, level="R")` to get residue list; store `distance_to_ligand_A` as min distance from residue heavy atoms to any ligand atom (depends on T013, T014)
- [x] T017 [P] [US2] Implement BioPython Superimposer alignment in `pipeline/stages/structalign/structalign.py`: take intersection of VRK1 and EGFR KLIFS positions with valid residues in both; build equal-length `fixed_atoms` (VRK1 Cα) and `moving_atoms` (EGFR Cα) lists; call `Bio.PDB.Superimposer().set_atoms(fixed_atoms, moving_atoms)`; record RMSD and count of shared Cα; write alignment quality to report (depends on T015, T016)
- [x] T018 [P] [US2] Implement `BindingSiteComparison` table in `pipeline/stages/structalign/structalign.py`: for each of 85 KLIFS positions produce row with `vrk1_aa`, `egfr_aa`, `identical`, `difference_type` (classify using rules from `data-model.md`: identical/conservative/steric/charge/h_bond/gap), `selectivity_candidate` (True if difference_type ≠ identical AND subpocket in `{"P-loop","Hinge","Gatekeeper","DFG"}`); write `data/results/<GENE>/binding_site_vrk1.csv`, `binding_site_egfr.csv`, `binding_site_comparison.csv` (depends on T015, T016)
- [x] T019 [US2] Implement `structural_selectivity_report.md` generation in `pipeline/stages/structalign/structalign.py`: report sections — best structure selected + criteria; EGFR reference; RMSD + shared Cα count; gatekeeper line (KLIFS position 45 residue for each kinase); total positions differing + count of selectivity candidates; auto-generated selectivity hypothesis paragraph identifying which VRK1-specific residues SCF-013 likely engages; print console output per contract format (depends on T017, T018)
- [x] T020 [US2] Add `structalign` subcommand to `pipeline/cli.py`: flags `--target/-t`, `--all`, `--include-vrk2` (default False), `--force`, `--cutoff` (default 6.0), `--data-dir`; call `run_structalign(target, settings, cache, force, include_vrk2, cutoff)`; exit 0 on success, exit 1 with Rich error on failure (depends on T019)

**Checkpoint**: US2 complete. `pipeline structalign --target VRK1` produces all 4 output files; report contains Met131 and Thr790 gatekeepers; second run uses cached PDB files.

---

## Phase 5: User Story 3 — VRK1 vs. VRK2 Selectivity Profile (Priority: P2)

**Goal**: `pipeline structalign --target VRK1 --include-vrk2` adds a VRK2 column to the binding site comparison, enabling a three-way analysis of VRK1/VRK2/EGFR selectivity.

**Independent Test**: Running with `--include-vrk2` produces a `binding_site_comparison.csv` containing a `vrk2_aa` column. The report includes a "VRK2 comparison" section. If no PDB structure exists, the AlphaFold model is downloaded and low-confidence positions (pLDDT < 70) are excluded with a warning.

### Implementation

- [x] T021 [US3] Implement VRK2 structure fetch in `pipeline/stages/structalign/structalign.py`: query KLIFS for VRK2 structures; if none found, download AlphaFold model for UniProt O95551 from `https://alphafold.ebi.ac.uk/api/prediction/O95551` → extract download URL → stream PDB to `settings.shared_structures_dir / "VRK2_AF.pdb"`; record model type (`crystal` or `alphafold`) and source in output
- [x] T022 [US3] Implement pLDDT confidence filtering in `pipeline/stages/structalign/structalign.py`: for AlphaFold models, B-factor field = pLDDT score; for each binding site residue check `min(atom.bfactor for atom in residue) >= 70`; exclude residues failing this threshold from comparison table with `vrk2_aa = "LOW_CONF"`; print count of excluded positions as a console warning (depends on T021)
- [x] T023 [US3] Extend `BindingSiteComparison` table with `vrk2_aa` column in `pipeline/stages/structalign/structalign.py`: rerun KLIFS mapping for VRK2; add `vrk2_aa` (1-letter code, `"GAP"`, `"LOW_CONF"`, or `"N/A"` if VRK2 not run) to each row; add "Three-Way Selectivity" section to `structural_selectivity_report.md` noting positions where VRK1 and VRK2 differ (potential VRK1-selective handle) (depends on T022)
- [x] T024 [US3] Wire `--include-vrk2` flag through `run_structalign()` call in `pipeline/cli.py` to activate VRK2 branch; update console output to show VRK2 structure source line when flag is active (depends on T020, T023)

**Checkpoint**: US3 complete. `pipeline structalign --target VRK1 --include-vrk2` adds the VRK2 column; report includes three-way section.

---

## Phase 6: User Story 4 — Integrated Oncology Report Update (Priority: P2)

**Goal**: After both analyses complete, `data/results/research_report.md` is updated to include DepMap lineage rankings and the structural selectivity interpretation in the VRK1 section.

**Independent Test**: After running both `pipeline depmap --target VRK1` and `pipeline structalign --target VRK1`, re-running the existing `pipeline report --target VRK1` (or both new commands) updates `data/results/research_report.md` so `grep "DepMap" data/results/research_report.md` returns a result and `grep "Met131" data/results/research_report.md` returns a result.

### Implementation

- [x] T025 [US4] Implement DepMap findings injection in `pipeline/stages/depmap/depmap.py`: after writing `depmap_report.md`, read `data/results/research_report.md`; locate the VRK1 oncology section (search for `## VRK1` + `### Oncology` markers or equivalent); replace or append a "DepMap Cancer Dependency" subsection with the top-3 lineage table and pan-essential finding; write the updated report back (depends on T011)
- [x] T026 [US4] Implement structural findings injection in `pipeline/stages/structalign/structalign.py`: after writing `structural_selectivity_report.md`, read `data/results/research_report.md`; locate the VRK1 chemistry strategy section; replace or append a "Structural Selectivity (VRK1 vs EGFR)" subsection with gatekeeper comparison and selectivity hypothesis paragraph; write the updated report back (depends on T019)
- [x] T027 [US4] Add report update calls to `pipeline/cli.py`: call the report injection functions at the end of both `depmap` and `structalign` command handlers, after all output files are written successfully; print a single status line `  VRK1: research_report.md updated` (depends on T012, T020, T025, T026)

**Checkpoint**: US4 complete. Running either new command updates the research report. Both injection markers are present in the output.

---

## Final Phase: Polish & Cross-Cutting Concerns

- [x] T028 [P] Add unit tests for `dependency_tier` boundary values (exactly −0.5, −0.3, −0.1, 0.0, values 0.0001 above/below each threshold) and `pan_essential_flag` threshold (70% boundary) in `tests/unit/test_depmap.py`
- [x] T029 [P] Add unit tests for `BindingSiteComparison.difference_type` classification covering all 6 types (identical/conservative/steric/charge/h_bond/gap) and `selectivity_candidate` flagging logic in `tests/unit/test_structalign.py`
- [x] T030 Run quickstart.md end-to-end validation: cold start for both commands, verify output files and console format match contracts, then verify cache-hit runs complete in < 5 seconds

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Setup (Phase 1) — BLOCKS structalign shared cache path
- **US1 (Phase 3)**: Depends on Foundational — independent of US2
- **US2 (Phase 4)**: Depends on Foundational — independent of US1; requires VRK1 PDB cache already populated by prior `pipeline fetch --target VRK1`
- **US3 (Phase 5)**: Depends on US2 (Phase 4) complete — extends structalign module
- **US4 (Phase 6)**: Depends on US1 (Phase 3) AND US2 (Phase 4) complete — injects findings from both
- **Polish (Final Phase)**: Depends on all desired stories complete

### User Story Dependencies

- **US1 (P1)**: Can start after Phase 2 — no dependency on US2
- **US2 (P1)**: Can start after Phase 2 — no dependency on US1; requires pre-existing PDB cache from `pipeline fetch`
- **US3 (P2)**: Depends on US2 — extends the structalign module
- **US4 (P2)**: Depends on US1 + US2 — reads their output files

### Within US1

T005 → (T006 [P], T007 [P]) → T008 → T009 → (T010, T011 in parallel within same module) → T012

### Within US2

(T013 [P], T014 [P]) → (T015 [P], T016 [P]) → (T017 [P], T018 [P]) → T019 → T020

---

## Parallel Example: User Story 2 (US2)

US2 has a diamond dependency pattern with three parallel layers:

```bash
# Layer 1 — fetch structures (run together):
Task: "Scan VRK1 PDB cache for best structure"          # T013
Task: "Fetch EGFR 1M17 reference to shared cache"       # T014

# Layer 2 — map + extract (run together once T013+T014 done):
Task: "KLIFS pocket residue mapping for VRK1 + EGFR"   # T015
Task: "BioPython NeighborSearch binding site extraction" # T016

# Layer 3 — align + compare (run together once T015+T016 done):
Task: "Superimposer structural alignment + RMSD"         # T017
Task: "BindingSiteComparison table + difference_type"    # T018

# Sequential finish:
Task: "Generate structural_selectivity_report.md"        # T019
Task: "Add structalign command to cli.py"                # T020
```

---

## Implementation Strategy

### MVP First (US1 Only — DepMap analysis)

1. Complete Phase 1: Setup (T001–T003)
2. Complete Phase 2: Foundational (T004)
3. Complete Phase 3: US1 — DepMap (T005–T012)
4. **STOP and VALIDATE**: `pipeline depmap --target VRK1` produces ranked lineage report
5. Ship MVP if cancer indication hypothesis is confirmed

### Incremental Delivery

1. Setup + Foundational → Environment ready
2. US1 complete → DepMap cancer dependency answer (MVP)
3. US2 complete → Structural selectivity answer (completes the oncology thesis)
4. US3 complete → VRK1 vs VRK2 selectivity profile (IP strategy input)
5. US4 complete → Consolidated research report updated for funding discussions

### Parallel Developer Strategy

After Phase 2:
- Developer A: US1 (depmap module + CLI command)
- Developer B: US2 (structalign module + CLI command)
Both are in different files with no shared code between them.

---

## Notes

- [P] tasks = different files, no dependencies on each other
- US1 and US2 are independently deliverable P1 stories — sequence is preference, not constraint
- US2 requires `data/cache/VRK1/structures/6AC9.pdb` to already exist (run `pipeline fetch --target VRK1` first)
- Both modules reuse the existing `CacheManager` from `pipeline/cache.py` — no changes to the cache layer
- US3 and US4 are P2 (enhancement) — skip these for MVP
