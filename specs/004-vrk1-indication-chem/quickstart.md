# Quickstart & Validation: VRK1 Indication Validation and Chemistry Deepening

## Prerequisites

```bash
# Confirm spec-003 outputs are present
ls data/results/VRK1/depmap_report.md          # DepMap cache required for biomarker
ls data/cache/VRK1/structures/6AC9.pdb         # Structure required for dock + cocrystal
ls data/results/VRK1/binding_site_comparison.csv  # Required for cocrystal brief

# Install new optional dependencies (docking only)
pip install "vina>=1.2.7" "meeko>=0.5" rdkit-pypi
```

---

## US1 — Biomarker Analysis

### Cold start
```bash
pipeline biomarker --target VRK1 --lineage "CNS/Brain"
```

**Expected console output** (approximately):
```
────────── VRK1 — biomarker ──────────
  VRK1: loading dependency cache (1,078 lines)...
  VRK1: downloading OmicsSomaticMutations.csv (~150 MB)...
  VRK1: analysing CNS/Brain (90 lines: N dependent / M non-dependent)...
  VRK1: NNN genes tested | N significant (p < 0.05) | top: GENE (OR=N.N, p=0.NNN)
  VRK1: biomarker_report.md written
  VRK1: research_report.md updated
```

**Validation checklist**:
- [ ] `data/results/VRK1/biomarker_report.md` exists and contains a ranked mutation table
- [ ] `data/results/VRK1/biomarker_results.csv` exists with columns: `gene`, `lineage`, `odds_ratio`, `p_value`, `n_dependent_with_mut`, `n_nondependent_with_mut`
- [ ] Report header states "CNS/Brain" and the number of dependent vs non-dependent lines
- [ ] At least one gene has `p_value < 0.05` (or report explicitly states "no significant biomarker found")
- [ ] `data/results/VRK1/research_report.md` has been updated with a biomarker section

### Cache hit
```bash
time pipeline biomarker --target VRK1 --lineage "CNS/Brain"
```
**Expected**: `SKIP (cached NNN results)` — completes in < 5 seconds.

### Error handling
```bash
pipeline biomarker --target VRK1 --lineage "InvalidLineage"
```
**Expected**: non-zero exit, error message names the invalid lineage and lists available lineages from the DepMap cache.

---

## US2 — VRK2 Three-Way Comparison

### Cold start
```bash
pipeline structalign --target VRK1 --include-vrk2 --force
```

**Expected additional console lines** (beyond spec-003 output):
```
  VRK1: fetching VRK2 structure...
  VRK1: VRK2 — [crystal/pdb_id OR alphafold/O95551] (N residues excluded if AlphaFold)
  VRK1: three-way classification — N VRK1-specific | N pan-VRK | N conserved
  VRK1: structural_selectivity_report.md updated with VRK2 section
```

**Validation checklist**:
- [ ] `data/results/VRK1/binding_site_comparison.csv` has new columns: `vrk2_aa`, `vrk1_vrk2_diff`, `selectivity_class`
- [ ] `selectivity_class` column contains only: "VRK1-specific", "pan-VRK vs EGFR", "VRK2 vs VRK1+EGFR", "conserved"
- [ ] All 85 KLIFS positions have a `selectivity_class` value (no nulls)
- [ ] `data/results/VRK1/structural_selectivity_report.md` contains a "VRK2 Three-Way Comparison" section
- [ ] Report states the number of VRK1-specific positions explicitly

### Regression check (no VRK2 flag)
```bash
pipeline structalign --target VRK1 --force
```
**Expected**: Output identical to spec-003 — no `vrk2_aa` or `selectivity_class` columns in CSV, no VRK2 section in report.

---

## US3 — SCF-013 Docking

### Preflight (Vina not installed)
```bash
pip uninstall vina -y 2>/dev/null; pipeline dock --target VRK1 --scaffold SCF-013
```
**Expected**: Exit code 1, message names AutoDock Vina as missing dependency with `pip install vina meeko rdkit-pypi` instruction.

### Cold start (Vina installed)
```bash
pip install "vina>=1.2.7" "meeko>=0.5" rdkit-pypi
pipeline dock --target VRK1 --scaffold SCF-013
```

**Expected console output**:
```
────────────── VRK1 — dock ──────────────
  VRK1: checking AutoDock Vina installation... ok
  VRK1: preparing receptor 6AC9...
  VRK1: running positive-control dock (ANP into 6AC9)...
  VRK1: control RMSD: N.NN Å [✓ validated / ⚠ unvalidated]
  VRK1: preparing ligand SCF-013...
  VRK1: docking SCF-013 (exhaustiveness=32)...
  VRK1: top pose: -N.N kcal/mol | contacts N/10 selectivity candidates
  VRK1: docking_report.md written
```

**Validation checklist**:
- [ ] `data/results/VRK1/docking_report.md` exists
- [ ] Report contains a positive-control RMSD value
- [ ] Control RMSD ≤ 2.0 Å (setup validated) — if > 2.0, report warns but does not fail
- [ ] Report lists top-3 SCF-013 poses with affinity (kcal/mol)
- [ ] Report maps contacts to KLIFS selectivity candidates from spec-003
- [ ] `data/results/VRK1/docking_results_SCF-013.csv` exists with columns: `pose_rank`, `affinity_kcal_mol`, `contacted_selectivity_candidates`
- [ ] `data/results/VRK1/docking_poses_SCF-013.pdbqt` exists and is non-empty

### Cache hit
```bash
time pipeline dock --target VRK1 --scaffold SCF-013
```
**Expected**: `SKIP (cached 9 poses)` — completes in < 5 seconds.

### Missing scaffold
```bash
pipeline dock --target VRK1 --scaffold SCF-999
```
**Expected**: Exit code 2, error names the missing scaffold ID and lists available scaffold IDs.

---

## US4 — Co-Crystal Brief

### Cold start
```bash
pipeline cocrystal --target VRK1 --scaffold SCF-013
```

**Expected console output**:
```
───────── VRK1 — cocrystal ──────────
  VRK1: loading structure 6AC9...
  VRK1: retrieving crystallisation conditions from RCSB...
  VRK1: N condition(s) found
  VRK1: analysing SCF-013 for crystal packing compatibility...
  VRK1: cocrystal_brief.md written
```

**Validation checklist**:
- [ ] `data/results/VRK1/cocrystal_brief.md` exists
- [ ] Brief names PDB 6AC9 as reference crystal form
- [ ] Brief cites confirmed conditions: 27.5% PEG 3350, 0.2M ammonium sulfate, 0.1M HEPES pH 7.0
- [ ] Brief specifies minimum resolution for each key subpocket (Gatekeeper, Hinge, DFG)
- [ ] Brief lists recommended minimum overall resolution (≤ 2.2 Å)
- [ ] Brief identifies atoms in SCF-013 that may contact crystal symmetry mates (or states "none flagged")

### Error: structalign not run first
```bash
rm -f data/results/VRK1/binding_site_comparison.csv
pipeline cocrystal --target VRK1 --scaffold SCF-013
```
**Expected**: Exit code 1, error states that structural selectivity analysis must be run first (`pipeline structalign --target VRK1`).

---

## Full Sequence (in order)

```bash
# 1. Biomarker (reads existing DepMap cache)
pipeline biomarker --target VRK1 --lineage "CNS/Brain"

# 2. VRK2 comparison (extends structalign)
pipeline structalign --target VRK1 --include-vrk2 --force

# 3. Docking
pipeline dock --target VRK1 --scaffold SCF-013

# 4. Co-crystal brief
pipeline cocrystal --target VRK1 --scaffold SCF-013

# 5. Verify all output files
ls data/results/VRK1/biomarker_report.md \
   data/results/VRK1/docking_report.md \
   data/results/VRK1/cocrystal_brief.md
```

All four commands should complete without error. Cache-hit re-runs of each should complete in < 5 seconds.
