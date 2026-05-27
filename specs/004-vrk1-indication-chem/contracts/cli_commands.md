# CLI Command Contracts: VRK1 Indication Validation and Chemistry Deepening

All commands follow the established `rxpipeline` CLI pattern: `pipeline <command> [OPTIONS]`.

---

## `pipeline biomarker`

Compute co-mutation enrichment for a target gene's dependency lineage.

```
pipeline biomarker [OPTIONS]

Options:
  --target / -t TEXT        Gene symbol (e.g. VRK1) [required]
  --lineage TEXT            OncotreeLineage to analyse (e.g. "CNS/Brain") [required]
  --all                     Run for all targets in config
  --force                   Re-run even if cached output exists
  --min-lines INT           Min cell lines with mutation for inclusion [default: 3]
  --data-dir PATH           Override data directory
```

**Exit codes**: 0 = success, 1 = target not in DepMap cache (run `pipeline depmap` first), 2 = mutation file download failed.

**Console output pattern**:
```
────────── VRK1 — biomarker ──────────
  VRK1: loading dependency cache (1,078 lines)...
  VRK1: downloading OmicsSomaticMutations.csv (~150 MB)...
  VRK1: analysing CNS/Brain (90 lines: 75 dependent / 15 non-dependent)...
  VRK1: 847 genes tested | 12 significant (p < 0.05) | top: TP53 (OR=4.2, p=0.003)
  VRK1: biomarker_report.md written
  VRK1: research_report.md updated
```

**Cache hit pattern**:
```
  VRK1       biomarker            SKIP  (cached 847 results)
```

**Output files**:
- `data/results/{target}/biomarker_report.md`
- `data/results/{target}/biomarker_results.csv`

---

## `pipeline structalign` (extended)

Existing command; `--include-vrk2` flag behaviour is completed by this spec.

```
pipeline structalign [OPTIONS]

Options:
  --target / -t TEXT        Gene symbol [required]
  --all                     Run for all targets
  --include-vrk2            Add VRK2 three-way comparison
  --force                   Re-run even if cached outputs exist
  --cutoff FLOAT            Binding site distance cutoff in Å [default: 6.0]
  --data-dir PATH           Override data directory
```

**With `--include-vrk2`, additional console lines**:
```
  VRK1: fetching VRK2 structure...
  VRK1: VRK2 — AlphaFold model (3 residues excluded, pLDDT < 70)
  VRK1: three-way classification — 12 VRK1-specific | 8 pan-VRK | 46 conserved
  VRK1: structural_selectivity_report.md updated with VRK2 section
```

**Output files** (additions to spec-003):
- `data/results/{target}/binding_site_comparison.csv` — extended with `vrk2_aa`, `vrk1_vrk2_diff`, `selectivity_class` columns
- `data/results/{target}/structural_selectivity_report.md` — extended with VRK2 section

---

## `pipeline dock`

Dock a scaffold into a target's binding site using AutoDock Vina.

```
pipeline dock [OPTIONS]

Options:
  --target / -t TEXT        Gene symbol (e.g. VRK1) [required]
  --scaffold TEXT           Scaffold ID from compound library (e.g. SCF-013) [required]
  --all-scaffolds           Dock all scaffolds in the compound library
  --force                   Re-run even if cached output exists
  --exhaustiveness INT      Vina exhaustiveness [default: 32]
  --data-dir PATH           Override data directory
```

**Exit codes**: 0 = success (including control RMSD warning), 1 = AutoDock Vina not installed, 2 = scaffold not found in compound library, 3 = no valid poses produced.

**Console output pattern**:
```
────────────── VRK1 — dock ──────────────
  VRK1: checking AutoDock Vina installation... ok (v1.2.7)
  VRK1: preparing receptor 6AC9...
  VRK1: running positive-control dock (ANP into 6AC9)...
  VRK1: control RMSD: 1.43 Å ✓ (validated setup)
  VRK1: preparing ligand SCF-013...
  VRK1: docking SCF-013 (exhaustiveness=32)...
  VRK1: top pose: -8.2 kcal/mol | contacts 7/10 selectivity candidates
  VRK1: docking_report.md written
```

**Vina not installed error**:
```
  ERROR: AutoDock Vina is not installed.
  Install with: pip install vina meeko rdkit-pypi
  Then re-run: pipeline dock --target VRK1 --scaffold SCF-013
```

**Cache hit pattern**:
```
  VRK1       dock (SCF-013)       SKIP  (cached 9 poses)
```

**Output files**:
- `data/results/{target}/docking_report.md`
- `data/results/{target}/docking_results_{scaffold}.csv`
- `data/results/{target}/docking_poses_{scaffold}.pdbqt`

---

## `pipeline cocrystal`

Generate a co-crystallisation experimental brief for a target/scaffold pair.

```
pipeline cocrystal [OPTIONS]

Options:
  --target / -t TEXT        Gene symbol (e.g. VRK1) [required]
  --scaffold TEXT           Scaffold ID (e.g. SCF-013) [required]
  --force                   Re-generate even if brief already exists
  --data-dir PATH           Override data directory
```

**Exit codes**: 0 = success, 1 = no structural selectivity report found (run `pipeline structalign` first).

**Console output pattern**:
```
───────── VRK1 — cocrystal ──────────
  VRK1: loading structure 6AC9 (2.07 Å, space group P 21 21 21)...
  VRK1: retrieving crystallisation conditions from RCSB...
  VRK1: 1 VRK1 condition found (6AC9)
  VRK1: analysing SCF-013 for crystal packing compatibility...
  VRK1: 2 atoms flagged for potential packing contacts
  VRK1: cocrystal_brief.md written
```

**Cache hit pattern**:
```
  VRK1       cocrystal (SCF-013)  SKIP  (brief exists; use --force to regenerate)
```

**Output files**:
- `data/results/{target}/cocrystal_brief.md`
