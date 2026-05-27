# Data Model: VRK1 Indication Validation and Chemistry Deepening

## Entities

### BiomarkerResult
One record per gene × lineage combination that passes the minimum cell-count filter.

| Field | Type | Description |
|-------|------|-------------|
| `gene` | str | Hugo gene symbol (e.g. "TP53") |
| `lineage` | str | OncotreeLineage (e.g. "CNS/Brain") |
| `variant_classification` | str | Mutation type (e.g. "Missense_Mutation") or "any_deleterious" for gene-level rollup |
| `n_dependent_with_mut` | int | Strongly-dependent lines with this mutation |
| `n_dependent_without_mut` | int | Strongly-dependent lines without this mutation |
| `n_nondependent_with_mut` | int | Non-dependent lines with this mutation |
| `n_nondependent_without_mut` | int | Non-dependent lines without this mutation |
| `odds_ratio` | float | Sample odds ratio from Fisher's exact test |
| `p_value` | float | Two-sided p-value |
| `enrichment_direction` | str | "enriched_in_dependent" / "depleted_in_dependent" / "ns" |

**Validation rules**:
- Minimum 3 dependent lines with mutation for inclusion (same threshold as lineage filter in spec-003)
- `p_value` < 0.05 for "significant" label in report header; full table always written regardless

**Relationships**: Joins to spec-003 DepMap cache on `ModelID`; source mutation data from `OmicsSomaticMutations.csv` keyed on `ModelID` + `Hugo_Symbol`.

---

### ThreeWayComparison
Extends spec-003's `BindingSiteComparison` with VRK2 amino acid and three-way selectivity class.

| Field | Type | Description |
|-------|------|-------------|
| `klifs_position` | int | KLIFS pocket position 1–85 |
| `subpocket` | str | KLIFS subpocket name |
| `vrk1_aa` | str | VRK1 amino acid one-letter code (or "GAP") |
| `vrk2_aa` | str | VRK2 amino acid (or "GAP" / "LOW_CONF") |
| `egfr_aa` | str | EGFR amino acid (or "GAP") |
| `vrk1_vrk2_diff` | str | difference_type for VRK1 vs VRK2 |
| `vrk1_egfr_diff` | str | difference_type for VRK1 vs EGFR (from spec-003) |
| `selectivity_class` | str | "VRK1-specific" / "pan-VRK vs EGFR" / "VRK2 vs VRK1+EGFR" / "conserved" |
| `is_gatekeeper` | bool | Position == 45 |
| `is_hinge` | bool | Position in {46, 47, 48} |
| `vrk2_low_confidence` | bool | True if VRK2 residue flagged by pLDDT filter |

**Relationships**: Supersedes spec-003's `binding_site_comparison.csv` when `--include-vrk2` is passed; otherwise the two-way file is written unchanged.

---

### DockingPose
One record per pose per ligand docking run.

| Field | Type | Description |
|-------|------|-------------|
| `run_id` | str | "control" (native ligand) or scaffold ID (e.g. "SCF-013") |
| `pose_rank` | int | 1-indexed rank by Vina affinity score |
| `affinity_kcal_mol` | float | Vina binding affinity (negative = more favourable) |
| `rmsd_lb` | float | Lower-bound RMSD from best pose (Vina internal) |
| `rmsd_ub` | float | Upper-bound RMSD from best pose (Vina internal) |
| `control_rmsd` | float | For control run: RMSD vs crystallographic pose (Å); null for SCF-013 |
| `contacted_klifs_positions` | list[int] | KLIFS positions with ≥ 1 heavy atom within 4Å of pose |
| `contacted_selectivity_candidates` | list[int] | Subset of above that are selectivity candidates from spec-003 |

**Validation rules**:
- Control run must have `control_rmsd ≤ 2.0 Å` for SCF-013 results to be marked "validated setup"
- Only poses with `affinity_kcal_mol < 0` are written to output

---

### CrystalBrief
One document per target/scaffold pair — written as Markdown, not tabular.

| Field | Type | Description |
|-------|------|-------------|
| `target` | str | Gene symbol (e.g. "VRK1") |
| `scaffold_id` | str | Compound ID (e.g. "SCF-013") |
| `reference_pdb` | str | PDB ID used as reference crystal form (e.g. "6AC9") |
| `space_group` | str | Crystal space group from PDB header |
| `crystallisation_method` | str | E.g. "Vapor diffusion, hanging drop" |
| `crystallisation_conditions` | str | Full condition string from RCSB |
| `soaking_recommendations` | list[str] | Ordered list of recommended soaking conditions |
| `resolution_requirements` | dict[str, float] | Subpocket → minimum resolution (Å) |
| `flagged_scaffold_atoms` | list[str] | Atom descriptions likely to cause crystal packing clash |
| `selectivity_positions_to_resolve` | list[int] | KLIFS positions ranked by resolution requirement |

---

## File Outputs (per gene, in `data/results/{gene}/`)

| File | Content | Stage |
|------|---------|-------|
| `biomarker_report.md` | Ranked co-mutation table + top candidates | biomarker |
| `biomarker_results.csv` | Full BiomarkerResult table | biomarker |
| `binding_site_comparison.csv` | ThreeWayComparison (replaces 2-way when --include-vrk2) | structalign |
| `structural_selectivity_report.md` | Updated with VRK2 section (when --include-vrk2) | structalign |
| `docking_report.md` | Pose summary + selectivity contact map | dock |
| `docking_poses_{scaffold}.pdbqt` | All poses from Vina | dock |
| `docking_results_{scaffold}.csv` | DockingPose table | dock |
| `cocrystal_brief.md` | CrystalBrief document | cocrystal |

## Cache Files (in `data/cache/{gene}/`)

| File | Content |
|------|---------|
| `mutations_cache.json` | Parsed OmicsSomaticMutations for this gene's lineages |
| `receptor_{pdb_id}.pdbqt` | Prepared receptor for docking |
| `ligand_{scaffold_id}.pdbqt` | Prepared ligand for docking |
| `docking_control_rmsd.json` | Control run RMSD result (invalidated by --force) |
