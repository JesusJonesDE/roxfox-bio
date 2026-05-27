# Data Model: VRK1 Oncology Evidence

## Entities

### DepMapRecord
One row per cell line for a given gene.

| Field | Type | Source | Notes |
|-------|------|--------|-------|
| `model_id` | string | CRISPRGeneEffect.csv index | `ACH-XXXXXX` format |
| `gene_symbol` | string | column header | `VRK1` |
| `gene_effect` | float | CRISPRGeneEffect.csv | Chronos score; ≤ −0.5 = essential |
| `lineage` | string | Model.csv `OncotreeLineage` | e.g. `"Lung"`, `"Breast"` |
| `primary_disease` | string | Model.csv `OncotreePrimaryDisease` | e.g. `"Non-Small Cell Lung Cancer"` |
| `subtype` | string | Model.csv `OncotreeSubtype` | most specific OncotreeCode |
| `cell_line_name` | string | Model.csv `CCLEName` | human-readable |
| `dependency_tier` | enum | computed | `strongly_dependent` / `moderately_dependent` / `weakly_dependent` / `not_essential` |

**Validation rules:**
- `gene_effect` must be numeric; null values skipped (cell lines not screened for this gene)
- `lineage` nulls grouped as `"Unknown"` and excluded from lineage summaries requiring ≥ 3 lines
- `dependency_tier` derived from `gene_effect` thresholds: ≤ −0.5 → strongly_dependent; −0.5 to −0.3 → moderately_dependent; −0.3 to −0.1 → weakly_dependent; > −0.1 → not_essential

---

### LineageSummary
Aggregated per cancer lineage. One row per lineage.

| Field | Type | Derivation |
|-------|------|------------|
| `lineage` | string | groupby key |
| `n_lines` | int | count of cell lines with gene effect data |
| `median_effect` | float | median gene effect across lines |
| `mean_effect` | float | mean gene effect |
| `pct_strongly_dependent` | float | % lines with gene_effect ≤ −0.5 |
| `n_strongly_dependent` | int | count of lines with gene_effect ≤ −0.5 |
| `dependency_tier` | enum | derived from median_effect using same thresholds |
| `pan_essential_flag` | bool | True if pct_strongly_dependent > 70% across ALL lineages combined |

**Validation rules:**
- Only lineages with `n_lines ≥ 3` are included in the ranked output
- `pan_essential_flag` is a global flag computed once for the gene, not per lineage

---

### BindingSiteResidue
One row per KLIFS pocket position per kinase.

| Field | Type | Source | Notes |
|-------|------|--------|-------|
| `klifs_position` | int (1–85) | KLIFS API | Canonical pocket numbering |
| `subpocket` | string | KLIFS position map | e.g. `"Gatekeeper"`, `"Hinge"`, `"P-loop"`, `"DFG"` |
| `kinase` | string | — | `"VRK1"` or `"EGFR"` or `"VRK2"` |
| `pdb_id` | string | — | Source PDB structure |
| `chain` | string | — | Chain identifier |
| `pdb_residue_id` | int | KLIFS API | Residue number in PDB numbering |
| `amino_acid` | string (1-letter) | KLIFS API | Standard amino acid code; `"GAP"` if absent |
| `is_gatekeeper` | bool | klifs_position == 45 | |
| `is_hinge` | bool | klifs_position in 46–48 | |
| `distance_to_ligand_A` | float | computed | Min distance (Å) from residue heavy atoms to ligand; null if no ligand |

---

### BindingSiteComparison
One row per KLIFS position — the cross-kinase comparison table.

| Field | Type | Notes |
|-------|------|-------|
| `klifs_position` | int | |
| `subpocket` | string | |
| `vrk1_aa` | string | 1-letter code or `"GAP"` |
| `egfr_aa` | string | 1-letter code or `"GAP"` |
| `vrk2_aa` | string | 1-letter code, `"GAP"`, or `"N/A"` if VRK2 not run |
| `identical_vrk1_egfr` | bool | |
| `difference_type` | enum | `"identical"` / `"conservative"` / `"steric"` / `"charge"` / `"h_bond"` / `"gap"` |
| `selectivity_candidate` | bool | True if difference_type != identical AND subpocket in key positions |
| `is_gatekeeper` | bool | |
| `is_hinge` | bool | |
| `notes` | string | e.g. `"VRK1 Met131 vs EGFR Thr790 — steric difference exploitable"` |

**`difference_type` classification rules:**
- `identical`: same amino acid
- `conservative`: same physicochemical group (e.g. Val→Ile, Asp→Glu)
- `steric`: size difference > 2 carbon equivalents (e.g. Thr→Met, Gly→Phe)
- `charge`: charge difference (e.g. Lys→Glu)
- `h_bond`: H-bond donor/acceptor difference (e.g. Thr→Val)
- `gap`: one kinase has no residue at this position

---

## State Transitions

### DepMap analysis state
```
NOT_CACHED → FETCHING → CACHED → ANALYSED
                          ↑
                    (freshness check passes)
```

### Structural analysis state
```
NO_STRUCTURES → VRK1_READY → EGFR_FETCHED → KLIFS_MAPPED → ALIGNED → REPORT_WRITTEN
```

---

## Output Files

| File | Entity | Location |
|------|--------|----------|
| `depmap_vrk1_raw.csv` | DepMapRecord (all lines) | `data/cache/VRK1/` |
| `depmap_lineage_summary.csv` | LineageSummary | `data/results/VRK1/` |
| `depmap_report.md` | — | `data/results/VRK1/` |
| `binding_site_vrk1.csv` | BindingSiteResidue (VRK1) | `data/results/VRK1/` |
| `binding_site_egfr.csv` | BindingSiteResidue (EGFR) | `data/results/VRK1/` |
| `binding_site_comparison.csv` | BindingSiteComparison | `data/results/VRK1/` |
| `structural_selectivity_report.md` | — | `data/results/VRK1/` |
| `structures/6AC9.pdb` | — | `data/cache/VRK1/structures/` |
| `structures/1M17.pdb` | — | `data/cache/shared/structures/` |
