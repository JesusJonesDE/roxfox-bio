# Data Model: Computational Validation Gates

## Entities

### GateStatus (enum)
```
PASS | FAIL | ERROR | PENDING | NOT_RUN
```

### GateResult
Represents the outcome of a single validation gate for one scaffold × target pair.

| Field | Type | Description |
|---|---|---|
| gate_name | str | "admet" / "mmgbsa" / "selectivity" / "md" |
| status | GateStatus | PASS / FAIL / ERROR / PENDING / NOT_RUN |
| score | float | Primary numeric result (probability, kcal/mol, SI ratio, RMSD Å) |
| reason | str | Human-readable explanation of pass/fail decision |
| details | dict | Gate-specific sub-scores (e.g. BBB=0.82, CYP2D6=0.11 for ADMET) |
| report_path | Path | Path to the gate's markdown report file |
| duration_s | float | Wall-clock seconds the gate took to run |
| timestamp | str | ISO8601 when the gate ran |

**Cache key**: `validate_{gate_name}_{scaffold_id}` (per target in CacheManager)

---

### ValidationResult
All four gate results for one scaffold × target pair.

| Field | Type | Description |
|---|---|---|
| gene_symbol | str | e.g. "VRK1" |
| scaffold_id | str | e.g. "SCF-009" |
| smiles | str | Canonical SMILES used |
| gates | dict[str, GateResult] | Keyed by gate name |
| overall_pass | bool | True only if all run gates passed |
| handoff_ready | bool | True if all 4 gates passed |
| created_at | str | ISO8601 |

**Storage**: `data/results/{gene}/validation_result_{scaffold_id}.json`

---

### SelectivityPanel
Configurable set of off-target structures for one target class.

| Field | Type | Description |
|---|---|---|
| target_class | str | "kinase" |
| off_targets | list[OffTargetEntry] | Ordered list of off-target definitions |

### OffTargetEntry

| Field | Type | Description |
|---|---|---|
| gene | str | e.g. "EGFR" |
| pdb_id | str | e.g. "1M17" (or "AF2" for AlphaFold) |
| source | str | "pdb" or "alphafold" |
| ligand_id | str | e.g. "AQ4" — used to define docking box centroid |
| chain_id | str | e.g. "A" |
| warning | str | Optional; shown when using AF2 or low-res structure |

**Default kinase panel** (hardcoded, overridable):
- VRK2 → AF2 (O95551), chain A, warning: "AlphaFold2 model — lower confidence"
- EGFR → 1M17, chain A, ligand AQ4
- CDK2 → 1E9H, chain A, ligand ATP
- PLK1 → 2OKR, chain A, ligand ADP

---

### GateDashboard
Cross-scaffold view for a target.

| Field | Type | Description |
|---|---|---|
| gene_symbol | str | Target |
| scaffolds | list[str] | All scaffold IDs with any gate result |
| results | dict[str, ValidationResult] | Keyed by scaffold_id |
| generated_at | str | ISO8601 |

**Storage**: `data/results/{gene}/validation_dashboard.md` (markdown) + `validation_dashboard.json` (machine-readable)

---

### WetLabHandoffReport
Generated only for scaffolds where `handoff_ready == True`.

| Field | Type | Description |
|---|---|---|
| gene_symbol | str | Target |
| scaffold_id | str | Scaffold |
| smiles | str | Canonical SMILES |
| admet_summary | dict | Key ADMET scores |
| mmgbsa_dg | float | MM-GBSA ΔG in kcal/mol |
| selectivity_index | float | SI ratio |
| md_rmsd_mean | float | Mean RMSD Å over final simulation window |
| docking_report_path | Path | Link to existing docking report |
| gate_reports | dict[str, Path] | Links to each gate report |

**Storage**: `data/results/{gene}/wetlab_handoff_{scaffold_id}.md`

---

## State Transitions

```
NOT_RUN → PENDING (gate starts)
PENDING → PASS    (gate completes, threshold met)
PENDING → FAIL    (gate completes, threshold not met)
PENDING → ERROR   (gate crashed or could not converge)
Any     → NOT_RUN (cache cleared / --force)
```

Gate ordering enforced: ADMET → MM-GBSA → Selectivity → MD (MD only if previous 3 = PASS)
