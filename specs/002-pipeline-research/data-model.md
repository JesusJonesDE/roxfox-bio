# Data Model: Drug Discovery Pipeline Research Tool

**Branch**: `002-pipeline-research` | **Date**: 2026-05-26

---

## Entities

### Target

Represents a protein drug target with its program mapping and disease context.

| Field | Type | Description |
|-------|------|-------------|
| gene_name | str | Human gene symbol (e.g., `VRK1`) |
| uniprot_id | str | UniProt accession (e.g., `Q99986`) |
| chembl_id | str | ChEMBL target ID (resolved at fetch time) |
| program_code | str | Internal program code (e.g., `RXF-001`) |
| indications | list[str] | Disease indications (e.g., `["SMA", "Oncology"]`) |
| display_name | str | Human-readable name for reports |

**Pre-configured targets**:
- VRK1 / Q99986 / RXF-001 / SMA + Oncology
- IGHMBP2 / P38935 / RXF-002 / SMARD1
- VCP / P55072 / RXF-003 / FTD

**Validation**: `gene_name` and `uniprot_id` required; `chembl_id` resolved during fetch stage.

---

### BioactivityRecord

Raw bioactivity measurement from ChEMBL — one record per compound-assay pair.

| Field | Type | Description |
|-------|------|-------------|
| molecule_chembl_id | str | ChEMBL compound identifier |
| canonical_smiles | str | SMILES structure string |
| standard_type | str | Assay type: `IC50`, `Ki`, or `Kd` |
| standard_value | float | Potency value in reported units |
| standard_units | str | Units (normalised to `nM`) |
| value_nm | float | Potency in nM (normalised) |
| pchembl_value | float \| None | −log10(potency in M), if available |
| target_chembl_id | str | ChEMBL target ID |
| assay_chembl_id | str | ChEMBL assay ID |

**Validation**: Records with `standard_value` null or `canonical_smiles` null are discarded at fetch time.

---

### Compound

A processed compound passing the potency filter, with computed properties.

| Field | Type | Description |
|-------|------|-------------|
| molecule_chembl_id | str | ChEMBL identifier (primary key) |
| smiles | str | Canonical SMILES |
| best_value_nm | float | Best (lowest) potency across all assay types |
| best_assay_type | str | Assay type of the best value |
| molecular_weight | float | Calculated MW (Da) |
| logp | float | Calculated Crippen logP |
| hbd | int | Hydrogen bond donor count |
| hba | int | Hydrogen bond acceptor count |
| rotatable_bonds | int | Rotatable bond count |
| ro5_violations | int | Number of Lipinski violations (0–4) |
| passes_ro5 | bool | True if violations < 2 |
| scaffold_id | str \| None | Murcko scaffold canonical SMILES (FK → Scaffold) |
| off_target_flags | int | Count of unrelated targets with activity ≤ 1µM |
| selectivity_flag | bool | True if off_target_flags > 3 |

**Derived from**: one or more `BioactivityRecord` entries (best value per compound selected).

---

### Scaffold

A Murcko framework grouping structurally related compounds.

| Field | Type | Description |
|-------|------|-------------|
| scaffold_smiles | str | Canonical SMILES of the Murcko framework (primary key) |
| scaffold_id | str | Hash-based short identifier (e.g., `SCF-001`) |
| compound_count | int | Number of compounds sharing this scaffold |
| median_potency_nm | float | Median `best_value_nm` across member compounds |
| best_potency_nm | float | Lowest `best_value_nm` across member compounds |
| target_gene | str | Gene name of the target this scaffold belongs to |

---

### Structure

An experimental or predicted protein structure entry.

| Field | Type | Description |
|-------|------|-------------|
| structure_id | str | PDB ID or `AF-{uniprot_id}-F1` for AlphaFold |
| source | str | `PDB` or `AlphaFold` |
| resolution_angstrom | float \| None | Resolution in Å (null for AlphaFold/NMR) |
| method | str | `X-ray`, `Cryo-EM`, `NMR`, or `Predicted` |
| has_ligand | bool | Whether a small-molecule ligand is co-crystallised |
| ligand_ids | list[str] | CCD codes of bound ligands |
| chain_ids | list[str] | Polymer chain identifiers |
| mean_plddt | float \| None | Mean AlphaFold pLDDT score (null for PDB) |
| deposition_date | str \| None | ISO date of PDB deposition |
| target_uniprot | str | UniProt ID of the target |

---

### CacheEntry

Tracks the state of a cached API response on disk.

| Field | Type | Description |
|-------|------|-------------|
| target_gene | str | Gene name the cache entry belongs to |
| stage | str | Pipeline stage: `fetch`, `analyze`, `report` |
| source | str | Data source name (e.g., `chembl`, `pdb`) |
| fetched_at | str | ISO 8601 timestamp of when data was fetched |
| record_count | int | Number of records in the cached response |
| file_path | str | Relative path to the JSON file on disk |
| is_valid | bool | Whether completeness check passed |

**Manifest file** (`pipeline_manifest.json`): maps `"{stage}:{target}"` → `CacheEntry` metadata. Used by the skip-or-execute logic.

---

### Dossier

Logical representation of the per-target markdown research document.

| Section | Source data |
|---------|-------------|
| Overview | Target config |
| Genetic Evidence | Open Targets scores |
| Bioactivity Summary | Compound count, potency distribution |
| Scaffold Highlights | Top 5 scaffolds by compound count |
| Structural Data | Structure inventory summary |
| Selectivity Profile | Off-target flags summary |
| Competitive Landscape | ChEMBL clinical candidates + ClinicalTrials studies |
| Data Gaps / Limitations | Any stages that returned no data |

---

## Directory Layout

```
data/
├── cache/
│   ├── VRK1/
│   │   ├── chembl_2026-05-26T143000.json
│   │   ├── open_targets_2026-05-26T143005.json
│   │   ├── pdb_2026-05-26T143010.json
│   │   ├── alphafold_2026-05-26T143012.json
│   │   └── clinical_trials_2026-05-26T143015.json
│   ├── IGHMBP2/
│   └── VCP/
├── results/
│   ├── VRK1/
│   │   ├── compounds_filtered.csv
│   │   ├── scaffolds.csv
│   │   ├── structures.csv
│   │   └── dossier.md
│   ├── IGHMBP2/
│   └── VCP/
└── pipeline_manifest.json
```

---

## State Transitions

```
Stage status per {stage}:{target} in manifest:

[absent] → fetch → [fetch:complete]
                      ↓
                  analyze → [analyze:complete]
                               ↓
                           report → [report:complete]

--force fetch    → resets fetch:complete, analyze:complete, report:complete
--force analyze  → resets analyze:complete, report:complete only
--force report   → resets report:complete only
```
