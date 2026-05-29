# Data Model: IGHMBP2 Fragment-Based Virtual Screening

## Entities

### Pocket
A druggable cavity identified in the IGHMBP2 AlphaFold2 structure.

| Field | Type | Description |
|---|---|---|
| pocket_id | int | Rank from fpocket (1 = top) |
| gene_symbol | str | Target gene |
| score | float | fpocket composite score |
| druggability_score | float | fpocket druggability (0–1) |
| volume_A3 | float | Cavity volume in Å³ |
| centroid_x/y/z | float | Pocket centre coordinates |
| box_size_A | float | Docking box side length |
| plddt_mean | float | Mean AlphaFold pLDDT in pocket region |
| selected | bool | True = used for docking |

**Storage**: `data/results/IGHMBP2/pocket_analysis.json`

---

### Fragment
A small molecule from the ZINC fragment library.

| Field | Type | Description |
|---|---|---|
| fragment_id | str | e.g. "FRAG-00001" |
| smiles | str | Canonical SMILES (Ro3-filtered) |
| molecular_weight | float | MW in Da (≤ 250) |
| logp | float | cLogP (≤ 3) |
| hbd | int | H-bond donors (≤ 3) |
| hba | int | H-bond acceptors (≤ 3) |
| rotatable_bonds | int | Rotatable bonds (≤ 3) |
| zinc_id | str | ZINC database ID (if available) |
| scaffold | str | Murcko scaffold SMILES for deduplication |

**Storage**: `data/cache/shared/fragment_library/fragments_ro3.smi` (SMILES file)

---

### FragmentHit
A fragment with a docking result against the IGHMBP2 pocket.

| Field | Type | Description |
|---|---|---|
| fragment_id | str | Reference to Fragment |
| smiles | str | Fragment SMILES |
| affinity_kcal_mol | float | Vina top-pose affinity |
| n_poses | int | Number of valid poses |
| pose_file | Path | Path to docking_poses PDBQT |
| cluster_id | int | Assigned cluster (after clustering step) |
| is_representative | bool | True = cluster representative |
| rank | int | Global rank by affinity |

**Storage**: `data/results/IGHMBP2/fragment_hits.csv`

---

### Cluster
A group of chemically similar FragmentHits.

| Field | Type | Description |
|---|---|---|
| cluster_id | int | Sequential ID |
| representative_fragment_id | str | Best-scoring fragment in cluster |
| size | int | Number of fragments in cluster |
| centroid_smiles | str | SMILES of cluster representative |

**Storage**: `data/results/IGHMBP2/fragment_clusters.csv`

---

### GrownCandidate
A drug-like molecule derived from a fragment hit by structure-based growing.

| Field | Type | Description |
|---|---|---|
| candidate_id | str | e.g. "IGHMBP2-SCF-001" |
| parent_fragment_id | str | Source FragmentHit |
| smiles | str | Grown molecule SMILES |
| molecular_weight | float | MW 300–450 Da |
| logp | float | cLogP |
| hbd | int | H-bond donors |
| hba | int | H-bond acceptors |
| rotatable_bonds | int | Rotatable bonds |
| ro5_violations | int | Lipinski violations |
| passes_ro5 | bool | Ro5 compliant |
| sa_score | float | Synthetic accessibility (1=easy, 10=hard) |
| grow_method | str | "BRICS" or "SMARTS_rxn" |
| reaction_name | str | SMARTS transformation applied |

**Storage**: `data/results/IGHMBP2/grown_candidates.csv`

---

### FragmentScreenResult (pipeline state)
Tracks completion of each step for cache-based resumption.

| Field | Type | Description |
|---|---|---|
| gene_symbol | str | Target |
| step_pocket | bool | Pocket identification complete |
| step_library | bool | Fragment library ready |
| step_dock | bool | Fragment screen complete |
| step_cluster | bool | Clustering complete |
| step_grow | bool | Growing complete |
| step_admet | bool | ADMET screening complete |
| n_fragments_docked | int | Total fragments docked |
| n_hits | int | Top-N hits selected |
| n_clusters | int | Clusters found |
| n_grown | int | Grown candidates |
| n_candidates_final | int | Candidates in compounds_filtered.csv |
| completed_at | str | ISO8601 timestamp |

**Storage**: `data/cache/IGHMBP2/fragment_state.json`

---

## State Machine (step sequencing)

```
START
  ↓
pocket (US1) ← fpocket on AF2 structure
  ↓
library (US2) ← ZINC download + Ro3 filter
  ↓
dock (US3) ← batch Vina, exhaustiveness=4
  ↓
cluster (US4) ← Butina, ECFP4, cutoff=0.4
  ↓
grow (US5) ← BRICS + SMARTS reactions
  ↓
admet (US6) ← ADMET-AI + compounds_filtered.csv
  ↓
DONE (US7: full pipeline)
```

Each step caches output; `--step X` jumps to that step only; `--force` clears all cache.
