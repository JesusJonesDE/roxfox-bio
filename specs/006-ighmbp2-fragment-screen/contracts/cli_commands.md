# CLI Contract: pipeline fragment

## Command: `pipeline fragment`

Run fragment-based virtual screening for a target with no published small molecule data.

### Signature
```
pipeline fragment [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--target` / `-t` | str | required | Target gene (e.g. IGHMBP2) |
| `--step` | str | None (all) | Run single step: pocket / library / dock / cluster / grow / admet |
| `--force` | flag | False | Re-run all steps ignoring cached results |
| `--top-n` | int | 50 | Number of top fragment hits to carry forward to clustering |
| `--exhaustiveness` | int | 4 | Vina exhaustiveness for fragment docking |
| `--library-size` | int | 8000 | Max fragments to download/use |
| `--data-dir` | path | data/ | Override data directory |

### Behaviour

**Default (no --step)**: Runs all 6 steps in sequence. Stops on first failure. Each completed step is cached and skipped on re-run.

**Single step (--step pocket)**: Runs only the specified step regardless of other step state. Useful for inspecting intermediate results.

**Cache logic**: Each step checks its own cache key before running. If cached and `--force` not set, prints SKIP and continues to next step.

### Exit codes
- `0`: All steps completed successfully
- `1`: A step failed (step name and error reported)

### Output files

| File | Step | Description |
|---|---|---|
| `data/results/{gene}/pocket_analysis.json` | pocket | Pocket centroid, volume, score |
| `data/cache/shared/fragment_library/fragments_ro3.smi` | library | Filtered fragment SMILES |
| `data/results/{gene}/fragment_hits.csv` | dock | All fragment docking results |
| `data/results/{gene}/fragment_clusters.csv` | cluster | Cluster assignments |
| `data/results/{gene}/grown_candidates.csv` | grow | Drug-like grown molecules |
| `data/results/{gene}/compounds_filtered.csv` | admet | Final candidates in VRK1 schema |
| `data/results/{gene}/fragment_screen_report.md` | admet | Full pipeline summary |

### Examples

```bash
# Run full pipeline (overnight)
pipeline fragment --target IGHMBP2

# Run pocket identification only
pipeline fragment --target IGHMBP2 --step pocket

# Download and prepare fragment library only
pipeline fragment --target IGHMBP2 --step library

# Run docking screen with smaller library for testing
pipeline fragment --target IGHMBP2 --step dock --library-size 500 --exhaustiveness 4

# Re-run everything from scratch
pipeline fragment --target IGHMBP2 --force

# Use larger top-N for clustering
pipeline fragment --target IGHMBP2 --step cluster --top-n 100
```

---

## Step: pocket

**Inputs**: AlphaFold2 PDB from `data/cache/IGHMBP2/`
**Outputs**: `pocket_analysis.json`

Runs fpocket on the IGHMBP2 structure, parses `_info.txt`, selects the top pocket by druggability score with volume > 200 Å³. Falls back to 200 Å³ threshold if nothing found above it.

Reports: pocket rank, score, volume, centroid (x,y,z), pLDDT mean.

---

## Step: library

**Inputs**: None (downloads from ZINC22)
**Outputs**: `data/cache/shared/fragment_library/fragments_ro3.smi`

Downloads ZINC fragment tranches, applies Ro3 filter (MW ≤ 250, HBD ≤ 3, HBA ≤ 3, logP ≤ 3, RotB ≤ 3), deduplicates by Murcko scaffold, randomly samples to `--library-size`. Falls back to bundled 500-fragment file if download fails.

---

## Step: dock

**Inputs**: `fragments_ro3.smi`, `pocket_analysis.json`, receptor PDB
**Outputs**: `fragment_hits.csv`

Docks each fragment using the existing Vina infrastructure (receptor PDBQT reused from dock stage cache). Pocket centroid from step 1 defines the box (20 Å cube). Exhaustiveness=4. Results cached per fragment — interruption-safe.

Top-N hits (default 50) written to `fragment_hits.csv`.

---

## Step: cluster

**Inputs**: `fragment_hits.csv`
**Outputs**: `fragment_clusters.csv`

Butina clustering at Tanimoto 0.6 threshold using Morgan ECFP4 fingerprints. Assigns cluster_id and is_representative flag to each hit.

---

## Step: grow

**Inputs**: `fragment_clusters.csv`, docking poses for cluster representatives
**Outputs**: `grown_candidates.csv`

Grows cluster representatives using BRICS combinatorial building and SMARTS reaction library. Filters: MW 300–450, Ro5, SA score < 4. Up to 3 grown candidates per representative retained.

---

## Step: admet

**Inputs**: `grown_candidates.csv`
**Outputs**: `compounds_filtered.csv`, `fragment_screen_report.md`

Runs ADMET-AI on all grown candidates. Applies BBB threshold of 0.3 (relaxed from 0.5 due to neurological disease context). Writes `compounds_filtered.csv` in exact VRK1 schema. Writes full pipeline report.
