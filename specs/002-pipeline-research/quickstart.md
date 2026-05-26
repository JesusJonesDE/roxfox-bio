# Quickstart: Drug Discovery Pipeline Research Tool

**Branch**: `002-pipeline-research`

---

## Prerequisites

- conda (Miniconda or Anaconda) — required for RDKit
- Python 3.11+
- Internet access (public APIs, no auth required)

---

## Setup

```bash
# 1. Create conda environment
conda create -n rxpipeline python=3.11 -y
conda activate rxpipeline

# 2. Install RDKit (must be via conda)
conda install -c conda-forge rdkit -y

# 3. Install remaining Python dependencies
pip install -e .

# 4. Verify installation
pipeline --help
```

Expected output:
```
Usage: pipeline [OPTIONS] COMMAND [ARGS]...

  RoxFox Bio drug discovery research pipeline.

Commands:
  fetch    Fetch raw data from all external sources.
  analyze  Run analysis on cached data.
  report   Generate markdown dossier from analysis results.
  run      Execute all stages in sequence.
  status   Show cache and manifest state.
```

---

## First Run — VRK1 (RXF-001)

```bash
# Full pipeline, single target
pipeline run --target VRK1
```

Expected runtime: ~3–5 minutes (first run, all fetches live)

Expected output directory:
```
data/
├── cache/VRK1/
│   ├── chembl_2026-05-26T143000.json
│   ├── open_targets_2026-05-26T143005.json
│   ├── pdb_2026-05-26T143010.json
│   ├── alphafold_2026-05-26T143012.json
│   └── clinical_trials_2026-05-26T143015.json
└── results/VRK1/
    ├── compounds_filtered.csv
    ├── scaffolds.csv
    ├── structures.csv
    └── dossier.md
```

---

## Second Run — Cached (Fast)

```bash
pipeline run --target VRK1
```

Expected runtime: < 5 seconds. All stages report `SKIP (cached)`.

---

## All Three Targets

```bash
pipeline run --all
```

Expected runtime: ~10–15 minutes (first run across all targets)

---

## Force Refresh

```bash
# Re-fetch all data for VRK1 (data older than default 30 days)
pipeline run --target VRK1 --force fetch

# Re-run only the report stage
pipeline report --target VRK1 --force report

# Full refresh of everything
pipeline run --all --force all
```

---

## Check Status

```bash
pipeline status
```

Shows per-target freshness and completion state.

---

## Key Scenarios to Validate After First Run

| Scenario | Command | Expected result |
|----------|---------|-----------------|
| VRK1 compounds CSV non-empty | `wc -l data/results/VRK1/compounds_filtered.csv` | > 1 row (header + data) |
| VRK1 scaffold CSV has clusters | `head -5 data/results/VRK1/scaffolds.csv` | Multiple scaffolds with counts |
| VRK1 structures CSV lists PDB entries | `head -5 data/results/VRK1/structures.csv` | Rows with PDB IDs and resolution |
| VRK1 dossier has all sections | `grep "^##" data/results/VRK1/dossier.md` | 8 section headers |
| Re-run uses cache | `pipeline run --target VRK1` (2nd time) | All stages SKIP, < 5s total |
| Force refresh works | `pipeline fetch --target VRK1 --force fetch` | New timestamp in cache filenames |
| IGHMBP2 runs without compounds | `pipeline run --target IGHMBP2` | Completes; dossier notes data gap if no ChEMBL data |

---

## Troubleshooting

**RDKit import error**: Ensure you're in the `rxpipeline` conda environment (`conda activate rxpipeline`). RDKit installed via pip will not work reliably.

**Network timeout on ChEMBL**: ChEMBL has no hard rate limit but large targets (VRK1 with ~1200 records) take 30–60s to paginate fully. This is normal.

**`pipeline: command not found`**: Run `pip install -e .` from the project root with the conda env active.

**Stale manifest after corrupted cache**: Run `pipeline fetch --target {TARGET} --force fetch` to regenerate the cache for that target.
