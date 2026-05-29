# Quickstart: IGHMBP2 Fragment Virtual Screening

## Prerequisites

```bash
# fpocket — already installed
fpocket --version

# sa-score for synthetic accessibility
pip install sa-score

# All other dependencies already installed: vina, meeko, rdkit, admet-ai, httpx, biopython
```

---

## Scenario 1: Full pipeline (overnight run)

```bash
# Start the full pipeline before going home
pipeline fragment --target IGHMBP2

# Next morning: check results
cat data/results/IGHMBP2/fragment_screen_report.md
pipeline validate --target IGHMBP2 --dashboard
```

Expected runtime: 8–14 hours. Expected output: `compounds_filtered.csv` with 10–50 candidates ready for `pipeline dock`.

---

## Scenario 2: Validate setup first (30-minute test run)

```bash
# Step 1: Find the binding pocket (< 1 min)
pipeline fragment --target IGHMBP2 --step pocket

# Step 2: Download fragment library (< 10 min)
pipeline fragment --target IGHMBP2 --step library

# Step 3: Test dock on 100 fragments only (~15 min)
pipeline fragment --target IGHMBP2 --step dock --library-size 100

# Check results
cat data/results/IGHMBP2/fragment_hits.csv | head -10
```

---

## Scenario 3: Continue from interrupted run

The pipeline is interruption-safe. If it stopped mid-docking:

```bash
# Simply re-run — completed fragments are cached and skipped
pipeline fragment --target IGHMBP2
```

Output will show `SKIP (cached)` for completed steps and resume from where it stopped.

---

## Scenario 4: After pipeline completes — run full validation pipeline

```bash
# The compounds_filtered.csv is in VRK1 schema — use the same pipeline
pipeline dock --target IGHMBP2 --all-scaffolds --top-n 10
pipeline validate --target IGHMBP2 --all-scaffolds --top-n 10 --gate admet
pipeline validate --target IGHMBP2 --dashboard
```

---

## Output files reference

| File | What it contains |
|---|---|
| `pocket_analysis.json` | Pocket centroid, volume, druggability score, box size |
| `fragments_ro3.smi` | 8,000 Ro3-filtered ZINC fragments (shared cache) |
| `fragment_hits.csv` | Top-50 docked fragments with affinity scores |
| `fragment_clusters.csv` | Cluster assignments, representative flags |
| `grown_candidates.csv` | Drug-like grown molecules, SA scores |
| `compounds_filtered.csv` | Final candidates in VRK1 schema for pipeline dock |
| `fragment_screen_report.md` | Full summary of all steps and findings |

---

## Interpreting results

**Fragment hits**: Affinity ≤ −5.0 kcal/mol is the threshold for meaningful fragment binding. Fragments at −4 to −5 may still be valid — use cluster diversity as the primary filter.

**Grown candidates**: SA score < 4 means synthetically accessible. SA score 4–6 is borderline. Above 6 = difficult to synthesise.

**ADMET**: BBB threshold relaxed to 0.3 for IGHMBP2 (neurological context). CYP thresholds unchanged.

**What to do with the output**: Take `compounds_filtered.csv` → run `pipeline dock` and `pipeline validate` → shortlist 3–5 candidates with best combined docking + ADMET profile → purchase from Enamine or Sigma (many ZINC fragments are commercially available) → run thermal shift assay to confirm binding.
