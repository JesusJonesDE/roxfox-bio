# Quickstart: VRK1 Oncology Evidence — DepMap + Structural Alignment

**Branch**: `003-vrk1-oncology-evidence`

---

## Prerequisites

- conda environment `rxpipeline` with Python 3.11, RDKit, and existing pipeline deps
- Internet access (DepMap manifest API, RCSB PDB, KLIFS REST API — all public, no auth)
- VRK1 PDB structures already cached from a prior `pipeline fetch --target VRK1` run
- fpocket not required for these two commands

---

## New Dependencies

Install into the existing `rxpipeline` conda environment before first run:

```bash
conda activate rxpipeline
pip install biopython>=1.83 opencadd>=0.8
```

Verify:

```bash
python -c "import Bio; import opencadd; print('OK')"
```

---

## DepMap Cancer Dependency Analysis

### First Run (cold — downloads ~200 MB)

```bash
pipeline depmap --target VRK1
```

Expected runtime: 3–5 minutes (downloads CRISPRGeneEffect.csv and Model.csv).

Expected console output:
```
──────────── VRK1 — depmap ─────────────
  VRK1: 1,082 cell lines across 27 lineages
  VRK1: pan-essential: No
  VRK1: top lineage: Lung (median −0.62, 68% strongly dependent)
  VRK1: depmap_report.md written
```

Expected output files:
```
data/
├── cache/VRK1/
│   └── depmap_<TIMESTAMP>.json        # raw per-cell-line gene effect data
└── results/VRK1/
    ├── depmap_lineage_summary.csv     # ranked lineage table (≥3 lines)
    └── depmap_report.md               # markdown report with interpretation
```

### Second Run (cache hit — fast)

```bash
pipeline depmap --target VRK1
```

Expected runtime: < 2 seconds.

Expected console output:
```
  VRK1       depmap               SKIP  (cached 1082 records)
```

### Force Refresh

```bash
pipeline depmap --target VRK1 --force
```

Re-fetches from DepMap manifest even if cache is fresh.

---

## Structural Selectivity Alignment

### First Run (downloads EGFR reference 1M17)

```bash
pipeline structalign --target VRK1
```

Expected runtime: 2–4 minutes (one-time RCSB fetch for 1M17, KLIFS API calls).

Expected console output:
```
──────────── VRK1 — structalign ────────────
  VRK1: best structure 6AC9 (2.07 Å, ligand ANP, chain A)
  VRK1: EGFR reference 1M17 (2.60 Å, ligand AQ4, chain A)
  VRK1: binding site — VRK1: 38 residues | EGFR: 41 residues
  VRK1: KLIFS alignment RMSD: 1.23 Å (34 shared Cα)
  VRK1: gatekeeper — VRK1: Met131 | EGFR: Thr790
  VRK1: 18 positions differ | 6 selectivity candidates
  VRK1: structural_selectivity_report.md written
```

Expected output files:
```
data/
├── cache/VRK1/
│   └── structures/
│       └── 6AC9.pdb                       # already present from prior fetch
├── cache/shared/
│   └── structures/
│       └── 1M17.pdb                       # one-time RCSB download
└── results/VRK1/
    ├── binding_site_vrk1.csv              # 85-position KLIFS pocket for VRK1
    ├── binding_site_egfr.csv              # 85-position KLIFS pocket for EGFR
    ├── binding_site_comparison.csv        # residue-level diff table
    └── structural_selectivity_report.md  # report with selectivity hypothesis
```

### With VRK2 Three-Way Comparison

```bash
pipeline structalign --target VRK1 --include-vrk2
```

Adds VRK2 column to `binding_site_comparison.csv`. Falls back to AlphaFold AF-O95551-F1 if no PDB structure is found.

### Adjust Binding Site Cutoff

```bash
pipeline structalign --target VRK1 --cutoff 5.0
```

Tighter binding site definition (5 Å vs. default 6 Å). Fewer residues, higher confidence assignments.

---

## Running Both Analyses Together

The two commands are independent. Run sequentially or in parallel:

```bash
pipeline depmap --target VRK1 && pipeline structalign --target VRK1
```

Total expected runtime (cold): ~7–9 minutes. Both write to separate output files with no conflicts.

---

## Key Validation Checks After First Run

| Check | Command | Expected |
|-------|---------|----------|
| Lineage summary exists and non-empty | `wc -l data/results/VRK1/depmap_lineage_summary.csv` | ≥ 20 rows (header + ≥19 lineages) |
| Top lineage has negative median | `head -2 data/results/VRK1/depmap_lineage_summary.csv` | `median_effect` ≤ −0.4 |
| Pan-essential flag correct | `grep pan_essential data/results/VRK1/depmap_report.md` | Matches known VRK1 biology |
| Gatekeeper identified | `grep "Met131" data/results/VRK1/structural_selectivity_report.md` | Present |
| EGFR gatekeeper identified | `grep "Thr790" data/results/VRK1/structural_selectivity_report.md` | Present |
| Comparison table has rows | `wc -l data/results/VRK1/binding_site_comparison.csv` | > 30 rows |
| Cache hit is fast | Re-run either command | < 5s, prints SKIP |

---

## Troubleshooting

**`opencadd` KLIFS timeout**: KLIFS REST API is occasionally slow. `opencadd` defaults to 60s timeout; retry is automatic. If it fails repeatedly, run with `--force` after 5 minutes.

**`6AC9.pdb` not found**: Run `pipeline fetch --target VRK1` first to populate the PDB cache, then retry `structalign`.

**DepMap manifest returns 0 files**: The manifest URL (`https://depmap.org/portal/api/download/files`) returns a CSV of available files. If empty, DepMap may have changed their API. Check the URL manually and update `MANIFEST_URL` in the depmap stage.

**VRK1 column not found in CRISPRGeneEffect.csv**: DepMap encodes gene columns as `SYMBOL (EntrezID)`. VRK1's Entrez ID is 7443, so the column is `VRK1 (7443)`. If the column is absent, the release may have changed the format — check the header row.

**BioPython Superimposer unequal-length list error**: Ensure shared KLIFS positions are filtered to the intersection before calling `set_atoms`. Both fixed and moving atom lists must be equal length.
