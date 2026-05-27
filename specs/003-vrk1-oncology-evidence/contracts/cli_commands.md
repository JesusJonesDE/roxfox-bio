# CLI Command Contracts

## `pipeline depmap`

Fetch DepMap CRISPR gene effect data for the target gene and produce a ranked cancer lineage report.

```
pipeline depmap [OPTIONS]

Options:
  -t, --target TEXT     Target gene name (e.g. VRK1)  [required or --all]
  --all                 Run for all configured targets
  --force               Re-fetch even if cache is fresh
  --max-age INTEGER     Cache freshness threshold in days  [default: 30]
  --data-dir PATH       Override data directory
  --help
```

**Exit codes:**
- `0` — analysis completed successfully, all outputs written
- `1` — analysis failed (fetch error, parse error, or target not in DepMap)

**Outputs written on success:**
- `data/cache/<GENE>/depmap_<TIMESTAMP>.json` — raw per-cell-line gene effect data
- `data/results/<GENE>/depmap_lineage_summary.csv` — ranked lineage table
- `data/results/<GENE>/depmap_report.md` — markdown report

**Console output format (success):**
```
──────────── VRK1 — depmap ─────────────
  VRK1: 1,082 cell lines across 27 lineages
  VRK1: pan-essential: No
  VRK1: top lineage: Lung (median −0.62, 68% strongly dependent)
  VRK1: depmap_report.md written
```

**Console output format (cache hit):**
```
  VRK1       depmap               SKIP  (cached 1082 records)
```

---

## `pipeline structalign`

Download reference structures, extract ATP binding sites, align VRK1 vs. EGFR (and optionally VRK2), produce a residue-level selectivity comparison.

```
pipeline structalign [OPTIONS]

Options:
  -t, --target TEXT     Target gene name (e.g. VRK1)  [required or --all]
  --all                 Run for all configured targets
  --include-vrk2        Include VRK2 in three-way comparison  [default: True]
  --force               Re-run even if output files exist
  --cutoff FLOAT        Binding site distance cutoff in Å  [default: 6.0]
  --data-dir PATH       Override data directory
  --help
```

**Exit codes:**
- `0` — alignment completed, all outputs written
- `1` — structural analysis failed (no suitable PDB structure, KLIFS unavailable, or BioPython error)

**Outputs written on success:**
- `data/cache/<GENE>/structures/<PDB_ID>.pdb` — downloaded PDB files
- `data/results/<GENE>/binding_site_<GENE>.csv` — binding site residues for target
- `data/results/<GENE>/binding_site_egfr.csv` — binding site residues for EGFR reference
- `data/results/<GENE>/binding_site_comparison.csv` — residue-level comparison table
- `data/results/<GENE>/structural_selectivity_report.md` — markdown report

**Console output format (success):**
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
