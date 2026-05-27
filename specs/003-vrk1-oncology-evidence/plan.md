# Implementation Plan: VRK1 Oncology Evidence — DepMap + Structural Alignment

**Branch**: `003-vrk1-oncology-evidence` | **Date**: 2026-05-27 | **Spec**: [spec.md](spec.md)

**Input**: Feature specification from `/specs/003-vrk1-oncology-evidence/spec.md`

## Summary

Add two new CLI commands to the rxpipeline: `pipeline depmap` (DepMap CRISPR Chronos cancer dependency analysis) and `pipeline structalign` (VRK1 vs. EGFR binding site structural alignment via KLIFS + BioPython). Both extend the existing pipeline pattern — new `stages/` modules, new `cli.py` commands, file-based cache, markdown report output.

---

## Technical Context

**Language/Version**: Python 3.11

**Primary Dependencies**:
- Existing: `typer`, `rich`, `pandas`, `httpx`, `tenacity`
- New: `biopython>=1.83` (PDB parsing, NeighborSearch, Superimposer), `opencadd>=0.8` (KLIFS REST API wrapper)

**Storage**: Files only — `data/cache/<GENE>/` for raw data, `data/results/<GENE>/` for analysis outputs. Reuses existing `CacheManager`.

**Testing**: pytest (existing `tests/` directory)

**Target Platform**: macOS / Linux laptop, headless, internet access for cold runs

**Project Type**: CLI extensions to existing pipeline

**Performance Goals**:
- DepMap cold run ≤ 5 min (manifest fetch + 200 MB CSV download + pandas parse)
- Structalign cold run ≤ 15 min (RCSB PDB fetch + KLIFS API + BioPython alignment)
- Cache hit ≤ 30 seconds for both commands

**Constraints**:
- DepMap CRISPRGeneEffect.csv is ~200 MB; use `pandas usecols` to read only the VRK1 gene column + index — do not load all 18,000 columns into memory
- KLIFS API is rate-limited; use opencadd's built-in retry; do not hammer with parallel requests
- BioPython `Superimposer.set_atoms` requires equal-length atom lists — filter to intersection of shared KLIFS positions before calling

**Scale/Scope**: 3 pipeline targets (VRK1, IGHMBP2, VCP); these two new commands operate on a single `--target` at a time

---

## Constitution Check

Constitution file is a template (not project-specific). No active gates apply. Both commands follow the existing pipeline pattern: single-responsibility module, file-based output, no new storage tier, no authentication required.

---

## Project Structure

### Documentation (this feature)

```text
specs/003-vrk1-oncology-evidence/
├── plan.md              ← this file
├── spec.md
├── research.md
├── data-model.md
├── quickstart.md
├── contracts/
│   └── cli_commands.md
└── tasks.md             ← created by /speckit.tasks
```

### Source Code (changes to repository root)

```text
pipeline/
├── stages/
│   ├── depmap/                          # NEW
│   │   ├── __init__.py
│   │   └── depmap.py                    # fetch manifest → download CSV → lineage analysis
│   └── structalign/                     # NEW
│       ├── __init__.py
│       └── structalign.py               # PDB parse → KLIFS mapping → Superimposer → report
├── cli.py                               # MODIFIED — add `depmap` and `structalign` commands
└── (all other existing files unchanged)

pyproject.toml                           # MODIFIED — add biopython>=1.83, opencadd>=0.8

tests/
├── unit/
│   ├── test_depmap.py                   # NEW — lineage aggregation, tier classification, cache
│   └── test_structalign.py             # NEW — residue extraction, difference classification
└── (existing tests unchanged)

data/
├── cache/
│   ├── VRK1/
│   │   ├── depmap_<TIMESTAMP>.json      # raw per-cell-line gene effect data
│   │   └── structures/6AC9.pdb          # already present from prior fetch
│   └── shared/
│       └── structures/1M17.pdb          # one-time RCSB download (EGFR reference)
└── results/
    └── VRK1/
        ├── depmap_lineage_summary.csv
        ├── depmap_report.md
        ├── binding_site_vrk1.csv
        ├── binding_site_egfr.csv
        ├── binding_site_comparison.csv
        └── structural_selectivity_report.md
```

**Structure Decision**: Extends the existing `pipeline/stages/<stage>/` pattern with two new stage directories. No new top-level project or package; the existing `rxpipeline` CLI is extended with two subcommands. Output paths follow the established `data/cache/<GENE>/` + `data/results/<GENE>/` convention.

---

## Implementation Phases

### Phase 1 — Setup

1. Add `biopython>=1.83` and `opencadd>=0.8` to `pyproject.toml` dependencies
2. Create `pipeline/stages/depmap/__init__.py` and `pipeline/stages/structalign/__init__.py` (empty)

### Phase 2 — DepMap Module (`pipeline/stages/depmap/depmap.py`)

Implements `run_depmap(target, settings, cache, force)`:

1. **Manifest fetch**: `GET https://depmap.org/portal/api/download/files` → parse CSV → find latest `CRISPRGeneEffect.csv` and `Model.csv` signed GCS URLs
2. **Cache check**: if `data/cache/<GENE>/depmap_*.json` exists and is fresh, return cached records
3. **Download**: stream `CRISPRGeneEffect.csv` with `pd.read_csv(url, usecols=[0, target_col])` where `target_col` matches `VRK1 (7443)`; stream `Model.csv` for `ModelID`, `OncotreeLineage`, `OncotreePrimaryDisease`, `OncotreeSubtype`, `CCLEName`
4. **Join + classify**: merge on `ModelID`; derive `dependency_tier` from Chronos thresholds; handle null lineages as `"Unknown"` (excluded from summaries requiring ≥3 lines)
5. **Lineage aggregation**: groupby lineage → `LineageSummary` records; compute `pan_essential_flag` globally (> 70% of all screened lines strongly dependent)
6. **Output**: write `depmap_<TIMESTAMP>.json` (cache), `depmap_lineage_summary.csv` (ranked), `depmap_report.md`
7. **Console**: print header + top lineage summary per contract format; print SKIP on cache hit

**Key implementation notes**:
- VRK1 column header is `VRK1 (7443)` — detect by `str.startswith(f"{gene_symbol} (")` in column index
- Signed GCS URLs are ephemeral — fetch manifest fresh before each download attempt
- Use `tenacity.retry` with exponential backoff (3 attempts) for manifest and download requests
- Exclude lineages with `n_lines < 3` from ranked output but still track them in the raw JSON cache

### Phase 3 — Structural Alignment Module (`pipeline/stages/structalign/structalign.py`)

Implements `run_structalign(target, settings, cache, force, include_vrk2, cutoff)`:

1. **Select VRK1 structure**: scan `data/cache/VRK1/structures/` for `.pdb` files; pick highest resolution ligand-bound structure (prefer 6AC9); parse resolution and ligand from HEADER/REMARK records
2. **Fetch EGFR reference**: if `data/cache/shared/structures/1M17.pdb` absent, download from `https://files.rcsb.org/download/1M17.pdb`; store in `data/cache/shared/structures/`
3. **KLIFS mapping**: use `opencadd.databases.klifs` to fetch the 85-position KLIFS pocket residue table for VRK1 (structure 6AC9) and EGFR (structure 1M17); map KLIFS position → PDB residue number + amino acid
4. **Binding site extraction**: `Bio.PDB.NeighborSearch` with `cutoff` (default 6.0 Å) from ligand heavy atoms; filter HETATM by `res.id[0].startswith("H_")`; exclude water (`"W"`)
5. **Superimposition**: build equal-length Cα atom lists for shared KLIFS positions; `Bio.PDB.Superimposer.set_atoms(fixed, moving)`; report RMSD
6. **Comparison table**: for each of 85 KLIFS positions, record VRK1 amino acid, EGFR amino acid; classify `difference_type` (identical / conservative / steric / charge / h_bond / gap); flag `selectivity_candidate` where type ≠ identical AND subpocket is adenine-binding region
7. **VRK2 (optional)**: if `--include-vrk2`, fetch VRK2 KLIFS structure or download AlphaFold AF-O95551-F1; add `vrk2_aa` column; exclude pLDDT < 70 positions with warning
8. **Output**: write `binding_site_vrk1.csv`, `binding_site_egfr.csv`, `binding_site_comparison.csv`, `structural_selectivity_report.md`
9. **Console**: print per contract format

**Key implementation notes**:
- Disordered atoms: use `atom.get_altloc() in ("", "A")` to select conformer A before distance calculations
- HETATM check: `res.id[0]` is `"H_ANP"` not `"H_"` — match with `startswith("H_")` and exclude `"W"`
- `Superimposer.set_atoms` requires lists of equal length — take `set(vrk1_klifs_pos) & set(egfr_klifs_pos)` before building atom lists
- VRK1 is atypical kinase — KLIFS coverage may be partial; fall back to sequence alignment vs. PKA (1ATP) if KLIFS returns < 30 positions

### Phase 4 — CLI Integration (`pipeline/cli.py`)

Add two new Typer commands following the existing `fetch` / `analyze` / `report` pattern:

```python
@app.command()
def depmap(target, all_targets, force, max_age, data_dir): ...

@app.command()
def structalign(target, all_targets, include_vrk2, force, cutoff, data_dir): ...
```

Flags match the contracts in `contracts/cli_commands.md` exactly. Both commands exit 0 on success, 1 on failure (with rich error message).

### Phase 5 — Tests

**`tests/unit/test_depmap.py`**:
- Tier classification edge cases (boundary values: −0.5 exactly, −0.3, −0.1, 0.0)
- Lineage aggregation with < 3 lines excluded
- Pan-essential flag: True when > 70%, False when ≤ 70%
- Column detection for `VRK1 (7443)` vs. `VRK1 (99999)` pattern

**`tests/unit/test_structalign.py`**:
- `difference_type` classification (identical, conservative Val→Ile, steric Thr→Met, charge Lys→Glu, h_bond Thr→Val, gap)
- Selectivity candidate flagging (adenine-binding subpocket positions)
- Equal-length atom list enforcement before Superimposer

---

## Complexity Tracking

*No constitution violations. Two new modules following established pipeline pattern. No new storage tier, no authentication, no additional services.*
