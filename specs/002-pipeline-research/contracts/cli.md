# CLI Contract: Drug Discovery Pipeline Research Tool

**Branch**: `002-pipeline-research` | **Date**: 2026-05-26

The tool is installed as a `pipeline` shell command. All subcommands share the `--target` and `--force` flags.

---

## Global Flags

| Flag | Values | Default | Description |
|------|--------|---------|-------------|
| `--target` | `VRK1`, `IGHMBP2`, `VCP`, or any gene name | required (unless `--all`) | Restrict to one target |
| `--all` | flag | false | Run across all configured targets |
| `--force` | `fetch`, `analyze`, `report`, `all` | none | Bypass cache for specified stage(s) |
| `--max-age` | integer (days) | 30 | Override cache freshness threshold |
| `--data-dir` | path | `./data` | Override default data directory |

---

## Subcommands

### `pipeline fetch`

Fetch raw data from all external sources for the specified target(s). Writes JSON files to `data/cache/{target}/`. Updates manifest.

```
pipeline fetch --target VRK1
pipeline fetch --all
pipeline fetch --target VRK1 --force fetch
```

**stdout on cache hit**:
```
[VRK1] fetch  chembl          SKIP  (cached 3 days ago, 1243 records)
[VRK1] fetch  open_targets    SKIP  (cached 3 days ago)
[VRK1] fetch  pdb             SKIP  (cached 3 days ago, 26 structures)
[VRK1] fetch  alphafold       SKIP  (cached 3 days ago)
[VRK1] fetch  clinical_trials SKIP  (cached 3 days ago, 0 trials)
```

**stdout on fresh fetch**:
```
[VRK1] fetch  chembl          OK    (1243 records, 4.2s)
[VRK1] fetch  open_targets    OK    (score: 0.84, 2.1s)
[VRK1] fetch  pdb             OK    (26 structures, 1.8s)
[VRK1] fetch  alphafold       OK    (model AF-Q99986-F1, 1.1s)
[VRK1] fetch  clinical_trials OK    (0 trials, 0.9s)
```

**Exit codes**: 0 = success, 1 = partial failure (some sources failed), 2 = complete failure

---

### `pipeline analyze`

Run analysis stages on cached data: bioactivity filtering, Lipinski, scaffold clustering, structural inventory, selectivity profiling. Writes CSVs to `data/results/{target}/`.

```
pipeline analyze --target VRK1
pipeline analyze --all --force analyze
```

**stdout**:
```
[VRK1] analyze  bioactivity   OK    (847 compounds → 312 pass 10µM filter → 289 pass Ro5)
[VRK1] analyze  scaffolds     OK    (289 compounds → 94 unique scaffolds)
[VRK1] analyze  structures    OK    (26 PDB + 1 AlphaFold model)
[VRK1] analyze  selectivity   OK    (12 compounds flagged for off-target activity)
```

**Requires**: `fetch` stage complete for the target (errors if not).

---

### `pipeline report`

Generate the markdown dossier from analysis results. Writes `data/results/{target}/dossier.md`.

```
pipeline report --target VRK1
pipeline report --all --force report
```

**stdout**:
```
[VRK1] report  dossier        OK    → data/results/VRK1/dossier.md
```

**Requires**: `analyze` stage complete for the target (errors if not).

---

### `pipeline run`

Execute all three stages in sequence: fetch → analyze → report. Equivalent to calling each subcommand in order with the same flags.

```
pipeline run --target VRK1
pipeline run --all
pipeline run --target IGHMBP2 --force all
pipeline run --all --force fetch
```

**stdout**: Combined output of all three stages, separated by target header.

---

### `pipeline status`

Show the current cache and manifest state for all configured targets.

```
pipeline status
```

**stdout**:
```
Target    Fetch              Analyze            Report
VRK1      2026-05-23 14:30  2026-05-23 14:35  2026-05-23 14:38  ✓ current
IGHMBP2   2026-04-10 09:00  2026-04-10 09:05  —                 ⚠ stale (46d)
VCP       —                 —                 —                 ✗ not run
```

---

## Output Files Contract

### `compounds_filtered.csv`

| Column | Type | Description |
|--------|------|-------------|
| molecule_chembl_id | str | ChEMBL compound identifier |
| smiles | str | Canonical SMILES |
| best_value_nm | float | Best potency in nM |
| best_assay_type | str | IC50 / Ki / Kd |
| molecular_weight | float | MW in Da |
| logp | float | Crippen logP |
| hbd | int | H-bond donors |
| hba | int | H-bond acceptors |
| rotatable_bonds | int | Rotatable bonds |
| ro5_violations | int | 0–4 |
| passes_ro5 | bool | True if violations < 2 |
| scaffold_id | str | Murcko scaffold identifier |
| selectivity_flag | bool | True if > 3 off-target hits ≤ 1µM |

### `scaffolds.csv`

| Column | Type | Description |
|--------|------|-------------|
| scaffold_id | str | Short identifier |
| scaffold_smiles | str | Murcko scaffold SMILES |
| compound_count | int | Number of member compounds |
| median_potency_nm | float | Median potency of cluster |
| best_potency_nm | float | Best (lowest) potency in cluster |

### `structures.csv`

| Column | Type | Description |
|--------|------|-------------|
| structure_id | str | PDB ID or AlphaFold accession |
| source | str | PDB or AlphaFold |
| resolution_angstrom | float \| — | Resolution in Å |
| method | str | X-ray / Cryo-EM / NMR / Predicted |
| has_ligand | bool | Co-crystallised small molecule |
| ligand_ids | str | Comma-separated CCD codes |
| mean_plddt | float \| — | AlphaFold confidence score |
| deposition_date | str \| — | ISO date |

### `dossier.md` sections (required, in order)

1. `# {Target} — Research Dossier`
2. `## Overview` — target description, program code, indications
3. `## Genetic Evidence` — Open Targets scores, top disease associations
4. `## Bioactivity Summary` — compound counts at each filter stage, potency distribution
5. `## Scaffold Highlights` — top 5 scaffolds by size, with potency summary
6. `## Structural Data` — structure count, best structure recommendation
7. `## Selectivity Profile` — off-target liability summary
8. `## Competitive Landscape` — approved drugs, clinical candidates, active trials
9. `## Data Gaps & Limitations` — any missing or absent data with explanation
