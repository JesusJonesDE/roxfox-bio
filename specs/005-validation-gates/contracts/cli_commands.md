# CLI Contracts: pipeline validate

## Command: `pipeline validate`

Run one or all validation gates for a scaffold × target pair.

### Signature
```
pipeline validate [OPTIONS]
```

### Options

| Option | Type | Required | Default | Description |
|---|---|---|---|---|
| `--target` / `-t` | str | yes* | — | Target gene (e.g. VRK1) |
| `--all` | flag | yes* | False | Run for all configured targets |
| `--scaffold` | str | no | — | Single scaffold ID (e.g. SCF-009) |
| `--all-scaffolds` | flag | no | False | Run all Ro5-passing, non-flagged scaffolds |
| `--top-n` | int | no | None | With --all-scaffolds: limit to top N by docking affinity |
| `--gate` | str | no | None | Run single gate: admet / mmgbsa / selectivity / md |
| `--dashboard` | flag | no | False | Print gate dashboard and exit (no new gates run) |
| `--force` | flag | no | False | Re-run gates even if cached results exist |
| `--md-metal` | flag | no | False | With --gate md: use Metal plugin for 50 ns on M1 Max GPU (requires source-built plugin) |
| `--md-cloud` | flag | no | False | With --gate md: dispatch 50 ns simulation to cloud GPU (requires RUNPOD_API_KEY env var) |
| `--data-dir` | path | no | data/ | Override data directory |

*Either --target or --all required.

### Behaviour

**Default (no --gate flag)**: Run gates in order ADMET → MM-GBSA → Selectivity → MD.
MD only runs if the first three gates all passed. If any gate fails, subsequent gates are skipped and their status is set to NOT_RUN.

**Single gate (--gate admet)**: Run only the specified gate regardless of other gate states.

**Dashboard (--dashboard)**: Read cached gate results, print terminal table, write `validation_dashboard.md`. Does not run any gates.

### Exit codes
- `0`: All gates ran (pass or fail are valid outcomes; exit 0)
- `1`: A gate encountered an ERROR (not FAIL — crash or parametrisation failure)

### Output files written

| File | When |
|---|---|
| `data/results/{gene}/validate_admet_{scaffold}.md` | After ADMET gate |
| `data/results/{gene}/validate_mmgbsa_{scaffold}.md` | After MM-GBSA gate |
| `data/results/{gene}/validate_selectivity_{scaffold}.md` | After selectivity gate |
| `data/results/{gene}/validate_md_{scaffold}.md` | After MD gate |
| `data/results/{gene}/validation_result_{scaffold}.json` | After any gate completes |
| `data/results/{gene}/validation_dashboard.md` | With --dashboard |
| `data/results/{gene}/validation_dashboard.json` | With --dashboard |
| `data/results/{gene}/wetlab_handoff_{scaffold}.md` | If scaffold passes all gates |

### Examples

```bash
# Run all gates for one scaffold
pipeline validate --target VRK1 --scaffold SCF-009

# Run only ADMET gate
pipeline validate --target VRK1 --scaffold SCF-009 --gate admet

# Run ADMET + MM-GBSA for top 5 scaffolds (skip MD)
pipeline validate --target VRK1 --all-scaffolds --top-n 5 --gate admet
pipeline validate --target VRK1 --all-scaffolds --top-n 5 --gate mmgbsa

# Show dashboard for all validated VRK1 scaffolds
pipeline validate --target VRK1 --dashboard

# Force re-run selectivity gate ignoring cache
pipeline validate --target VRK1 --scaffold SCF-009 --gate selectivity --force

# Run 50 ns MD via Metal plugin (local GPU, requires source-built plugin)
pipeline validate --target VRK1 --scaffold SCF-009 --gate md --md-metal

# Run 50 ns MD via cloud GPU (~$5-10, requires RUNPOD_API_KEY)
pipeline validate --target VRK1 --scaffold SCF-009 --gate md --md-cloud
```

---

## Gate: ADMET

### Inputs required
- Scaffold SMILES (from `compounds_filtered.csv`)

### Outputs
```
validate_admet_{scaffold_id}.md
```

### Pass/fail criteria
| Property | Field name | Pass threshold |
|---|---|---|
| BBB penetration | BBB_Martini | > 0.5 |
| CYP1A2 inhibition | CYP1A2_Inhibitor | < 0.3 |
| CYP2D6 inhibition | CYP2D6_Inhibitor | < 0.3 |
| CYP3A4 inhibition | CYP3A4_Inhibitor | < 0.3 |
| Solubility | Solubility | logS > −4 |
| Oral bioavailability | HIA_Hou | > 0.3 |

Gate PASS = all properties pass. Gate FAIL = any property fails (report lists which).

---

## Gate: MM-GBSA

### Inputs required
- Scaffold SMILES
- Cached receptor PDB (`data/cache/{gene}/structures/*.pdb`)
- Docking result PDBQT (`data/results/{gene}/docking_poses_{scaffold}.pdbqt`)

### Outputs
```
validate_mmgbsa_{scaffold_id}.md
```

### Pass/fail criteria
- ΔG (MM-GBSA) ≤ −7.0 kcal/mol → PASS
- ΔG > −7.0 kcal/mol → FAIL

---

## Gate: Selectivity

### Inputs required
- Scaffold SMILES
- Primary target docking ΔG (from MM-GBSA gate or Vina score if MM-GBSA not run)
- Off-target PDB structures (downloaded to `data/cache/shared/selectivity_panel/`)

### Outputs
```
validate_selectivity_{scaffold_id}.md
```

### Pass/fail criteria
- Selectivity Index = |target ΔG| / max(|off-target ΔG|)
- SI ≥ 10 → PASS
- SI < 10 → FAIL (report identifies the closest off-target)

---

## Gate: MD

### Inputs required
- Scaffold SMILES (for parametrisation)
- Cached receptor PDB
- Top docking pose PDBQT

### Modes
- **Tier 1 — Fast mode** (default): 2 ns implicit solvent (OBC/GB), ~1–2 h on M1 Max CPU. No extra dependencies.
- **Tier 2 — Metal plugin** (`--md-metal`): 50 ns explicit solvent (TIP3P), ~20–40 ns/day on M1 Max Metal GPU (1–3 days). Requires source-built `openmm-metal` plugin. Pipeline auto-detects if installed.
- **Tier 3 — Cloud** (`--md-cloud`): 50 ns explicit solvent, dispatched to RunPod A100. Requires `RUNPOD_API_KEY` env var. Completes in < 6 h at ~$5–10/run.

### Outputs
```
validate_md_{scaffold_id}.md
validate_md_{scaffold_id}_rmsd.csv   # time vs RMSD data
```

### Pass/fail criteria
- Fast mode: mean RMSD over final 1 ns ≤ 3.0 Å → PASS
- Full mode: mean RMSD over final 25 ns ≤ 3.0 Å → PASS

---

## Dashboard output format

### Terminal (Rich table)
```
VRK1 — Validation Dashboard
┌──────────┬────────┬─────────┬─────────────┬──────┬──────────┐
│ Scaffold │ ADMET  │ MM-GBSA │ Selectivity │  MD  │ Handoff? │
├──────────┼────────┼─────────┼─────────────┼──────┼──────────┤
│ SCF-009  │  PASS  │  PASS   │    PASS     │ PASS │    ✓     │
│ SCF-156  │  PASS  │  PASS   │    FAIL     │ —    │    ✗     │
│ SCF-130  │  FAIL  │  —      │     —       │  —   │    ✗     │
└──────────┴────────┴─────────┴─────────────┴──────┴──────────┘
```

### Markdown (validation_dashboard.md)
Same table in markdown format plus links to each gate report and handoff report.
