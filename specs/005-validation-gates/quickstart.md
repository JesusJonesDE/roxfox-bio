# Quickstart: Computational Validation Gates

## Prerequisites

```bash
# Install new dependencies
pip install admet-ai
conda install -c conda-forge gmx_mmpbsa gromacs ambertools openmm \
  openmmforcefields openff-toolkit mdtraj pdbfixer
```

Existing dependencies (already installed): `vina`, `meeko`, `biopython`, `rdkit`, `pandas`, `rich`, `httpx`.

---

## Scenario 1: Validate SCF-009 through all gates (recommended first run)

```bash
# Step 1: Confirm docking result exists
pipeline dock --target VRK1 --scaffold SCF-009

# Step 2: Run ADMET gate (fastest, ~5 seconds)
pipeline validate --target VRK1 --scaffold SCF-009 --gate admet

# Step 3: Run MM-GBSA rescoring (~5 minutes)
pipeline validate --target VRK1 --scaffold SCF-009 --gate mmgbsa

# Step 4: Run selectivity panel (~30 minutes)
pipeline validate --target VRK1 --scaffold SCF-009 --gate selectivity

# Step 5: Run MD fast mode (~1 hour)
pipeline validate --target VRK1 --scaffold SCF-009 --gate md

# Step 6: View dashboard
pipeline validate --target VRK1 --dashboard
```

Expected output from dashboard:
```
SCF-009 | ADMET: PASS | MM-GBSA: PASS | Selectivity: ? | MD: ?
```

---

## Scenario 2: Screen top 5 VRK1 scaffolds through ADMET + MM-GBSA

```bash
pipeline validate --target VRK1 --all-scaffolds --top-n 5 --gate admet
pipeline validate --target VRK1 --all-scaffolds --top-n 5 --gate mmgbsa
pipeline validate --target VRK1 --dashboard
```

This takes ~30 minutes and gives you a ranked view of which scaffolds survive the two cheapest gates before committing to selectivity docking.

---

## Scenario 3: Full 50 ns MD — Cloud GPU (~$5–10, < 6 hours)

```bash
# Set your RunPod API key (one-time):
export RUNPOD_API_KEY=your_key_here

# Run the MD gate:
pipeline validate --target VRK1 --scaffold SCF-009 --gate md
```

The pipeline prepares the system locally (PDBFixer + OpenMM parametrisation, ~5 min), submits to RunPod A100, polls for completion, and downloads results automatically. Expect 1–6 hours and ~$5–10 per scaffold.

---

## Scenario 4: Generate wet-lab handoff report

The handoff report is generated automatically when a scaffold passes all four gates. Check:

```bash
ls data/results/VRK1/wetlab_handoff_*.md
```

If the file exists, the scaffold is ready for CRO submission.

---

## Gate pass thresholds (quick reference)

| Gate | Pass condition |
|---|---|
| ADMET | BBB > 0.5, all CYP < 0.3, logS > −4, HIA > 0.3 |
| MM-GBSA | ΔG ≤ −7.0 kcal/mol |
| Selectivity | Target affinity ≥ 10× better than best off-target |
| MD (fast) | Mean RMSD ≤ 3.0 Å over final 1 ns |
| MD (full) | Mean RMSD ≤ 3.0 Å over final 25 ns |

---

## Troubleshooting

**"No docking result found"** → Run `pipeline dock --target VRK1 --scaffold SCF-009` first.

**MM-GBSA fails to converge** → Gate is marked ERROR (not FAIL). Re-run with `--force`. If it persists, the ligand may not be parametrisable with GAFF2.

**MD runs out of memory** → The system automatically falls back to implicit solvent (OBC/GB) with a warning. Explicit solvent requires ~8 GB for a kinase system; 64 GB is sufficient.

**ADMET-AI import error** → Run `pip install admet-ai` (not conda — use pip directly into the active environment).

**VRK2 selectivity warning** → VRK2 uses an AlphaFold2 model (no crystal structure exists). Docking into AF2 structures is less reliable; this is flagged in the report but does not block the gate.
