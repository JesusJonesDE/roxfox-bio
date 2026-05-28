# Research: Computational Validation Gates

**Branch**: `005-validation-gates` | **Date**: 2026-05-28

---

## Gate 1 — ADMET Prediction

### Decision: ADMET-AI (`pip install admet-ai`)

**Rationale**: Single pip install, no API key, covers all 5 required properties (BBB, metabolic stability, CYP1A2/2D6/3A4, solubility, oral bioavailability) in one call. Built on Chemprop v2 models trained on TDC benchmark datasets. No known Apple Silicon issues.

**Alternatives considered**:
- Chemprop v2 directly: same models but has a known MPS backend issue on Apple Silicon (GitHub issue #1351, March 2026). Workaround exists (`PYTORCH_ENABLE_MPS_FALLBACK=1`) but adds friction.
- SwissADME: web-only, requires HTTP call, no offline use.
- ADMETlab: web API only.
- DeepChem: heavy dependency tree, overkill for ADMET-only use.

**Usage**:
```python
from admet_ai import ADMETModel
model = ADMETModel()
results = model.predict(smiles="Nc1ncc(...)cc1...")
# returns dict: {"BBB_Martini": 0.82, "CYP2D6_Substrate": 0.12, ...}
```

**Expected runtime**: < 5 seconds per scaffold (CPU inference, small model).

**Pass thresholds** (from spec):
- BBB_Martini (penetration probability) > 0.5
- CYP1A2/2D6/3A4 inhibition probability < 0.3 each
- Solubility (LogS) > −4 (≈ 10 µg/mL for MW ~400)
- Oral bioavailability: HIA (Human Intestinal Absorption) > 0.3

---

## Gate 2 — MM-GBSA Rescoring

### Decision: gmx_MMPBSA on CPU via conda-forge

**Rationale**: Battle-tested tool with full workflow support (PDBQT → PDB conversion, GAFF2 parametrisation, MM-GBSA energy calculation). Conda packages available for osx-arm64. Runs in 2–5 minutes per pose on CPU (sufficient — this is rescoring, not MD).

**Alternatives considered**:
- OpenMM custom MM-GBSA: possible but requires bespoke implementation; Metal GPU would help here but official Metal support does not exist in OpenMM (see MD section below).
- AmberTools mm_pbsa: CPU-only, runs under Rosetta x86 emulation on M1 — acceptable but slower.
- MDAnalysis + parmed: analysis tools only, not energy solvers.

**Install**:
```bash
conda create -n mmpbsa -c conda-forge gmx_mmpbsa gromacs ambertools rdkit openmm openmmforcefields
```

**Forcefield**: GAFF2 (via AmberTools) for small molecule; ff19SB for protein. Standard, well-validated for kinase-ligand systems.

**Expected accuracy**: 1.5–2.5 kcal/mol RMSE vs. experimental ΔG for kinase-ligand sets; R² ≈ 0.6–0.75. Adequate for ranking, not for absolute predictions.

**Expected runtime**: 2–5 min per scaffold on M1 Max CPU. Well within the 2-hour spec requirement.

**Pass threshold**: ΔG ≤ −7.0 kcal/mol (from spec).

---

## Gate 3 — Selectivity Docking Panel

### Decision: Reuse existing AutoDock Vina pipeline with 4 off-target structures

**Rationale**: The pipeline already has a working dock stage. Running it against additional structures requires only fetching the right PDB files and reusing the same `run_dock()` function.

**Off-target panel** (from research):

| Kinase | Structure | Resolution | Ligand | Chain | Notes |
|--------|-----------|-----------|--------|-------|-------|
| VRK2 | AlphaFold2 AF-O95551-F1 | — | — | A | No crystal structure in PDB; use AF2 with pLDDT ≥ 70 warning |
| EGFR | 1M17 (already cached) | 2.60 Å | AQ4 | A | Erlotinib-bound; already used in structalign |
| CDK2 | 1E9H | 2.10 Å | ATP | A | ATP-bound; clean hinge definition |
| PLK1 | 2OKR | 2.80 Å | ADP | A | ATP-site bound; canonical kinase domain |

**Selectivity index**: `SI = target_affinity_kcal / best_offtarget_affinity_kcal` (both negative; higher magnitude = tighter binding). Pass if SI ≥ 10.

**Note on VRK2**: Since KLIFS has no VRK2 structure, we use the AlphaFold2 model. Docking into AF2 structures is less reliable than crystal structures — this is flagged in the report as a caveat, not a blocker.

**Expected runtime**: ~30 min per scaffold for 4 off-target dockings on M1 Max.

---

## Gate 4 — MD Pose Stability

### Decision: Three-tier MD strategy — fast local / Metal plugin / cloud GPU

**Thorough research findings (May 2026)**:

| Engine | Apple Silicon GPU | Notes |
|---|---|---|
| OpenMM (official) | CPU only | No Metal backend in 8.x official releases |
| OpenMM + Metal plugin | ~20–40 ns/day | Community plugin, unmaintained since Aug 2024, pinned to OpenMM 8.1, source build only |
| GROMACS | CPU only | conda-forge osx-arm64 is CPU-only; OpenCL path deprecated |
| NAMD, AMBER, JAX-MD | CPU only | CUDA-only GPU paths; no Apple Silicon GPU support |
| TorchMD-NET / MACE-OFF | CPU only on M1 | No validated OpenMM-Torch MPS path; experimental only |

**OpenMM Metal plugin details** (`philipturner/openmm-metal`):
- Uses Apple's `cl2Metal` OpenCL-to-Metal translation layer
- Benchmarks: ~20–40 ns/day for a 60k-atom kinase system on M1 Max
- Last commit: August 2024 — effectively unmaintained
- Pinned to OpenMM 8.1 (not compatible with 8.5+)
- Must build from source; no conda binary
- Labelled internally as "HIP" due to an OpenMM minimiser workaround
- Risk: medium — unmaintained, but source-buildable and functional as of last report

**Cloud GPU benchmark** (for reference):
- RunPod / Lambda Labs A100: 100–800 ns/day for 60k-atom system; ~$5–10 per 50 ns run

### Decision: Cloud GPU only (RunPod)

Local MD on Apple Silicon is not viable for production 50 ns simulations — the official OpenMM CPU path is too slow (0.5–2 ns/day) and the only Metal GPU path is an unmaintained community plugin pinned to a 2-year-old OpenMM version. Cloud is the cleanest, most reliable, and cheapest-per-run option.

**Cloud mode (RunPod A100, community cloud)**:
- User provides `RUNPOD_API_KEY` environment variable
- Pipeline prepares the system locally (PDBFixer + OpenMM parametrisation + HMR), serialises to OpenMM XML
- Submits job to RunPod community cloud A100 with a 90-minute hard timeout
- Before submission: pipeline estimates cost from instance type × estimated duration; refuses to submit if estimate exceeds `--md-max-cost` (default $5)
- ~$1–3 per 20 ns run on community cloud A100; completes in < 2 hours
- Pass/fail: mean ligand heavy-atom RMSD ≤ 3.0 Å over final 10 ns

**Four optimisations over naive approach**:
1. **20 ns instead of 50 ns** — sufficient for go/no-go pose stability; most pose drift occurs in first 5–10 ns
2. **Hydrogen Mass Repartitioning (HMR) → 4 fs timestep** — doubles throughput at no accuracy cost for stability assessment; well-supported in OpenMM
3. **10 Å solvation shell instead of 12 Å** — reduces atom count by ~20% (~50k vs ~60k atoms)
4. **RunPod community cloud** — ~$0.89–1.39/hr vs ~$2.49/hr for secure cloud; ~40–50% cheaper

**Combined effect**:
| Config | Simulation | Timestep | Atoms | Speed | Wall time | Cost |
|---|---|---|---|---|---|---|
| Naive | 50 ns | 2 fs | ~60k | ~200 ns/day | ~6 h | ~$12 |
| Optimised | 20 ns | 4 fs (HMR) | ~50k | ~500 ns/day | **~1 h** | **~$1–3** |

**Hard cost ceiling**:
- RunPod job `timeout=90min` — job is killed after 90 minutes regardless; pipeline analyses available trajectory
- Pipeline `--md-max-cost` flag (default $5) — estimates cost before submission; raises ERROR if estimate exceeds cap rather than submitting silently

**Local preparation** (always runs locally, ~5 min):
PDB + top pose PDBQT → PDBFixer (cap termini, add missing atoms) → OpenMM Modeller + GAFF2/ff14SB parametrisation → HMR application → energy minimisation → serialise system XML → submit to RunPod

**Install**:
```bash
# Local preparation tools
conda install -c conda-forge openmm openmmforcefields openff-toolkit pdbfixer

# Cloud submission
pip install runpod
```

**Updated spec success criterion SC-004**: "20 ns MD completes in under 2 hours via RunPod community cloud A100 at a cost of ~$1–3 per scaffold. Hard cost cap defaults to $5 per job."

---

## Gate 5 — Dashboard

**Decision**: Rich terminal table + markdown file written to `data/results/{target}/validation_dashboard.md`. No additional dependencies — Rich is already used in the pipeline. Gate states stored as JSON in the cache system (reuses existing `CacheManager`).

---

## Dependency Installation Summary

```bash
# ADMET gate
pip install admet-ai

# MM-GBSA gate
conda install -c conda-forge gmx_mmpbsa gromacs ambertools rdkit

# MD gate — local preparation + cloud submission
conda install -c conda-forge openmm openmmforcefields openff-toolkit pdbfixer
pip install runpod

# Already installed: vina, meeko, biopython, pandas, rich, httpx, mdtraj
```

**Total new packages**: ~7 conda + 2 pip. All ARM64-native or tested under Rosetta.
MD cloud submission requires `RUNPOD_API_KEY` environment variable.

---

## Architecture Decision: New `pipeline/stages/validate/` module

The validation pipeline is implemented as a new stage directory following existing patterns:

```
pipeline/stages/validate/
├── __init__.py
├── validate.py          # orchestrator: run_validate(), gate registry, dashboard
├── gates/
│   ├── admet.py         # Gate 1: ADMET-AI predictions
│   ├── mmgbsa.py        # Gate 2: gmx_MMPBSA rescoring
│   ├── selectivity.py   # Gate 3: off-target docking panel
│   └── md.py            # Gate 4: OpenMM pose stability
└── report.py            # WetLabHandoffReport generator
```

Each gate module exposes a single function `run_gate(gene_symbol, scaffold_id, settings, cache, console) -> GateResult` where `GateResult` is a dataclass with `status`, `score`, `reason`, `report_path`.

A new `pipeline validate` CLI command is added to `cli.py` following the existing pattern.
