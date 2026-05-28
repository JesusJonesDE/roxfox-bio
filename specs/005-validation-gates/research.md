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

### Decision: OpenMM CPU platform with TIP3P + ff14SB + GAFF2 — **with revised runtime expectations**

**Critical finding**: OpenMM has **no Metal GPU backend** on Apple Silicon. The Metal plugin is not in the official release. This means all MD runs on CPU only.

**Revised runtime expectation on M1 Max**:
- CPU platform: 0.5–2 ns/day (all cores)
- 50 ns simulation: **25–100 days** — far beyond the 6-hour spec target
- **The spec's 6-hour MD target is not achievable locally for 50 ns explicit solvent**

**Revised approach**:

Option A — **Implicit solvent (OBC/GB) short run (10 ns)**: Run 10 ns in implicit solvent (no water box, ~25k atoms). Achieves ~5–15 ns/day on CPU. 10 ns completes in ~1–2 days. Still too slow for interactive use.

Option B — **Cloud GPU (recommended for production MD)**: AWS p3.2xlarge (NVIDIA V100) or Lambda Labs A10 achieves 50–100 ns/day for a kinase complex. Cost: ~$5 per 50 ns run. This is the practical path for the MD gate.

Option C — **Reduce MD gate scope to energy minimisation + short equilibration (1–2 ns)**: Sufficient to detect grossly unstable poses. Completes in ~1 hour on CPU.

### Decision for implementation:
- **Local fast mode** (default): 2 ns implicit solvent, RMSD check, completes ~1 hour on CPU. Gate pass/fail based on whether the pose remains in the pocket.
- **Full mode** (`--md-full` flag): Queue a 50 ns explicit solvent run; user is informed it will take ~2 days on CPU or can be dispatched to cloud. Results written when complete.
- **Update spec success criterion SC-004**: Revise from "6 hours for 50 ns" to "1 hour for 2 ns local fast mode; 50 ns available via --md-full flag".

**Install**:
```bash
conda install -c conda-forge openmm openmmforcefields openff-toolkit mdtraj pdbfixer
```

**Workflow**: PDB + top pose PDBQT → PDBFixer (cap termini, add missing atoms) → OpenMM Modeller (solvate/implicit) → minimise → 2 ns NVT → MDTraj RMSD analysis.

---

## Gate 5 — Dashboard

**Decision**: Rich terminal table + markdown file written to `data/results/{target}/validation_dashboard.md`. No additional dependencies — Rich is already used in the pipeline. Gate states stored as JSON in the cache system (reuses existing `CacheManager`).

---

## Dependency Installation Summary

```bash
# ADMET gate
pip install admet-ai

# MM-GBSA gate + MD gate
conda install -c conda-forge gmx_mmpbsa gromacs ambertools openmm openmmforcefields openff-toolkit mdtraj pdbfixer rdkit

# Already installed: vina, meeko, biopython, pandas, rich, httpx
```

**Total new packages**: ~8 conda + 1 pip. All ARM64-native or tested under Rosetta.

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
