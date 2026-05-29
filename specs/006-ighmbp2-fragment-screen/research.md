# Research: IGHMBP2 Fragment-Based Virtual Screening

**Branch**: `006-ighmbp2-fragment-screen` | **Date**: 2026-05-29

---

## Pocket Identification — fpocket

**Decision**: subprocess + manual `_info.txt` parsing

**Rationale**: No official Python binding for fpocket. The `_info.txt` file is whitespace-delimited and contains all pocket descriptors needed. Parse with pandas `read_csv(sep=r'\s+')`.

**Key output files**:
- `<stem>_out/<stem>_info.txt` — per-pocket: Score, Druggability_Score, Volume, Mean_local_hydrophobic_density, Polarity_score
- `<stem>_out/pockets/pocket<N>_atm.pdb` — individual pocket atoms for centroid calculation

**Pocket selection**: Use Druggability_Score as primary sort; select top-ranked pocket with Volume > 200 Å³. If pLDDT is available (from AlphaFold JSON), filter out pockets whose centroid falls in regions with pLDDT < 70.

**Centroid calculation**: Average Cα coordinates of pocket-lining residues from `pocket<N>_atm.pdb`.

**Alternatives considered**: pypocketome (niche, poor maintenance), P2Rank (ML-based, heavier dependency), SiteMap (commercial). fpocket is the standard open-source choice and is already installed.

---

## Fragment Library

**Decision**: ZINC fragment-like subset, direct tranche download, ~10,000 compounds

**Rationale**: ZINC22 provides public fragment-like tranches without API key. For a focused pilot screen, 10,000 Ro3-filtered compounds (sampled from the fragment-like subset) provides good coverage of fragment chemical space while keeping docking time to ~6–8 hours on M1 Max at exhaustiveness=4.

**Download approach**:
- Primary: ZINC22 tranche files at `https://files.docking.org/2D/` (SMILES format, no auth)
- Filter: MW ≤ 250, HBD ≤ 3, HBA ≤ 3, logP ≤ 3, rotatable bonds ≤ 3 (Rule of Three)
- Sample: random 10,000 after filtering (or all if fewer)
- Fallback: bundled 500-fragment SMILES file (committed to repo) if download fails

**Deduplication**: Scaffold-based using Murcko scaffold (RDKit `MurckoScaffold.GetScaffoldForMol`); one compound per scaffold.

**Alternatives considered**: ChEMBL fragment subset (smaller, less diverse), commercial fragment libraries (require purchase), manually curated academic sets (too small).

---

## Fragment Docking

**Decision**: Reuse existing `run_dock()` infrastructure with exhaustiveness=4

**Rationale**: The dock stage already handles receptor PDBQT prep, ligand PDBQT prep (meeko), box definition, and Vina execution. Setting exhaustiveness=4 is standard for fragment screens (fragments have smaller search spaces due to fewer rotatable bonds). This gives ~4× speedup over the lead-docking exhaustiveness=32.

**Estimated throughput on M1 Max**: ~1–2 min/fragment → 10,000 fragments ≈ 10,000–20,000 min ≈ 7–14 hours. To stay within the 12-hour spec target, use 8,000 fragments maximum.

**Pocket box**: Use fpocket centroid + 20 Å box (ATP pockets in helicases are ~800–1200 Å³, requiring a larger box than kinase ATP sites).

**Resumption**: Cache each fragment result with key `fragment_dock_{fragment_id}` — completed fragments skip on re-run.

**Failure handling**: meeko parametrisation failures (common for fragments with unusual valences) are caught, logged, and skipped. Target: < 5% failure rate.

---

## Fragment Clustering

**Decision**: RDKit Butina clustering, Morgan ECFP4 fingerprints, cutoff=0.4 (60% similarity = same cluster)

**Rationale**: Butina is the standard RDKit clustering algorithm for this use case. Morgan radius=2 (ECFP4) is preferred for fragments because scaffold topology is more informative than pharmacophore features at fragment MW. Cutoff=0.4 (60% Tanimoto similarity threshold) gives meaningful chemotype separation at fragment size.

**Implementation**:
```python
from rdkit.ML.Cluster import Butina
from rdkit import DataStructs
from rdkit.Chem import AllChem
# Morgan r=2, 1024 bits (fragments have fewer heavy atoms)
fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols]
dists = []
for i in range(1, len(fps)):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend([1 - s for s in sims])
clusters = Butina.ClusterData(dists, len(fps), cutoff=0.4, isDistData=True)
```

**Representative selection**: Best-scoring fragment per cluster.

---

## Fragment Growing

**Decision**: RDKit BRICS + curated SMARTS reaction library

**Two-track approach**:

1. **BRICS combinatorial building** (`rdkit.Chem.BRICS`): Decompose top fragments at BRICS bonds, then reassemble with drug-like building blocks from a BRICS fragment pool (pre-built from ChEMBL drug fragments). Fast, combinatorial, chemically valid.

2. **SMARTS reaction growing**: Apply ~15 medicinal chemistry transformations to add drug-like substituents at exit vectors:
   - Amide formation, sulfonamide, urea, N-alkylation
   - Aromatic ring decoration (methyl, halogens, CF3)
   - Ring closure (4-6 membered rings)
   - Ether/thioether formation

**Filter after growing**: MW 300–450, Ro5-compliant, rotatable bonds ≤ 8.

**Synthetic accessibility**: Use `sascorer` from RDKit Contrib (`pip install sa-score` provides standalone package). SA score < 4 = synthetically accessible.

**Output**: ~50–200 grown candidates per top fragment representative; filter to top 3 per representative by predicted affinity estimate (Vina score after growing).

---

## ADMET Screening

**Decision**: Reuse existing `run_admet_gate()` from spec-005

**Rationale**: Already implemented, tested, and uses ADMET-AI with the correct field names. No changes needed.

**BBB threshold**: Apply relaxed threshold of 0.3 for IGHMBP2 candidates (same reasoning as VRK1: SMARD1 may involve neurological component; and ADMET-AI BBB predictions are conservative). Document this in the report.

---

## Final Output Schema

**Decision**: Exact match to VRK1 `compounds_filtered.csv` schema

Columns: `molecule_id, smiles, best_value_nm, best_assay_type, molecular_weight, logp, hbd, hba, rotatable_bonds, ro5_violations, passes_ro5, scaffold_id, source, off_target_flags, selectivity_flag`

For grown candidates:
- `best_value_nm` = None (no experimental data yet)
- `best_assay_type` = "fragment_screen_predicted"
- `scaffold_id` = `IGHMBP2-SCF-{N:03d}` (matching hitid convention)
- `source` = "fragment_virtual_screen"
- `off_target_flags` = 0 (not yet assessed)
- `selectivity_flag` = False

---

## Architecture: New `pipeline/stages/fragment/` module

```
pipeline/stages/fragment/
├── __init__.py
├── fragment.py      # orchestrator: run_fragment(), step registry, report
├── pocket.py        # Step 1: fpocket → pocket selection
├── library.py       # Step 2: ZINC download + Ro3 filter + deduplication
├── screen.py        # Step 3: batch Vina docking (reuses dock stage)
├── cluster.py       # Step 4: Butina clustering
├── grow.py          # Step 5: BRICS + SMARTS growing
└── output.py        # Step 6: ADMET + compounds_filtered.csv + report
```

New CLI command: `pipeline fragment --target IGHMBP2 [--step pocket|library|dock|cluster|grow|admet] [--force] [--top-n 50]`
