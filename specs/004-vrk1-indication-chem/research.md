# Research: VRK1 Indication Validation and Chemistry Deepening

## US1 — CNS/Brain Biomarker Analysis

### Decision: CCLE mutation file
**Chosen**: `OmicsSomaticMutations.csv` from the DepMap 26Q1 manifest (same manifest URL as spec-003).
**Rationale**: Confirmed as the current file name. MAF-like format, joins to existing dependency cache on `ModelID`.
**Key columns**: `ModelID`, `Hugo_Symbol`, `Variant_Classification`, `isDeleterious`, `isCOSMIChotspot`.

### Decision: Enrichment test
**Chosen**: Fisher's exact test (`scipy.stats.fisher_exact`, `alternative='greater'`) + sample odds ratio.
**Rationale**: Standard for 2×2 contingency tables with small cell counts (individual mutations per lineage). No parametric assumptions. Odds ratio gives effect size independent of significance.
**Filtering**: Drop `Variant_Classification == 'Silent'`; keep `isDeleterious == True` to focus on likely-functional mutations. Group by `Hugo_Symbol` (gene level) for primary report; optionally by `(Hugo_Symbol, Variant_Classification)` for subtype detail.
**Gotchas**:
- File is ~100–200 MB — stream with `pd.read_csv(chunksize=...)` or read once and cache.
- Multiple rows per cell line per locus (multi-allelic) — deduplicate to one row per (ModelID, Hugo_Symbol) before building contingency table.
- Join key is `ModelID` — identical to the dependency cache format from spec-003.

---

## US2 — VRK2 Paralog Comparison

### Decision: VRK2 structure source
**Chosen**: KLIFS REST API first (crystal structure if available); AlphaFold AF-O95551-F1 fallback. Already wired in CLI from spec-003 — this story only extends the analytical output.
**pLDDT threshold**: 70 (inherited from spec-003 `_PLDT_CONFIDENCE_THRESHOLD`).

### Decision: Three-way selectivity classification
**Logic** (per KLIFS position):
- VRK1 == VRK2 AND VRK1 != EGFR → `"pan-VRK vs EGFR"` (avoid EGFR, hit both VRK paralogs)
- VRK1 != VRK2 AND VRK1 != EGFR → `"VRK1-specific"` (avoid both EGFR and VRK2)
- VRK1 != VRK2 AND VRK1 == EGFR → `"VRK2 vs VRK1+EGFR"` (VRK2 unique — lower priority)
- VRK1 == VRK2 == EGFR → `"conserved"` (no selectivity value)

---

## US3 — SCF-013 Docking

### Decision: Docking tool
**Chosen**: AutoDock Vina 1.2.7 via `pip install vina`. Free, local, reproducible, industry-recognised.
**Alternatives considered**: Web APIs (DiffDock, SwissTargetPrediction) — rejected for rate-limit fragility and non-reproducibility.

### Decision: Ligand preparation
**Chosen**: RDKit (SMILES → 3D conformer via `EmbedMolecule` + `UFFOptimizeMolecule`) → meeko (`MoleculePreparation`) → PDBQT.
**Receptor preparation**: meeko `mk_prepare_receptor.py` CLI or equivalent Python call for PDB → PDBQT.

### Decision: Docking box definition
**Chosen**: Centroid of binding-site residue Cα coordinates (from the 6Å NeighborSearch already computed in spec-003 structalign), plus 3Å margin each side.
**Box size**: `2 × (max_dist_from_centroid + 3.0)` Å per axis — approximately 22–28 Å for VRK1's ATP site.
**Exhaustiveness**: 32 (balanced accuracy/speed; ~1–2 min per ligand on a laptop).

### Decision: Positive-control validation
**Chosen**: Re-dock ANP (native ligand from 6AC9) into 6AC9. Success threshold: RMSD ≤ 2.0 Å vs crystallographic pose. If control fails, SCF-013 results are flagged as unvalidated.
**RMSD calculation**: BioPython `Superimposer` on matched heavy atoms (ligand-only, not Cα).

### Decision: Score parsing
**Chosen**: Parse `REMARK VINA RESULT:` lines from output PDBQT (standard Vina output format). Returns `(affinity_kcal_mol, rmsd_lb, rmsd_ub)` per pose.

### Decision: Contact mapping
**Chosen**: For each docked pose, count heavy atoms of the pose within 4Å of each KLIFS selectivity-candidate residue Cα (from spec-003 binding_site_comparison.csv). Annotate which of the 10 identified selectivity candidates are contacted.

### Additional dependencies
- `rdkit` (already in scientific Python environments; add to pyproject.toml if absent)
- `meeko>=0.5` for SMILES → PDBQT
- `vina>=1.2.7` for docking execution
- Both marked as optional extras in pyproject.toml; pipeline checks at startup and fails cleanly if missing.

---

## US4 — Co-Crystal Structure Brief

### Decision: Crystallisation data source
**Chosen**: RCSB PDB GraphQL API (`https://data.rcsb.org/graphql`) for programmatic retrieval of `exptl_crystal_grow.method` and `exptl_crystal_grow.pdbx_details`.

**Confirmed 6AC9 conditions** (live query):
- Method: Vapor diffusion, hanging drop
- Conditions: 27.5% w/v PEG 3350, 0.2 M ammonium sulfate, 0.1 M HEPES pH 7.0

### Decision: Homolog fallback
**Chosen**: If no VRK1-specific conditions found beyond 6AC9, search KLIFS for the next-closest VRK family structure by RMSD and retrieve its conditions as a secondary recommendation.

### Decision: Resolution requirements per subpocket
Based on kinase crystallography standards:
| Subpocket | Required resolution |
|-----------|-------------------|
| Gatekeeper (pos 45) | ≤ 2.5 Å (side-chain rotamer) |
| Hinge (pos 46–48) | ≤ 2.2 Å (backbone H-bonds with ligand) |
| DFG motif (pos 72–76) | ≤ 2.0 Å (DFG-in vs DFG-out conformation) |
| P-loop (pos 1–8) | ≤ 2.5 Å (flexible loop positioning) |

**Recommended target**: ≤ 2.2 Å overall to confidently resolve hinge and gatekeeper simultaneously. 6AC9 at 2.07 Å already meets this; a co-crystal at the same resolution or better is achievable.

### Decision: Scaffold atom flagging
**Chosen**: Identify rotatable bonds in SCF-013 (RDKit `rdMolDescriptors.CalcNumRotatableBonds`) and flag atoms that could adopt conformations clashing with crystal contacts (symmetry mates within 4 Å of binding site, using BioPython CRYST1 record).
