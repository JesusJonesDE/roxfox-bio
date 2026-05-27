# Research: VRK1 Oncology Evidence

## DepMap Access

**Decision:** Use the DepMap manifest REST API to fetch signed GCS download URLs, then stream `CRISPRGeneEffect.csv` and `Model.csv` from those URLs. Parse using pandas with `usecols` to keep memory usage small (target gene column only).

**Rationale:** There is no single-gene query API on the Broad DepMap portal for CRISPR Chronos data. The manifest endpoint (`https://depmap.org/portal/api/download/files`) returns a CSV of all available files with pre-signed GCS download URLs per release. This is the officially documented access method. Signed URLs are ephemeral — they must be fetched fresh before each download.

**Alternatives considered:**
- `depmap-downloader` (PyPI): wraps the same manifest API but still downloads the full file; adds a dependency for no benefit in this use case.
- Sanger DepMap API (`depmap.sanger.ac.uk`): has a proper REST gene-slice API but is a separate dataset (CRISPR Avana, not Chronos); different cell line set and effect scoring.
- Bulk static GCS URLs: return 403 without a signed token — not viable.

**Key data facts:**
- Latest release: DepMap Public 26Q1 (2026-04-01)
- CRISPRGeneEffect.csv: ~1,100 rows (cell lines, indexed by `ModelID` = `ACH-XXXXXX`), ~18,000 gene columns
- VRK1 column header: `VRK1 (7443)` (gene symbol + Entrez ID)
- Lineage join via `Model.csv` on `ModelID` column
- Lineage columns: `OncotreeLineage` (broadest, e.g. `"Lung"`), `OncotreePrimaryDisease`, `OncotreeSubtype`
- No authentication required for manifest or signed URLs
- Memory: using `usecols` to select only VRK1 column reduces parse footprint to ~a few MB

**Dependency tier thresholds (established DepMap conventions):**
- Strongly dependent: gene effect ≤ −0.5
- Moderately dependent: −0.5 to −0.3
- Weakly dependent: −0.3 to −0.1
- Not essential: > −0.1
- Pan-essential flag: VRK1 strongly dependent in > 70% of all screened lines

---

## Structural Alignment

**Decision:** BioPython for PDB parsing and binding site extraction + KLIFS REST API for canonical kinase pocket residue mapping + BioPython Superimposer for structural alignment. Use `opencadd` as a thin KLIFS wrapper.

**Rationale:** BioPython is the standard Python library for PDB file manipulation. It provides `PDBParser`, `NeighborSearch` (KD-tree backed, efficient for 6 Å queries), and `Superimposer`. KLIFS defines a canonical 85-residue binding site numbering system for all kinases — position 45 is universally the gatekeeper residue. Using KLIFS ensures the VRK1/EGFR comparison is made at structurally equivalent positions rather than sequence-position equivalents (which diverge significantly between distantly related kinases). `opencadd` (`pip install opencadd`) provides a clean DataFrame interface to the KLIFS REST API.

**Alternatives considered:**
- Sequence-only alignment (BioPython `pairwise2`): insufficient precision for structurally distant kinases like VRK1 (atypical kinase) vs. EGFR (receptor tyrosine kinase); KLIFS structural alignment is more reliable.
- PyMOL scripting: requires a GUI or PyMOL Python API license — not suitable for a headless pipeline.
- Schrodinger/MOE: commercial, not available.

**Reference structures:**

| Role | PDB ID | Chain | Ligand | Resolution | Notes |
|------|--------|-------|--------|------------|-------|
| VRK1 holo | **6AC9** | A | AMP-PNP (`ANP`) | 2.07 Å | Best resolution VRK1 ATP-analog complex |
| EGFR reference | **1M17** | A | Erlotinib (`AQ4`) | 2.60 Å | Wild-type EGFR, active conformation, most cited |
| VRK2 (if KLIFS has it) | TBD | — | — | — | Fall back to AlphaFold AF-O95551-F1 |

**Gatekeeper residues (confirmed from literature):**
- VRK1: **Met131** — confirmed in Couñago et al. 2017 (PMC6746079) and structure 6AC9
- EGFR: **Thr790** (wild-type) — canonical gatekeeper; T790M is the resistance mutation
- VRK2: Estimated **Cys213** based on sequence homology; confirm via KLIFS

**KLIFS API endpoints used:**
- `GET /structures_pdb_list?pdb-codes={PDB_ID}` → get KLIFS structure ID
- `GET /interactions_match_residues?structure_ID={id}` → 85-residue position mapping
- Position 45 = gatekeeper universally

**Binding site extraction approach:**
- Use `Bio.PDB.NeighborSearch` with cutoff 6.0 Å from any ligand heavy atom
- Filter HETATM by residue name (ligand HET code) and `res.id[0].startswith("H_")`
- Exclude water (`res.id[0] == "W"`) and standard solvent ions

**Superimposition approach:**
- Align on shared Cα atoms at KLIFS positions present in both structures
- Use `Bio.PDB.Superimposer.set_atoms(fixed_atoms, moving_atoms)`
- Report RMSD as alignment quality metric

**Key gotchas identified:**
1. Disordered atoms: explicitly select conformer A before distance calculations
2. HETATM encoding: `res.id[0]` is `"H_ANP"` not `"H"` — match by startswith
3. VRK1 is atypical kinase — KLIFS coverage may be limited; fall back to PKA (1ATP) sequence alignment
4. `Superimposer.set_atoms` requires exactly equal-length lists — drop unshared KLIFS positions from both before calling

---

## New Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `biopython` | ≥ 1.83 | PDB parsing, NeighborSearch, Superimposer |
| `opencadd` | ≥ 0.8 | KLIFS REST API wrapper (kinase pocket residue mapping) |
| `numpy` | ≥ 1.24 | Already available via pandas/rdkit; needed explicitly by BioPython |

Both packages are `pip install`-able into the existing environment. Add to `pyproject.toml` dependencies.

---

## NEEDS CLARIFICATION — All Resolved

No open questions remain. All technical decisions documented above.
