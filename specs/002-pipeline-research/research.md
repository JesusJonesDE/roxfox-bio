# Research: Drug Discovery Pipeline Research Tool

**Branch**: `002-pipeline-research` | **Date**: 2026-05-26

---

## Decision 1: CLI Framework

**Decision**: Typer + Rich

**Rationale**: Typer is the modern standard for Python 3.11+ CLIs — uses type annotations instead of decorators, generates `--help` automatically, and pairs natively with Rich for progress/logging output. The pipeline has multiple subcommands (`fetch`, `analyze`, `report`, `run`) with shared flags (`--target`, `--force`) which Typer handles cleanly with a top-level `app` and subcommand callbacks.

**Alternatives considered**:
- Click: Battle-tested but more boilerplate; no native type hint integration
- argparse: Too low-level for a multi-subcommand pipeline
- Makefile: Not composable, no cache-skip logic

---

## Decision 2: HTTP Client

**Decision**: `httpx` with retry logic via `tenacity`

**Rationale**: `httpx` is async-capable, has a cleaner API than `requests`, and supports streaming response bodies — needed for paginated ChEMBL calls returning 10k+ records. `tenacity` provides exponential backoff on 429/503 responses cleanly without manual retry loops.

**Alternatives considered**:
- requests: Synchronous only; no streaming for large JSON bodies
- aiohttp: Full async overkill for a CLI tool that runs stages sequentially

---

## Decision 3: ChEMBL Data Access

**Decision**: Direct REST API calls (not `chembl_webresource_client`)

**Rationale**: The Python client library wraps REST but adds serialization overhead on large result sets and doesn't support streaming. Direct API gives full control over pagination and response handling.

**Key endpoints**:
- Target lookup: `GET https://www.ebi.ac.uk/chembl/api/data/target?search={UNIPROT_ID}&format=json`
- Bioactivity: `GET https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id={ID}&assay_type__in=B,F&activity_type__in=IC50,Ki,Kd&limit=1000&offset={N}&format=json`
- Pagination: max 1000 per page; `page_meta.total_count` gives total
- Returns per record: `molecule_chembl_id`, `canonical_smiles`, `standard_value`, `standard_type`, `standard_units`, `pchembl_value`

**Rate limits**: No hard limit; add ~0.5s delay between pages; exponential backoff on 429.

---

## Decision 4: Open Targets Access

**Decision**: GraphQL via `httpx` POST to `https://api.platform.opentargets.org/api/v4/graphql`

**Rationale**: Used successfully in prior research session. GraphQL allows fetching exactly the fields needed (genetic evidence scores, disease associations) in a single query. No authentication required.

**Key query fields**: `target.id`, `target.approvedSymbol`, `associatedDiseases.rows[].disease.name`, `associatedDiseases.rows[].score`, `target.tractability`

---

## Decision 5: PDB Structure Inventory

**Decision**: RCSB PDB REST Search API (`https://search.rcsb.org/rcsbsearch/v1/query`) + summary endpoint

**Rationale**: RCSB's JSON search API supports querying by UniProt accession with filters on resolution, experimental method, and ligand presence. Returns PDB IDs which can then be resolved to metadata via `https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}`.

**Key fields returned**: PDB ID, resolution, experimental method (X-ray/Cryo-EM/NMR), R-free, ligand IDs (from `rcsb_binding_affinity` or `nonpolymer_entities`)

---

## Decision 6: AlphaFold Access

**Decision**: EBI AlphaFold REST API — `GET https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ID}`

**Rationale**: Simple REST call returns model metadata including pLDDT scores and download URLs. No authentication required. Returns array of models; take the first (canonical) entry.

---

## Decision 7: ClinicalTrials.gov Access

**Decision**: ClinicalTrials.gov API v2 — `GET https://clinicaltrials.gov/api/v2/studies`

**Rationale**: v2 is the current stable API (v1 deprecated). Query by gene/protein name using the `query.term` parameter. Returns study metadata including phase, status, interventions, and sponsor.

---

## Decision 8: Cheminformatics — RDKit

**Decision**: RDKit installed via conda-forge (`conda install -c conda-forge rdkit`)

**Rationale**: RDKit is the industry standard for cheminformatics in Python. Conda-forge provides the most stable, up-to-date builds. pip-based `rdkit` wheels exist but conda is preferred for dependency management.

**Lipinski Ro5 pattern**:
```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

mol = Chem.MolFromSmiles(smiles)
mw = Descriptors.MolWt(mol)
logp = Crippen.MolLogP(mol)
hbd = Descriptors.NumHDonors(mol)
hba = Descriptors.NumHAcceptors(mol)
violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
passes_ro5 = violations < 2
```

**Murcko scaffold pattern**:
```python
from rdkit.Chem.Scaffolds import MurckoScaffold

scaffold_smi = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smiles, includeChirality=False)
```

**Gotcha**: Some ChEMBL SMILES fail RDKit sanitization. Handle with try/except and log failures; skip the compound rather than crashing the pipeline.

---

## Decision 9: Caching Architecture

**Decision**: Per-source JSON files in `data/cache/{target}/{source}_{timestamp}.json` + a `pipeline_manifest.json` for stage completion tracking

**Rationale**: Simple, inspectable, no extra dependencies. Each cached file records its fetch timestamp. The manifest tracks which `{stage}:{target}` combinations have been completed and when. Freshness check: compare manifest timestamp to `--max-age` setting (default 30 days).

**Alternative considered**: `diskcache` library — adds a dependency for no meaningful benefit at this scale; SQLite overkill for a single-researcher tool.

---

## Decision 10: Packaging

**Decision**: `pyproject.toml` with `[project.scripts]` entry point (`pipeline = "pipeline.cli:app"`)

**Rationale**: Modern Python packaging standard. Install with `pip install -e .` in a conda/venv environment. No build step needed for development. Entry point makes `pipeline` available as a shell command after install.

---

## Decision 11: Data Output Format

**Decision**: CSV for tabular data (compounds, scaffolds, structures) + Markdown for dossiers

**Rationale**: CSV is universally readable in Excel, pandas, and any text editor — no special tooling needed. Markdown dossiers are readable in any editor and render on GitHub. No database needed.

---

## Summary: Full Tech Stack

| Layer | Choice | Version |
|-------|--------|---------|
| Language | Python | 3.11+ |
| CLI | Typer + Rich | latest stable |
| HTTP | httpx + tenacity | latest stable |
| Cheminformatics | RDKit | conda-forge latest |
| Data processing | pandas | latest stable |
| Caching | Filesystem JSON | built-in |
| Packaging | pyproject.toml | PEP 517 |
| Data sources | ChEMBL REST, Open Targets GraphQL, RCSB REST, AlphaFold EBI REST, ClinicalTrials v2 | public APIs |
