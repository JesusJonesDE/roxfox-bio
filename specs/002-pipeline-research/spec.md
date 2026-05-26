# Feature Specification: Drug Discovery Pipeline Research Tool

**Feature Branch**: `002-pipeline-research`

**Created**: 2026-05-26

**Status**: Draft

**Input**: Python CLI pipeline for deep computational research across three drug programs — RXF-001 (VRK1), RXF-002 (IGHMBP2), RXF-003 (VCP) — querying public biomedical databases, producing per-target research dossiers with compound data, structural inventory, and competitive landscape.

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Run Full Research Pipeline for a Target (Priority: P1)

A researcher runs the complete pipeline for one or all targets from the command line. The tool fetches all relevant public data, applies filters and analysis, and writes structured outputs to disk. On a re-run, already-fetched and still-valid data is reused without re-querying. A `--force` flag triggers fresh fetches for specific stages or targets.

**Why this priority**: This is the core value of the tool — replacing manual API queries with a single reproducible command. Everything else builds on this.

**Independent Test**: Run the pipeline against one target (e.g., VRK1). Confirm output directories are created, data files are written, and a markdown dossier exists. Re-run without `--force` and confirm no new network requests are made for cached stages.

**Acceptance Scenarios**:

1. **Given** no prior data exists, **When** the user runs `pipeline run --target VRK1`, **Then** all stages execute in order, data is fetched from all configured sources, outputs are written to the VRK1 directory tree, and a completion summary is printed.
2. **Given** a prior successful run exists for VRK1, **When** the user runs `pipeline run --target VRK1` again, **Then** all cached stages are skipped and the run completes without network calls for those stages.
3. **Given** a prior run exists, **When** the user runs `pipeline run --target VRK1 --force fetch`, **Then** only the fetch stage re-executes; analysis and report stages use the freshly fetched data.
4. **Given** no target is specified, **When** the user runs `pipeline run --all`, **Then** all three targets (VRK1, IGHMBP2, VCP) are processed sequentially.

---

### User Story 2 — Inspect Compound Data and Scaffold Analysis (Priority: P2)

After running the pipeline, the researcher reviews the filtered compound list and scaffold clusters for a target to identify the best chemical starting points. They open a CSV of drug-like compounds ranked by potency, and a scaffold summary showing which Murcko frameworks appear most frequently and in which potency range.

**Why this priority**: The compound and scaffold output is the primary scientific deliverable — it informs what chemistry exists and what IP space may be open.

**Independent Test**: After running the pipeline for VRK1, open `results/VRK1/compounds_filtered.csv` and confirm it contains compounds with IC50 < 10µM, molecular weight, logP, HBD, HBA, rotatable bonds, and Lipinski pass/fail. Open `results/VRK1/scaffolds.csv` and confirm Murcko scaffolds are listed with compound counts and median potency.

**Acceptance Scenarios**:

1. **Given** the fetch and analysis stages have completed for a target, **When** the user opens the compounds CSV, **Then** it contains only compounds passing the potency threshold (IC50/Ki/Kd ≤ 10µM), with all Lipinski properties calculated and a pass/fail flag per compound.
2. **Given** the analysis stage has completed, **When** the user opens the scaffold summary, **Then** each unique Murcko scaffold is listed with the number of compounds sharing it and the median potency of that cluster.
3. **Given** a target with no bioactivity data in the database, **When** the pipeline runs, **Then** the compounds CSV is written with headers only and the dossier documents the data gap explicitly.

---

### User Story 3 — Review Structural Data Inventory (Priority: P2)

The researcher checks which experimentally determined protein structures and predicted models are available for each target, including resolution quality and whether a ligand is co-crystallised (indicating a known binding pocket). This informs whether structure-based approaches are feasible.

**Why this priority**: Structural availability directly impacts the next research phase and is a key credibility signal for VC and SAB conversations.

**Independent Test**: After running the pipeline for IGHMBP2, open `results/IGHMBP2/structures.csv` and confirm it lists PDB entries with resolution, experimental method, and ligand presence, plus the AlphaFold model accession and confidence score summary.

**Acceptance Scenarios**:

1. **Given** experimental structures exist for a target, **When** the structural stage runs, **Then** the structures CSV lists each entry with: structure ID, resolution (Å), experimental method, whether a small-molecule ligand is co-crystallised, and chain identifiers.
2. **Given** no experimental structures exist for a target, **When** the structural stage runs, **Then** the predicted model is retrieved and recorded with confidence metadata; the dossier notes the absence of experimental data.
3. **Given** both experimental structures and a predicted model exist, **When** the structural stage runs, **Then** both are recorded and the dossier identifies the highest-quality structure available for downstream work.

---

### User Story 4 — Read the Research Dossier (Priority: P3)

After the pipeline completes, the researcher reads a single markdown file per target that summarises all findings: genetic evidence, compound landscape, scaffold highlights, structural availability, selectivity liabilities, and competitive context. This document is shareable as an internal research memo and suitable for investor data rooms.

**Why this priority**: The dossier synthesises all data into a human-readable, shareable output — what gets presented to collaborators, SAB members, and investors.

**Independent Test**: After a full pipeline run for VCP, open `results/VCP/dossier.md` and confirm it contains all sections: Overview, Genetic Evidence, Bioactivity Summary, Scaffold Highlights, Structural Data, Selectivity Profile, Competitive Landscape, and Data Gaps / Limitations.

**Acceptance Scenarios**:

1. **Given** all pipeline stages have completed for a target, **When** the report stage runs, **Then** `dossier.md` is generated with all sections populated from the computed data.
2. **Given** a stage produced no data, **When** the dossier is generated, **Then** the relevant section states the absence of data with a brief explanation, rather than being omitted.
3. **Given** the dossier is generated, **When** it is read, **Then** it contains no raw API payloads, no code, and no implementation details — only scientific findings written in plain language.

---

### Edge Cases

- What happens when an external data source is unreachable? → The stage fails gracefully with a clear error message; any valid cached data from a prior run is preserved and the pipeline continues other stages if possible.
- What happens when a target returns zero compounds? → Pipeline completes; dossier documents the data gap; CSV is written with headers only.
- What happens when cached data is corrupt or incomplete? → Freshness/validity check detects the issue; the affected stage re-runs automatically and overwrites the corrupted files.
- What if the same compound appears in multiple assay types (IC50, Ki, Kd)? → All records are retained; the best (lowest) value per compound per assay type is used for filtering and ranking.
- What if a target has hundreds of thousands of bioactivity records? → Fetching is paginated; results are streamed to disk rather than held in memory; the run completes without memory issues.

---

## Requirements *(mandatory)*

### Functional Requirements

**Data Fetching**

- **FR-001**: The system MUST fetch all bioactivity records for each target from a public compound–target interaction database, covering IC50, Ki, and Kd assay types.
- **FR-002**: The system MUST retrieve genetic association evidence scores for each target from a population-level disease–gene resource.
- **FR-003**: The system MUST retrieve all experimentally determined protein structures for each target, including resolution, experimental method, and ligand co-crystallisation status.
- **FR-004**: The system MUST retrieve the AI-predicted structural model for each target when experimental structures are absent or low quality.
- **FR-005**: The system MUST retrieve competitive landscape data: approved drugs, clinical candidates, and active trials associated with each target.

**Caching & Idempotency**

- **FR-006**: The system MUST cache all fetched API responses to disk with a fetch timestamp; on subsequent runs, cached data MUST be reused if it is within the configured freshness window and passes a completeness check.
- **FR-007**: The system MUST provide a `--force` flag that accepts a stage name (`fetch`, `analyze`, `report`) or `all` to bypass the cache and re-execute for the specified scope.
- **FR-008**: The system MUST NOT overwrite existing processed result files unless the corresponding stage is explicitly re-run or forced.
- **FR-009**: The default cache freshness window MUST be 30 days and MUST be user-configurable.

**Analysis**

- **FR-010**: The system MUST filter bioactivity records to retain only those with a potency value ≤ 10,000 nM (10µM) for IC50, Ki, or Kd.
- **FR-011**: The system MUST calculate Lipinski Rule of Five properties for each compound passing the potency filter: molecular weight, calculated logP, hydrogen bond donor count, hydrogen bond acceptor count. Compounds failing two or more rules MUST be flagged.
- **FR-012**: The system MUST extract the Murcko scaffold for each drug-like compound and group compounds by scaffold, reporting cluster size and median potency per scaffold.
- **FR-013**: The system MUST profile selectivity by identifying known activity of each filtered compound against other human protein targets and flagging compounds with potent off-target activity (≤ 1µM) against more than three unrelated targets.

**Reporting**

- **FR-014**: The system MUST produce a filtered compound CSV per target containing: compound identifier, SMILES, potency value, assay type, Lipinski properties, Lipinski pass/fail flag, and Murcko scaffold identifier.
- **FR-015**: The system MUST produce a scaffold summary CSV per target containing: scaffold SMILES, scaffold identifier, compound count, and median and best potency values.
- **FR-016**: The system MUST produce a structures CSV per target listing experimental and predicted structures with quality metadata.
- **FR-017**: The system MUST produce a markdown research dossier per target synthesising findings from all stages.

**CLI Interface**

- **FR-018**: The system MUST be operable as a CLI with subcommands: `fetch`, `analyze`, `report`, and `run` (executes all stages in sequence).
- **FR-019**: The system MUST accept `--target <NAME_OR_ID>` to restrict execution to one target and `--all` to run all configured targets.
- **FR-020**: The system MUST print a human-readable progress log to stdout during execution, including stage name, target, cache hit/miss status, record counts, and outcome.
- **FR-021**: The system MUST support adding a new target via a configuration entry without requiring code changes.

### Key Entities

- **Target**: A protein drug target identified by gene name and UniProt ID, mapped to a program code (RXF-001/002/003) and one or more disease indications.
- **BioactivityRecord**: A measured interaction between a compound and a target — assay type (IC50/Ki/Kd), value in nM, source record identifier, and assay conditions.
- **Compound**: A small molecule with a database identifier and SMILES structure, with computed physicochemical properties and potency values.
- **Scaffold**: A Murcko framework extracted from a compound, used to group chemically related compounds into clusters.
- **Structure**: A protein structure entry (experimental or predicted) with resolution, experimental method, and ligand-binding metadata.
- **CacheEntry**: A stored API response including fetch timestamp, data source, target scope, stage name, record count, and a validity flag.
- **Dossier**: A per-target markdown document aggregating all analysis outputs into a human-readable research summary.

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A full pipeline run across all three targets completes in under 15 minutes on a standard laptop with a stable internet connection.
- **SC-002**: A re-run with all stages fully cached completes in under 30 seconds per target with zero external network calls.
- **SC-003**: The filtered compound list for VRK1 is non-empty and contains compounds consistent with known published ChEMBL data for that target.
- **SC-004**: Each markdown dossier is fully self-contained — all key findings are readable without opening any CSV or raw data file.
- **SC-005**: Adding a new research target requires only a single configuration entry; no code changes are needed.
- **SC-006**: A corrupted or missing cache file is detected and recovered automatically without user intervention in 100% of tested cases.

---

## Assumptions

- The three initial targets (VRK1, IGHMBP2, VCP) are pre-configured with their UniProt IDs; additional targets can be added via configuration without code changes.
- All external data sources are accessed via their public APIs with no authentication or API keys required.
- The tool runs on a single researcher's machine; no distributed processing, server, or cloud infrastructure is needed.
- Cache freshness defaults to 30 days; data older than this is considered stale and triggers a re-fetch on the next run.
- Selectivity profiling is limited to compounds already in the fetched bioactivity dataset — no separate proteome-wide query is needed.
- The tool does not perform molecular docking, 3D conformer generation, or any computation requiring local molecular simulation software.
- Veber rules (rotatable bonds and polar surface area) are an optional secondary filter and do not affect the primary filtered compound output.
- Output is stored on the local filesystem; no cloud storage, remote sync, or database is required.
- The researcher is the sole user; no authentication, role management, or multi-user support is in scope.
