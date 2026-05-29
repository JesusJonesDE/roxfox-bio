# Feature Specification: IGHMBP2 Fragment-Based Virtual Screening

**Feature Branch**: `006-ighmbp2-fragment-screen`

**Created**: 2026-05-29

**Status**: Draft

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Pocket Identification (Priority: P1)

A researcher needs to know where small molecules could bind on the IGHMBP2 protein before running any screening. The system analyses the cached AlphaFold2 structure of IGHMBP2, identifies all druggable pockets using geometric cavity detection, ranks them by druggability score, and selects the ATP/helicase binding site as the primary docking pocket. The pocket coordinates, volume, and druggability score are reported.

**Why this priority**: Without a defined pocket, no docking can occur. All downstream steps depend on this. The AlphaFold2 structure is already cached — this step is fast and foundational.

**Independent Test**: Running `pipeline fragment --target IGHMBP2 --step pocket` identifies at least one pocket with volume > 300 Å³ and writes `pocket_analysis.json` with centroid coordinates within 60 seconds.

**Acceptance Scenarios**:

1. **Given** the AlphaFold2 IGHMBP2 structure is cached, **When** pocket detection runs, **Then** at least one pocket with volume > 300 Å³ is identified and ranked
2. **Given** multiple pockets are found, **When** the ATP site is identifiable by proximity to conserved helicase motifs (Walker A/B), **Then** it is selected as the primary pocket
3. **Given** no pocket exceeds the minimum volume threshold, **Then** the threshold is relaxed to 200 Å³ with a warning in the report
4. **Given** the pocket is identified, **Then** centroid coordinates, box dimensions, and druggability score are written to `pocket_analysis.json`

---

### User Story 2 — Fragment Library Download and Preparation (Priority: P1)

A researcher needs a curated set of fragment-like small molecules to screen against the pocket. The system downloads a subset of the ZINC fragment library (MW 100–250, logP ≤ 3, ≤ 5 heavy atoms non-hydrogen, passes Ro3), deduplicates by scaffold, and prepares a local SMILES file of 5,000–10,000 fragments ready for docking.

**Why this priority**: Without the fragment library, no screening can happen. This runs once and is cached — subsequent runs skip the download.

**Independent Test**: Running `pipeline fragment --target IGHMBP2 --step library` downloads and caches a fragment library of ≥ 5,000 compounds in under 10 minutes, with all compounds passing fragment-likeness filters (MW ≤ 250, ≤ 3 HBD, ≤ 3 HBA, logP ≤ 3).

**Acceptance Scenarios**:

1. **Given** the library is not cached, **When** the step runs, **Then** it downloads from ZINC and stores locally
2. **Given** the library is already cached, **When** the step runs, **Then** it skips the download and uses the cached file
3. **Given** downloaded compounds, **When** filtered by fragment-likeness (Rule of Three), **Then** all compounds in the output pass: MW ≤ 250, HBD ≤ 3, HBA ≤ 3, logP ≤ 3, rotatable bonds ≤ 3
4. **Given** the library contains duplicates by scaffold, **When** deduplication runs, **Then** only one representative per scaffold is retained

---

### User Story 3 — Fragment Docking Screen (Priority: P1)

A researcher runs all fragments from the library against the identified IGHMBP2 pocket and receives a ranked list of the top-scoring hits. The system docks each fragment into the pocket using the existing Vina docking infrastructure and returns the top 50 fragments by predicted binding affinity.

**Why this priority**: This is the core scientific value — identifying fragments that fit the IGHMBP2 pocket. The M1 Max can run ~5,000 fragments overnight at low exhaustiveness.

**Independent Test**: Running `pipeline fragment --target IGHMBP2 --step dock` docks ≥ 5,000 fragments and returns 50 hits with affinity ≤ −5.0 kcal/mol, completing within 12 hours on local hardware.

**Acceptance Scenarios**:

1. **Given** a fragment library and pocket definition, **When** docking runs, **Then** each fragment is docked and a score recorded
2. **Given** all fragments are docked, **When** results are ranked, **Then** the top 50 by affinity (most negative) are retained
3. **Given** a fragment fails to dock (SMILES cannot be parametrised), **Then** it is skipped with a log entry; the screen continues
4. **Given** the screen was partially completed and interrupted, **When** it resumes, **Then** already-docked fragments are not re-docked (cache-based resumption)
5. **Given** the screen completes, **Then** results are written to `fragment_hits.csv` with columns: fragment_id, smiles, affinity_kcal_mol, pose_file

---

### User Story 4 — Fragment Hit Clustering (Priority: P2)

A researcher wants to understand the chemical diversity of the top fragment hits — whether they all represent the same binding mode or distinct chemotypes. The system clusters the top 50 fragments by Tanimoto similarity, identifies cluster representatives, and annotates hits with their cluster membership.

**Why this priority**: Without clustering, a researcher might grow 10 analogs of effectively the same fragment and miss 3 distinct binding modes. Clustering is fast and significantly increases the value of the subsequent growing step.

**Independent Test**: Running `pipeline fragment --target IGHMBP2 --step cluster` produces a clustering table with ≥ 3 distinct clusters from the top 50 hits and writes `fragment_clusters.csv`.

**Acceptance Scenarios**:

1. **Given** top 50 fragment hits, **When** clustering runs at Tanimoto similarity threshold 0.6, **Then** each fragment is assigned a cluster ID
2. **Given** clusters, **When** representatives are selected, **Then** the fragment with the best affinity in each cluster is the representative
3. **Given** clustering results, **Then** `fragment_clusters.csv` contains: fragment_id, smiles, affinity, cluster_id, is_representative

---

### User Story 5 — Fragment Growing into Drug-Like Molecules (Priority: P2)

A researcher takes the top cluster representative fragments and grows them computationally into larger, drug-like molecules (MW 300–450) that are predicted to retain and improve upon the fragment's binding mode. The system uses structure-based growing: extending the fragment into unexplored regions of the pocket identified from the docking pose.

**Why this priority**: Fragments themselves are too small to be drugs. Growing produces candidates that can be synthesised and tested. This is the step that generates the actual hit list for IGHMBP2.

**Independent Test**: Running `pipeline fragment --target IGHMBP2 --step grow` takes the top 10 cluster representatives and produces ≥ 20 grown molecules with MW 300–450, passing Ro5, written to `grown_candidates.csv`.

**Acceptance Scenarios**:

1. **Given** top cluster representatives with docking poses, **When** growing runs, **Then** each fragment is extended with drug-like substituents into available pocket space
2. **Given** grown molecules, **When** filtered, **Then** all pass: MW 300–450, Ro5 compliant, synthetic accessibility score < 4 (estimated)
3. **Given** grown molecules, **Then** the parent fragment ID is recorded for traceability
4. **Given** fewer than 5 cluster representatives, **When** growing runs, **Then** it uses all available representatives without error

---

### User Story 6 — ADMET Screening and Final Ranking (Priority: P2)

A researcher receives a final ranked list of grown candidates with ADMET predictions, ready to drop directly into the full dock → validate pipeline. The system runs the existing ADMET gate on all grown candidates and produces `compounds_filtered.csv` in the same schema as VRK1, enabling immediate use of `pipeline dock` and `pipeline validate`.

**Why this priority**: ADMET filtering removes grown candidates with liabilities before any resource is invested in further computation or synthesis. The output in the VRK1 schema means zero additional pipeline work to continue.

**Independent Test**: After running the full fragment pipeline, `data/results/IGHMBP2/compounds_filtered.csv` exists with ≥ 10 rows in the same column schema as `data/results/VRK1/compounds_filtered.csv`, ready for `pipeline dock --target IGHMBP2 --all-scaffolds`.

**Acceptance Scenarios**:

1. **Given** grown candidates, **When** ADMET runs, **Then** each candidate receives BBB, CYP, solubility, and bioavailability scores
2. **Given** ADMET results, **When** final ranking is produced, **Then** candidates are sorted by composite score (docking affinity × ADMET pass rate)
3. **Given** the final ranked list, **Then** `compounds_filtered.csv` matches the VRK1 column schema exactly
4. **Given** `compounds_filtered.csv` is written, **Then** `pipeline dock --target IGHMBP2 --all-scaffolds` runs without error
5. **Given** the full pipeline completes, **Then** `fragment_screen_report.md` is written summarising pocket, top fragments, clusters, grown candidates, and ADMET results

---

### User Story 7 — End-to-End Pipeline Command (Priority: P1)

A researcher runs the entire fragment screening workflow with a single command and returns the next morning to a ranked candidate list. The system executes all steps in sequence (pocket → library → dock → cluster → grow → ADMET) with progress reporting and cache-based resumption.

**Why this priority**: The individual steps have scientific value, but the researcher's primary need is a single command that produces candidates without babysitting each step.

**Independent Test**: Running `pipeline fragment --target IGHMBP2` on a clean cache completes all steps overnight (< 14 hours) and writes all output files including `compounds_filtered.csv` with ≥ 10 candidates.

**Acceptance Scenarios**:

1. **Given** a clean cache, **When** `pipeline fragment --target IGHMBP2` runs, **Then** all 6 steps execute in sequence automatically
2. **Given** the pipeline is interrupted mid-way, **When** it is re-run, **Then** completed steps are skipped and it resumes from the first incomplete step
3. **Given** any step fails, **Then** the pipeline stops, reports which step failed and why, and does not silently continue to the next step
4. **Given** `--force` is passed, **Then** all steps re-run from scratch, ignoring cached results
5. **Given** `--step <name>` is passed, **Then** only that single step runs

---

### Edge Cases

- What if the AlphaFold2 structure has poor pLDDT scores in the pocket region (< 70)? → Flag with warning in report; docking results are flagged as lower confidence but not blocked
- What if ZINC is unavailable during library download? → Fall back to a bundled minimal fragment library (500 fragments) and warn the user
- What if all 5,000 fragments fail to dock? → Pipeline stops at step 3 with a clear error describing the failure mode
- What if the growing step produces 0 drug-like molecules? → Report the failure; provide top fragments directly in `compounds_filtered.csv` as fallback
- What if the pocket centroid falls in a disordered loop? → Use the second-ranked pocket automatically, report the substitution

---

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST identify at least one druggable binding pocket (volume > 200 Å³) in the IGHMBP2 AlphaFold2 structure
- **FR-002**: The system MUST download and cache a fragment library of ≥ 5,000 compounds passing fragment-likeness criteria (MW ≤ 250, logP ≤ 3, HBD ≤ 3, HBA ≤ 3, rotatable bonds ≤ 3)
- **FR-003**: The system MUST dock all library fragments against the identified pocket and rank by affinity
- **FR-004**: The system MUST cluster top-50 fragment hits by chemical similarity and identify cluster representatives
- **FR-005**: The system MUST grow top cluster representatives into drug-like molecules (MW 300–450, Ro5-compliant, SA score < 4)
- **FR-006**: The system MUST run ADMET predictions on all grown candidates
- **FR-007**: The system MUST write `compounds_filtered.csv` in the same schema as VRK1, ready for `pipeline dock`
- **FR-008**: The system MUST write `fragment_screen_report.md` summarising all steps, key findings, and candidate quality
- **FR-009**: Each step MUST cache its output and skip re-computation on re-run unless `--force` is passed
- **FR-010**: The pipeline MUST support `--step <name>` to run individual steps
- **FR-011**: The system MUST handle fragment docking failures gracefully — individual failures skip the fragment; pipeline continues
- **FR-012**: All computation MUST run locally on Apple M1 Max; no cloud services required
- **FR-013**: The complete pipeline MUST complete in under 14 hours on local hardware

### Key Entities

- **Pocket**: A druggable cavity in the IGHMBP2 structure with centroid coordinates, volume (Å³), surface area, and druggability score
- **Fragment**: A small molecule (MW ≤ 250) from the ZINC fragment library with SMILES, fragment ID, and Ro3 properties
- **FragmentHit**: A fragment with a docking result — affinity (kcal/mol), pose file path, and cluster assignment
- **Cluster**: A group of chemically similar fragment hits with a designated representative
- **GrownCandidate**: A drug-like molecule (MW 300–450) derived from a fragment hit by structure-based growing, with parent fragment ID, SMILES, affinity estimate, and ADMET scores
- **FragmentScreenResult**: The complete output of the pipeline — pocket definition, top hits, clusters, grown candidates — stored per target

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Pocket identification completes in under 60 seconds and identifies at least one pocket with volume > 200 Å³
- **SC-002**: Fragment library of ≥ 5,000 compounds is downloaded and filtered in under 10 minutes
- **SC-003**: 5,000 fragment docking runs complete in under 12 hours on Apple M1 Max
- **SC-004**: At least 50 fragments score ≤ −5.0 kcal/mol (threshold for meaningful fragment binding)
- **SC-005**: Fragment clustering identifies ≥ 3 distinct chemotype clusters from the top 50 hits
- **SC-006**: Fragment growing produces ≥ 20 grown candidates that pass Ro5 and SA score < 4
- **SC-007**: `compounds_filtered.csv` is produced with ≥ 10 rows in the exact VRK1 schema
- **SC-008**: `pipeline dock --target IGHMBP2 --all-scaffolds` runs without error after fragment pipeline completes
- **SC-009**: Complete pipeline runs end-to-end in under 14 hours on local hardware with no manual intervention

---

## Assumptions

- The AlphaFold2 IGHMBP2 structure (UniProt P38935) is already cached from `pipeline fetch`; the fragment pipeline reads it directly
- fpocket is already installed and accessible from the command line (confirmed from earlier pipeline work)
- The existing `pipeline dock` infrastructure (Vina, meeko, receptor PDBQT preparation) is reused for fragment docking
- The existing ADMET gate from spec-005 is reused unchanged for grown candidate scoring
- Fragment library sourced from ZINC fragments subset (freely downloadable, no registration required); a 500-fragment bundled fallback is included for offline use
- Fragment docking exhaustiveness is set to 4 (vs. 32 for leads) to enable ~5,000 runs in a reasonable time; this is standard practice for fragment screens
- Fragment growing uses RDKit reaction-based enumeration with a curated set of medicinal chemistry transformations (e.g., amide bond formation, N-alkylation, ring closures); not a de novo generative model
- Synthetic accessibility (SA) score is estimated computationally using RDKit SA score module — not experimentally validated
- The selectivity gate from spec-005 is not run on grown candidates at this stage (no off-target structures for IGHMBP2 helicase panel); selectivity assessment deferred to after experimental validation
- IGHMBP2 AlphaFold2 structure pLDDT reliability: the helicase core domain is typically high-confidence (pLDDT > 70); disordered linkers are excluded from pocket analysis
