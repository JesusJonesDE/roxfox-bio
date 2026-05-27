# Feature Specification: VRK1 Oncology Evidence — DepMap Dependency + Structural Selectivity

**Feature Branch**: `003-vrk1-oncology-evidence`

**Created**: 2026-05-27

**Status**: Draft

**Input**: Build two parallel computational analyses to ground the VRK1 oncology thesis: (1) DepMap cancer dependency analysis identifying which cancer types are specifically dependent on VRK1 for survival, and (2) structural comparison of VRK1 vs. EGFR ATP binding sites to explain why lead scaffold SCF-013 is selective.

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Cancer Dependency Landscape (Priority: P1)

A researcher running the pipeline wants to know which cancer types specifically depend on VRK1 for survival, so they can identify the right patient population and indication for the RXF-001 program.

They run a single command and receive a ranked list of cancer lineages by VRK1 dependency strength, with classification of each lineage as strongly dependent, moderately dependent, or not dependent. The output includes a plain-language interpretation that can go directly into the research dossier.

**Why this priority**: This is the single most important piece of evidence for the oncology thesis. Without it, the cancer indication is a hypothesis. With it, the program has a data-backed patient selection rationale. It directly determines the clinical development path and the investor pitch.

**Independent Test**: Running the depmap analysis for VRK1 produces a ranked cancer lineage report covering at least 15 lineages with quantified dependency scores. The report is complete and self-contained — it can be read and acted on without any other pipeline output.

**Acceptance Scenarios**:

1. **Given** the pipeline has been set up and internet access is available, **When** the user runs the depmap analysis for VRK1, **Then** a ranked table of cancer lineages with VRK1 dependency scores is produced showing which lineages have the strongest dependency
2. **Given** VRK1 dependency scores are retrieved, **When** the analysis completes, **Then** the output distinguishes between cancer types where VRK1 is selectively essential, universally essential (pan-essential), and not essential
3. **Given** the analysis completes, **When** the researcher reviews the report, **Then** the top 3 most dependent cancer lineages are highlighted with a plain-language interpretation
4. **Given** DepMap data is already cached, **When** the command is re-run within the cache freshness window, **Then** the analysis completes from cache without re-fetching any data
5. **Given** DepMap is temporarily unavailable, **When** the command is run, **Then** a clear error message is shown and no partial or corrupt output is written

---

### User Story 2 — Structural Selectivity Explanation (Priority: P1)

A medicinal chemist wants to understand why SCF-013 is selective for VRK1 over EGFR — the dominant off-target — so they can design new analogs that preserve or improve that selectivity.

They run a single command and receive a residue-level comparison of the VRK1 and EGFR ATP binding sites, the gatekeeper residue for both kinases, and a plain-language structural hypothesis about which feature of SCF-013 exploits the binding site difference.

**Why this priority**: The entire chemistry strategy for RXF-001 depends on this. Without knowing the structural basis for SCF-013's selectivity, analog design is trial-and-error. With it, the medicinal chemistry program has a rational hypothesis — we know which vector to modify and which interactions to preserve.

**Independent Test**: Running the structural alignment command produces a binding site residue comparison table, identifies the gatekeeper residue for both kinases, and generates a report with a structural interpretation. The report is actionable for a medicinal chemist without any additional analysis.

**Acceptance Scenarios**:

1. **Given** PDB structures for VRK1 are available in the pipeline cache, **When** the structural alignment runs, **Then** the ATP binding site residues of VRK1 and EGFR are extracted, aligned by structural position, and differences are listed in a table
2. **Given** the alignment completes, **When** selectivity determinants are identified, **Then** the gatekeeper residue for both VRK1 and EGFR is explicitly reported with its amino acid identity
3. **Given** selectivity determinants are identified, **When** the report is generated, **Then** it includes a hypothesis about which structural feature of SCF-013 engages VRK1-specific residues
4. **Given** VRK1 has multiple crystal structures in the cache, **When** the best structure is selected, **Then** the selection criteria (resolution, ligand-bound status) are documented in the report
5. **Given** the EGFR reference structure is not yet downloaded, **When** the command runs for the first time, **Then** it automatically fetches a well-characterised EGFR inhibitor co-crystal and caches it

---

### User Story 3 — VRK1 vs. VRK2 Selectivity Profile (Priority: P2)

A program scientist wants to understand whether new VRK1 analogs will also inhibit VRK2 (the closest paralog), since the 2019 source paper for SCF-013 was designed for dual VRK1/2 activity. They need to know whether VRK2 selectivity is achievable and whether it is strategically desirable before the medicinal chemistry program begins.

**Why this priority**: This directly affects IP strategy and clinical positioning. A VRK1-selective compound has a different patent landscape and safety profile than a dual VRK1/2 compound. The answer shapes which analogs to synthesise first.

**Independent Test**: Running the structural alignment with VRK2 included produces a three-way binding site comparison (VRK1 / VRK2 / EGFR) showing how the three kinases differ, enabling a decision on whether selective vs. dual inhibition is feasible.

**Acceptance Scenarios**:

1. **Given** a VRK2 structural model is available (crystal structure or AlphaFold), **When** the alignment runs with VRK2 included, **Then** a three-way residue comparison table is produced covering all three kinases
2. **Given** the three-way comparison is complete, **When** binding site differences are analysed, **Then** the report indicates whether VRK1/VRK2 differentiation is structurally feasible based on binding site divergence
3. **Given** no VRK2 crystal structure is available in PDB, **When** the command runs, **Then** the AlphaFold model is downloaded automatically and its use is documented with pLDDT confidence scores

---

### User Story 4 — Integrated Oncology Report Update (Priority: P2)

After both analyses complete, the researcher wants the master research report (`data/results/research_report.md`) updated to include the new oncology evidence, so all findings are consolidated for stakeholder sharing.

**Why this priority**: The research report is the primary deliverable for funding discussions. New data must flow into it without manual assembly.

**Independent Test**: After both analyses complete, re-running the pipeline report command for VRK1 produces an updated research report that includes the DepMap lineage rankings and the structural selectivity interpretation in the VRK1 section.

**Acceptance Scenarios**:

1. **Given** both analyses have completed, **When** the research report is regenerated, **Then** a new oncology evidence section appears with the DepMap top-lineage findings
2. **Given** the structural analysis is complete, **When** the report is generated, **Then** the gatekeeper residue comparison and selectivity hypothesis appear in the VRK1 chemistry strategy section

---

### Edge Cases

- What happens when DepMap does not contain VRK1 in its gene list? → Report clearly states data unavailability; no silent empty output; suggests checking gene symbol alias
- What happens when DepMap classifies VRK1 as pan-essential (kills >70% of all cell lines)? → This is explicitly handled as a negative finding for oncology selectivity; the report explains implications for a therapeutic window
- What happens when the best VRK1 crystal structure has no bound ligand? → Falls back to apo structure with a documented warning; analysis continues using the protein structure
- What happens when binding site residues cannot be confidently aligned due to low structural similarity? → Reports alignment quality score; low-confidence positions are flagged and excluded from selectivity interpretation
- What happens when the DepMap API rate-limits requests? → Exponential backoff retry; graceful failure with partial results if retries are exhausted after 3 attempts
- What happens when VRK2 has no crystal structure and AlphaFold confidence is low in the binding site region? → pLDDT scores are reported; low-confidence residues (pLDDT < 70) are excluded from binding site comparison with an explicit warning

---

## Requirements *(mandatory)*

### Functional Requirements

**DepMap Analysis**

- **FR-001**: The system MUST retrieve VRK1 CRISPR gene effect scores (Chronos method) across all available cancer cell lines from the public DepMap dataset
- **FR-002**: The system MUST group cell lines by cancer lineage and compute per-lineage summary statistics: median gene effect, mean gene effect, count of strongly dependent lines (gene effect ≤ −0.5), and percentage of lineage that is strongly dependent
- **FR-003**: The system MUST classify each lineage into a dependency tier: strongly dependent (median ≤ −0.5), moderately dependent (−0.5 to −0.3), weakly dependent (−0.3 to −0.1), or not essential (> −0.1)
- **FR-004**: The system MUST detect and explicitly report pan-essential behaviour, defined as VRK1 being strongly dependent in more than 70% of all screened cell lines regardless of lineage
- **FR-005**: The system MUST produce a ranked output ordered by median VRK1 gene effect (most negative first), covering all lineages with at least 3 cell lines in the dataset
- **FR-006**: The system MUST cache retrieved DepMap data using the existing pipeline cache mechanism, respecting the configured freshness window
- **FR-007**: The system MUST produce a per-target markdown report with the lineage ranking table, dependency tier summary, and a plain-language interpretation of which cancer types are most promising as indications

**Structural Comparison**

- **FR-008**: The system MUST select the highest-resolution ligand-bound VRK1 crystal structure from the PDB structures already fetched and cached in the pipeline
- **FR-009**: The system MUST download a reference EGFR inhibitor co-crystal structure (erlotinib or gefitinib complex) from the RCSB PDB as a one-time fetch; it does not need to be part of the standard EGFR pipeline target fetch
- **FR-010**: The system MUST extract binding site residues for both kinases as all residues with at least one heavy atom within 6 Å of the bound ligand
- **FR-011**: The system MUST produce a residue-level comparison table aligned by structural position, showing for each position: the VRK1 residue, the EGFR residue, whether they are identical, and the type of difference (steric / charge / hydrogen bonding / none)
- **FR-012**: The system MUST identify and explicitly label the gatekeeper residue for both VRK1 and EGFR by its standard position in the kinase domain
- **FR-013**: The system MUST highlight all positions where VRK1 and EGFR differ and classify each as a candidate selectivity determinant or a non-determinant based on its proximity to the adenine-binding region
- **FR-014**: The system MUST optionally include VRK2 in the comparison when its structure or AlphaFold model is available, producing a three-way alignment table
- **FR-015**: The system MUST produce a markdown report with the binding site comparison, gatekeeper identification, and a selectivity hypothesis paragraph that a medicinal chemist can act on
- **FR-016**: Both analyses MUST be exposed as named CLI commands following existing pipeline conventions: `--target`, `--all`, `--force`, `--data-dir`, `--max-age` flags
- **FR-017**: Both analyses MUST be runnable independently and simultaneously without conflict

### Key Entities

- **DepMap Gene Effect Score**: A numerical score per gene per cell line representing the growth effect of knockout (strongly negative = essential for survival, near zero = dispensable). Source: Broad Institute DepMap, CRISPR Chronos method.
- **Cancer Lineage**: A grouping of cell lines by tissue of origin (e.g., Breast, Lung, Colorectal) used to aggregate dependency scores across biologically related lines
- **Dependency Tier**: Classification of a lineage's reliance on VRK1 (strongly dependent / moderately dependent / weakly dependent / not essential / pan-essential)
- **ATP Binding Site**: The set of protein residues within 6 Å of the bound ATP or ATP-mimic ligand in a kinase crystal structure; the primary site targeted by kinase inhibitors
- **Gatekeeper Residue**: A single amino acid at the back of the kinase ATP site whose size and identity is the primary determinant of inhibitor selectivity between kinase family members
- **Binding Site Alignment**: A position-by-position correspondence between the ATP binding site residues of two kinases, enabling direct structural comparison
- **Selectivity Determinant**: A residue position where VRK1 and EGFR differ in a chemically meaningful way (size, charge, hydrogen bonding capacity) that a selective inhibitor can exploit

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The DepMap analysis covers all major cancer lineages with at least 3 cell lines (expected ≥ 20 lineages) and completes within 5 minutes from a cold start on a standard laptop with internet access
- **SC-002**: The DepMap output unambiguously identifies the top indication hypothesis: at least one lineage with median gene effect ≤ −0.5 is reported, OR pan-essential behaviour is confirmed — in either case the result is actionable with no ambiguity
- **SC-003**: The structural alignment report identifies the gatekeeper residue for both VRK1 and EGFR and lists at least 3 positions where the binding sites differ, enabling a structural selectivity hypothesis to be stated in plain language
- **SC-004**: Both analyses complete successfully from a clean run within 15 minutes on a standard internet connection
- **SC-005**: Re-running either analysis within the cache freshness window completes within 30 seconds
- **SC-006**: The markdown reports produced by both analyses are self-contained — a reader with a drug discovery background can understand the key findings without referencing external sources
- **SC-007**: The research report is automatically updated with the new findings when the report stage is re-run after both analyses complete

---

## Assumptions

- DepMap CRISPR Chronos gene effect scores for VRK1 are available in the public DepMap dataset (VRK1 is a protein-coding gene present in most DepMap releases)
- The existing pipeline cache system is used for DepMap data storage without structural changes
- The best-resolution VRK1 crystal structure is already in the pipeline PDB cache from a prior `pipeline fetch` run; no additional VRK1 structure downloads are required
- A standard erlotinib- or gefitinib-bound EGFR co-crystal is sufficient as a selectivity reference; EGFR is one of the best-characterised kinases and has hundreds of published structures
- BioPython is installable in the pipeline environment for PDB structure parsing
- VRK2 comparison is best-effort: if no crystal structure exists, the AlphaFold model is used and low-confidence regions are excluded
- Both analyses write to `data/results/<GENE>/` and do not modify any existing cached API data
- The analyses are read-only with respect to ChEMBL, Open Targets, and other existing pipeline data sources
- Both modules follow the same CLI flag conventions as existing pipeline commands (`--target`, `--all`, `--force`)
- DepMap API access does not require authentication for the public gene effect dataset
