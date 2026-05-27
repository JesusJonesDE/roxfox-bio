# Feature Specification: VRK1 Indication Validation and Chemistry Deepening

**Feature Branch**: `004-vrk1-indication-chem`

**Created**: 2026-05-27

**Status**: Draft

## User Scenarios & Testing *(mandatory)*

### User Story 1 — CNS/Brain Biomarker Analysis (Priority: P1)

A researcher runs a single command against the CNS/Brain cell lines that are strongly VRK1-dependent and gets back a report identifying which co-mutations (e.g. TP53, CDKN2A, ATRX) consistently appear in dependent lines but not in non-dependent ones. This becomes the proposed patient stratification biomarker for a first-in-human trial design.

**Why this priority**: CNS/Brain is the largest high-dependency lineage (90 cell lines, 83% strongly dependent). A clear biomarker transforms "VRK1 is important in brain cancer" into "patients with mutation X in glioblastoma should receive this drug" — the single most investor-legible output from this programme.

**Independent Test**: Can be run with just the DepMap cache already on disk plus a CCLE mutation data download. Delivers a standalone biomarker candidate table without any other story being complete.

**Acceptance Scenarios**:

1. **Given** the DepMap VRK1 dependency data is cached, **When** the user runs `pipeline biomarker --target VRK1 --lineage "CNS/Brain"`, **Then** a report is produced listing co-mutations ranked by enrichment score in strongly-dependent versus non-dependent lines.
2. **Given** the biomarker report exists, **When** no enriched mutations are found above threshold, **Then** the report clearly states "no significant co-mutation biomarker identified" rather than failing silently.
3. **Given** the CCLE mutation file is unavailable or has changed format, **When** the command runs, **Then** a clear error explains which data source is missing and how to obtain it.

---

### User Story 2 — VRK1 vs VRK2 Paralog Structural Comparison (Priority: P2)

A researcher runs the existing structural comparison command with a flag to include VRK2 and gets a three-way comparison table (VRK1 / VRK2 / EGFR). The report distinguishes positions that are "VRK1-specific handles" (unique to VRK1 relative to both VRK2 and EGFR) from "pan-VRK handles" (shared between VRK1 and VRK2, but different from EGFR).

**Why this priority**: Whether the drug hits only VRK1 or also VRK2 is a critical safety and pharmacology question that investors and scientific advisors will ask immediately after seeing the two-way data. The flag is already built; this story completes its analytical output.

**Independent Test**: Running the structural comparison with the VRK2 flag produces the three-way comparison file and an updated report. Fully testable without the biomarker or docking stories.

**Acceptance Scenarios**:

1. **Given** the VRK1 crystal structure is cached, **When** the user runs the structural comparison with the VRK2 flag, **Then** the output contains amino acid data for all three kinases at all 85 pocket positions.
2. **Given** no VRK2 crystal structure is available in the kinase database, **When** the command runs, **Then** an AlphaFold model is downloaded automatically and low-confidence residues are flagged clearly in the report.
3. **Given** the three-way table is generated, **Then** the report explicitly labels each selectivity position as "VRK1-specific," "pan-VRK," or "VRK-family vs EGFR."
4. **Given** the command is run without the VRK2 flag, **Then** the output is identical to the existing two-way report — no regression.

---

### User Story 3 — SCF-013 Computational Docking into 6AC9 (Priority: P3)

A researcher runs a docking command with a scaffold identifier and the VRK1 crystal structure, and receives a ranked list of predicted binding poses plus an assessment of which selectivity-candidate positions from the structural alignment are actually contacted by the scaffold. The output confirms or refutes the hypothesis that SCF-013 reaches the gatekeeper and hinge positions that differentiate VRK1 from EGFR.

**Why this priority**: Docking translates the structural selectivity analysis into actionable chemistry guidance — it answers whether SCF-013 physically fits the VRK1 back-pocket in a way EGFR cannot accommodate. This directs the next round of analog synthesis.

**Independent Test**: Validated by a positive-control run where the native ligand from the crystal structure is re-docked into its own binding site. The control must reproduce the known pose within 2 Å before SCF-013 results are considered credible.

**Acceptance Scenarios**:

1. **Given** the VRK1 crystal structure is cached and SCF-013 is in the compound library, **When** `pipeline dock --target VRK1 --scaffold SCF-013` is run, **Then** a docking report is produced with the top-ranked poses and their scores.
2. **Given** docking completes, **Then** the report maps each pose's atomic contacts onto the KLIFS selectivity candidates from the structural alignment.
3. **Given** the positive-control native ligand re-dock achieves RMSD ≤ 2.0 Å, **Then** the SCF-013 results are marked as coming from a validated setup; otherwise a warning is shown.
4. **Given** the scaffold identifier is not found in the compound library, **When** the command runs, **Then** an error names the missing scaffold and lists available scaffold IDs.

---

### User Story 4 — Co-Crystal Structure Preparation Guidance (Priority: P4)

A researcher runs a co-crystal planning command and receives a structured experimental brief: which crystal form of VRK1 to use, which soaking or co-crystallisation conditions to attempt first (based on published VRK1 crystal literature), which atoms in SCF-013 are most likely to need modification for crystallographic compatibility, and what resolution would be sufficient to resolve the key selectivity positions.

**Why this priority**: A co-crystal structure is the definitive validation experiment and the milestone typically required before Series A fundraising or partnership discussions. This story automates the literature synthesis that a structural biologist would otherwise spend weeks compiling.

**Independent Test**: Produces a standalone experimental brief document from publicly available crystallographic data, independent of any wet-lab execution.

**Acceptance Scenarios**:

1. **Given** the structural selectivity report exists, **When** `pipeline cocrystal --target VRK1 --scaffold SCF-013` is run, **Then** an experimental brief is written citing at least one published VRK1 crystallisation condition.
2. **Given** the brief is produced, **Then** it ranks the selectivity-candidate residues by the resolution required to resolve them and states the minimum target resolution for the experiment.
3. **Given** no published VRK1-specific soaking conditions are found, **Then** the brief defaults to conditions from the closest available structural homolog and notes the substitution explicitly.

---

### Edge Cases

- What if the CCLE mutation file format changes between DepMap releases — does the biomarker stage detect the schema change and fail with an actionable error?
- What if VRK2 has no kinase database entry and AlphaFold returns a low-confidence model covering the entire kinase domain — is the three-way comparison still produced with appropriate confidence caveats?
- What if SCF-013 cannot be parsed as a valid chemical structure — is the error message informative enough for a medicinal chemist to diagnose the problem?
- What if docking produces no poses above the scoring threshold — does the report explain why (e.g. steric clash, charge mismatch) rather than returning an empty file?
- What if the biomarker analysis finds many mutations with similar enrichment — how are ties broken and is the result stable across small sample fluctuations?

## Requirements *(mandatory)*

### Functional Requirements

**Biomarker Analysis (US1)**

- **FR-001**: The system MUST download and cache CCLE mutation data for the target lineage without requiring the researcher to handle files manually.
- **FR-002**: The system MUST compute enrichment of each somatic mutation in strongly-dependent versus non-dependent cell lines and rank results by statistical significance.
- **FR-003**: The biomarker report MUST include effect size (e.g. odds ratio or enrichment ratio) alongside p-value, so biological relevance can be judged independently of statistical significance.
- **FR-004**: The biomarker output MUST be injectable into the master research report using the same pattern as the existing depmap and structalign stages.
- **FR-005**: The lineage filter MUST accept any lineage name present in the DepMap dependency data, not only CNS/Brain.

**VRK2 Comparison (US2)**

- **FR-006**: The system MUST classify each of the 85 pocket positions as "VRK1-specific," "pan-VRK," or "VRK-family vs EGFR" based on the three-way amino acid comparison.
- **FR-007**: AlphaFold pLDDT filtering MUST be applied when a VRK2 crystal structure is unavailable, and the count of excluded low-confidence residues MUST appear in the report.
- **FR-008**: Running the structural comparison without the VRK2 flag MUST produce output identical to the existing two-way report — no regressions.

**Docking (US3)**

- **FR-009**: The system MUST accept a scaffold identifier from the existing compound library and resolve it to a chemical structure automatically, without requiring manual SMILES input.
- **FR-010**: Docking MUST use the binding site as defined by the 6 Å distance criterion already computed during the structural alignment stage.
- **FR-011**: The docking report MUST annotate each pose with which pocket selectivity candidates are within contact distance of the ligand.
- **FR-012**: A positive-control docking run using the native crystal ligand MUST be performed and its RMSD from the crystallographic pose reported before SCF-013 results are shown.
- **FR-013**: When AutoDock Vina is not installed, the command MUST fail immediately with a message naming the missing tool and providing the one-line installation command.

**Co-Crystal Brief (US4)**

- **FR-014**: The system MUST retrieve published crystallisation conditions for VRK1 from publicly accessible structural biology data sources.
- **FR-015**: The brief MUST identify which atoms in the scaffold are predicted to form the key selectivity contacts and flag any atoms likely to interfere with crystal packing.
- **FR-016**: The brief MUST specify the minimum resolution needed to resolve each key selectivity-candidate residue in the binding site.

### Key Entities

- **BiomarkerResult**: Co-mutation identifier, lineage, enrichment score, p-value, count of dependent lines with mutation, count of non-dependent lines with mutation.
- **ThreeWayComparison**: KLIFS position, subpocket, VRK1 amino acid, VRK2 amino acid, EGFR amino acid, VRK1–VRK2 difference type, VRK1–EGFR difference type, selectivity class label.
- **DockingPose**: Pose rank, docking score, RMSD to native ligand (positive-control run), list of contacted KLIFS positions within 4 Å.
- **CrystalBrief**: Crystal form reference, soaking or co-crystallisation conditions with literature source, minimum target resolution, flagged scaffold atoms, ranked contact residues.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A researcher produces the CNS/Brain biomarker report in a single command with no manual file handling, completing in under 10 minutes on a standard laptop.
- **SC-002**: The three-way VRK1/VRK2/EGFR comparison classifies 100% of the 85 pocket positions — no position is left without a selectivity label.
- **SC-003**: The docking positive-control (native ligand re-docked into 6AC9) achieves RMSD ≤ 2.0 Å, confirming the setup is reliable before SCF-013 results are interpreted.
- **SC-004**: The co-crystal brief cites at least one published VRK1 crystallisation condition and specifies the minimum resolution required to resolve each key selectivity-candidate position.
- **SC-005**: All four new pipeline commands follow the existing cache/force pattern — re-running any command without `--force` completes in under 5 seconds by reading cached outputs.

## Assumptions

- The DepMap VRK1 dependency cache from spec-003 is already present on disk; the biomarker stage reads from it rather than re-downloading gene effect data.
- CCLE somatic mutation data is publicly available from the Broad DepMap portal under the same download manifest used in spec-003.
- SCF-013 SMILES is already present in the compound library produced by the existing fetch stage and does not require manual entry.
- Docking uses AutoDock Vina (free, local, open-source). The pipeline checks for its presence at startup and prints installation instructions if it is missing, rather than crashing mid-run.
- The co-crystal brief is a synthesis of publicly accessible literature — it is a research planning document, not a wet-lab result.
- VRK1 6AC9 remains the reference structure for all four stories; no new PDB download is required beyond what spec-003 already produced.
- The VRK2 paralog comparison (US2) uses the `--include-vrk2` flag already wired into the CLI from spec-003; this story completes its analytical output, not its CLI plumbing.
