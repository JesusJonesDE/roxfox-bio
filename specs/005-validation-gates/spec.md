# Feature Specification: Computational Validation Gates

**Feature Branch**: `005-validation-gates`

**Created**: 2026-05-28

**Status**: Draft

## User Scenarios & Testing *(mandatory)*

### User Story 1 — ADMET Drug-likeness Gate (Priority: P1)

A researcher runs a scaffold through automated drug-likeness and absorption screening before investing in further computation. The system predicts BBB penetration, metabolic stability, CYP inhibition risk, aqueous solubility, and oral bioavailability from the compound structure alone — no wet lab required. A pass/fail decision is returned with a human-readable report and a numeric score per property.

**Why this priority**: BBB penetration is the single most critical filter for the CNS indication (VRK1 in brain tumors). A compound that cannot cross the blood-brain barrier cannot be a candidate regardless of binding affinity. This gate kills bad candidates earliest and cheapest.

**Independent Test**: Running `pipeline validate --target VRK1 --scaffold SCF-009 --gate admet` produces a report with pass/fail for each ADMET property and an overall gate decision within 60 seconds, without requiring wet-lab data.

**Acceptance Scenarios**:

1. **Given** a valid scaffold SMILES, **When** the ADMET gate runs, **Then** the system returns scores for BBB penetration, metabolic stability, CYP1A2/2D6/3A4 inhibition risk, aqueous solubility, and oral bioavailability
2. **Given** ADMET results, **When** any CNS-critical property fails its threshold (e.g. BBB penetration predicted negative), **Then** the gate is marked FAIL and the failing property is highlighted
3. **Given** a multi-fragment SMILES (salt form), **When** the gate runs, **Then** the largest fragment is used automatically without error
4. **Given** a scaffold that passes all thresholds, **When** the gate completes, **Then** it is marked PASS and the result is cached

---

### User Story 2 — MM-GBSA Binding Affinity Rescoring Gate (Priority: P2)

A researcher wants a more physically accurate binding affinity estimate than the raw docking score. The system takes the top docking pose, runs a physics-based free energy rescoring calculation entirely on the local machine, and returns a ΔG estimate in kcal/mol. This replaces the raw Vina score as the primary affinity ranking signal.

**Why this priority**: Raw docking scores have poor correlation with experimental IC50s. MM-GBSA is substantially more accurate and changes the scaffold ranking. It directly determines which scaffold gets prioritised for wet-lab work.

**Independent Test**: Running `pipeline validate --target VRK1 --scaffold SCF-009 --gate mmgbsa` produces a ΔG estimate and a PASS/FAIL based on a −7 kcal/mol threshold using only the locally cached PDB and docking pose.

**Acceptance Scenarios**:

1. **Given** a docked scaffold with at least one pose, **When** the MM-GBSA gate runs, **Then** it uses the top-ranked pose and returns ΔG in kcal/mol
2. **Given** MM-GBSA ΔG > −7.0 kcal/mol, **Then** the gate is marked FAIL
3. **Given** MM-GBSA ΔG ≤ −7.0 kcal/mol, **Then** the gate is marked PASS
4. **Given** no prior docking result exists, **When** the gate is requested, **Then** the system raises a clear error directing the user to run `pipeline dock` first

---

### User Story 3 — Selectivity Docking Panel Gate (Priority: P2)

A researcher needs to know if a scaffold is selective for the target kinase versus key off-targets before committing to synthesis. The system docks the scaffold into a panel of off-target structures (VRK2, EGFR, CDK2, PLK1) and computes a selectivity index: ratio of target affinity to best off-target affinity. A selectivity index ≥ 10× is a pass.

**Why this priority**: A scaffold that binds many kinases equally is not a drug. Selectivity must be assessed computationally before paying for wet-lab kinase panels ($5–15k per compound). This gate identifies whether a scaffold has inherent selectivity before spending money.

**Independent Test**: Running `pipeline validate --target VRK1 --scaffold SCF-009 --gate selectivity` downloads or reuses cached off-target structures, docks the scaffold into each, and returns a selectivity table and index within 2 hours locally.

**Acceptance Scenarios**:

1. **Given** a target scaffold and its top docking affinity, **When** the selectivity gate runs, **Then** it docks into VRK2, EGFR, CDK2, and PLK1 and returns affinity for each
2. **Given** selectivity index < 10×, **Then** gate is FAIL with a warning identifying the problematic off-target
3. **Given** selectivity index ≥ 10×, **Then** gate is PASS
4. **Given** a target that is one of the panel kinases, **Then** that kinase is excluded from the off-target panel automatically
5. **Given** an off-target PDB is not cached, **When** the gate runs, **Then** it downloads the best available structure automatically

---

### User Story 4 — MD Pose Stability Gate (Priority: P3)

A researcher wants to confirm the docked binding pose is physically stable and not a docking artefact. The system runs a short molecular dynamics simulation (50–100 ns) of the top docking pose in explicit solvent on the local machine, measures RMSD of the ligand over time, and reports whether the pose is stable or drifts out of the pocket.

**Why this priority**: MD is the highest-confidence computational filter but most expensive. It runs last, only on scaffolds that have passed gates 1–3. It prevents false positives from proceeding to wet lab.

**Independent Test**: Running `pipeline validate --target VRK1 --scaffold SCF-009 --gate md` on the M1 Max completes a 50 ns simulation in under 6 hours and returns a RMSD trajectory and a PASS/FAIL.

**Acceptance Scenarios**:

1. **Given** a scaffold passing gates 1–3, **When** the MD gate runs, **Then** it simulates 50 ns of the protein-ligand complex in explicit solvent
2. **Given** mean ligand RMSD > 3.0 Å over the last 25 ns, **Then** gate is FAIL
3. **Given** mean ligand RMSD ≤ 3.0 Å, **Then** gate is PASS
4. **Given** the simulation would exceed available memory, **Then** the system falls back to a 10 ns implicit solvent run with a warning
5. **Given** `pipeline validate` runs all gates, **Then** MD runs last and only if gates 1–3 all passed

---

### User Story 5 — Gate Dashboard and Wet-Lab Handoff Report (Priority: P2)

A researcher running multiple scaffolds needs a single view showing which candidates passed which gates. The system produces a gate dashboard (terminal table + markdown file) listing all scaffolds × gates with PASS/FAIL/PENDING status, and generates a wet-lab handoff report for any scaffold that passes all gates.

**Why this priority**: Without a consolidated view, researchers must manually check individual reports. The dashboard is the decision-making artefact that determines what goes to a CRO.

**Independent Test**: After validating multiple scaffolds, running `pipeline validate --target VRK1 --dashboard` prints a grid and writes `validation_dashboard.md` with one row per scaffold and one column per gate.

**Acceptance Scenarios**:

1. **Given** multiple scaffolds have been validated, **When** the dashboard command runs, **Then** it shows scaffolds as rows and gates (ADMET / MM-GBSA / Selectivity / MD) as columns with PASS / FAIL / PENDING / NOT RUN per cell
2. **Given** a scaffold passes all four gates, **Then** it is automatically included in a wet-lab handoff report containing docking report, ADMET summary, selectivity table, and MD RMSD summary
3. **Given** a scaffold fails any gate, **Then** the failing gate and reason are shown and the scaffold is excluded from the handoff report

---

### Edge Cases

- What happens when a scaffold has no docking result yet? → Gate errors clearly and directs to `pipeline dock` first; does not crash
- What if an off-target kinase has no published crystal structure? → AlphaFold2 structure is used with a visible warning flag in the report
- What if MM-GBSA fails to converge? → Gate is marked ERROR (not FAIL), numerical failure is reported, user can re-run with `--force`
- What if the MD simulation runs out of memory on 64 GB? → Falls back to implicit solvent 10 ns run with a warning in the report
- What if a scaffold SMILES cannot be parametrised for MD? → Gate fails immediately with a clear error; rest of pipeline unaffected
- What if a gate is re-run on a cached result? → Returns cached result immediately unless `--force` is passed

---

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST predict ADMET properties (BBB penetration, metabolic stability, CYP1A2/2D6/3A4 inhibition, aqueous solubility, oral bioavailability) from SMILES without wet-lab data
- **FR-002**: Each gate MUST produce a structured pass/fail decision, a numeric score, a human-readable reason, and a markdown report file
- **FR-003**: MM-GBSA rescoring MUST run entirely on local hardware, completing in under 2 hours per scaffold on Apple M1 Max 64 GB
- **FR-004**: The selectivity gate MUST dock the scaffold into at least VRK2, EGFR, CDK2, and PLK1 and compute a selectivity index vs. the primary target
- **FR-005**: The MD gate MUST simulate at least 50 ns of the protein-ligand complex and report mean ligand RMSD over the final half of the trajectory
- **FR-006**: The system MUST provide `pipeline validate --gate <name>` to run a single gate and `pipeline validate` (no gate flag) to run all gates in sequence
- **FR-007**: The system MUST provide `pipeline validate --dashboard` showing all scaffolds × gates as a pass/fail grid
- **FR-008**: The dashboard MUST auto-generate a wet-lab handoff report for any scaffold passing all gates
- **FR-009**: Gate results MUST be cached; re-running must return cached results unless `--force` is passed
- **FR-010**: All gates MUST run on local hardware; free public API calls for structure downloads only are permitted; no paid or account-required services
- **FR-011**: In sequential mode, the MD gate MUST only execute after gates 1–3 have passed
- **FR-012**: The system MUST automatically strip counterions from multi-fragment SMILES before any gate runs

### Key Entities

- **ValidationGate**: Single computational check (ADMET / MM-GBSA / Selectivity / MD) with status (PASS/FAIL/ERROR/PENDING), score, reason, and report path
- **ValidationResult**: All gate results for one scaffold × target pair; overall pass/fail; timestamp
- **GateDashboard**: Cross-scaffold grid aggregating ValidationResults for a target; drives wet-lab handoff decisions
- **WetLabHandoffReport**: Document for scaffolds passing all gates; contains all computational evidence for CRO submission
- **SelectivityPanel**: Configurable set of off-target structures used in the selectivity gate; default panel for kinases is VRK2/EGFR/CDK2/PLK1

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: ADMET gate completes in under 60 seconds per scaffold on local hardware
- **SC-002**: MM-GBSA rescoring completes in under 2 hours per scaffold on Apple M1 Max 64 GB
- **SC-003**: Selectivity panel docking (4 off-targets) completes in under 2 hours per scaffold locally
- **SC-004**: MD simulation (50 ns) completes in under 6 hours per scaffold on Apple M1 Max GPU
- **SC-005**: The gate dashboard correctly shows PASS/FAIL/PENDING for all validated scaffolds
- **SC-006**: At least one VRK1 scaffold (SCF-009 or SCF-156) produces a complete wet-lab handoff report after passing all gates
- **SC-007**: Running all gates except MD for the top 5 VRK1 scaffolds completes in under 5 hours on local hardware
- **SC-008**: Repeated runs on cached results return identical decisions (deterministic)

---

## Assumptions

- Docking results from `pipeline dock` must exist before any validation gate runs; this dependency is enforced
- The off-target selectivity panel (VRK2, EGFR, CDK2, PLK1) applies to kinase targets; VCP and IGHMBP2 panels are out of scope for this feature
- Apple M1 Max 64 GB unified memory is sufficient for 50 ns MD of a kinase-small molecule complex (~60,000 atoms explicit solvent); implicit solvent fallback is available
- ADMET pass thresholds: BBB penetration probability > 0.5; metabolic half-life > 30 min; CYP inhibition probability < 0.3 per isoform; solubility > 10 µg/mL; oral bioavailability > 30%
- MM-GBSA pass threshold: ΔG ≤ −7.0 kcal/mol
- Selectivity pass threshold: target affinity at least 10× better (more negative ΔG) than best off-target
- MD pass threshold: mean ligand heavy-atom RMSD ≤ 3.0 Å over final 25 ns
- Free locally-runnable tools are preferred; if unavailable, free public APIs (no account required) are the only permitted fallback
- IGHMBP2 is excluded from this feature (no compound library exists yet)
