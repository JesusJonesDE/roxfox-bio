# {{ target.gene_name }} — Research Dossier

**Program**: {{ target.program_code }} | **Generated**: {{ generated_date }} | **Pipeline version**: 0.1

---

## Overview

**Target**: {{ target.display_name }}
**UniProt**: {{ target.uniprot_id }}{% if target.chembl_id %} | **ChEMBL**: {{ target.chembl_id }}{% endif %}
**Indications**: {{ target.indications | join(", ") }}

{{ target.gene_name }} is being investigated as part of the RoxFox Bio drug discovery programme for {{ target.indications | join(" and ") }}.

---

## Genetic Evidence

{% if genetic.has_data %}
**Open Targets association score**: {{ genetic.top_score }} (scale 0–1; higher = stronger evidence)

**Top disease associations**:
{% for assoc in genetic.top_diseases %}
- **{{ assoc.name }}** — score: {{ "%.3f" | format(assoc.score) }}
{% endfor %}

**Tractability**:
{% for t in genetic.tractability %}
- {{ t.modality }}: {{ t.label }}
{% endfor %}
{% else %}
No genetic evidence data available from Open Targets. This may indicate the target has not yet been systematically evaluated in GWAS or functional genomics studies.
{% endif %}

---

## Bioactivity Summary

{% if compounds.has_data %}
- **Total ChEMBL bioactivity records**: {{ compounds.total_raw }}
- **Compounds with activity ≤ 10µM**: {{ compounds.total_filtered }}
- **Compounds passing Lipinski Ro5**: {{ compounds.ro5_pass }}
- **Potency range**: {{ compounds.best_nm }} nM – {{ compounds.worst_nm }} nM
- **Median potency**: {{ compounds.median_nm }} nM
{% else %}
No bioactivity data found in ChEMBL for this target. This may reflect that {{ target.gene_name }} is an early-stage or under-explored target with limited published chemical matter.
{% endif %}

---

## Scaffold Highlights

{% if scaffolds.has_data %}
Top scaffolds by number of active compounds (IC50/Ki/Kd ≤ 10µM):

| Rank | Scaffold ID | Compounds | Best potency (nM) | Median potency (nM) |
|------|-------------|-----------|-------------------|---------------------|
{% for s in scaffolds.top %}| {{ loop.index }} | {{ s.scaffold_id }} | {{ s.compound_count }} | {{ s.best_potency_nm }} | {{ s.median_potency_nm }} |
{% endfor %}

{{ scaffolds.total }} unique Murcko scaffolds identified across all active compounds.
{% else %}
No scaffold data available (requires active compounds in ChEMBL).
{% endif %}

---

## Structural Data

{% if structures.has_data %}
**Experimental structures available**: {{ structures.pdb_count }}
{% if structures.best_pdb %}
**Best resolution structure**: {{ structures.best_pdb.structure_id }} ({{ structures.best_pdb.resolution_angstrom }} Å, {{ structures.best_pdb.method }}){% if structures.best_pdb.has_ligand %} — ligand co-crystallised ✓{% endif %}

{% endif %}
{% if structures.ligand_bound_count > 0 %}
**Ligand-bound structures**: {{ structures.ligand_bound_count }} of {{ structures.pdb_count }} have a co-crystallised small molecule, confirming an accessible binding pocket.
{% else %}
No ligand-bound structures available — binding pocket confirmation will require computational prediction or crystallography campaigns.
{% endif %}

{% if structures.alphafold %}
**AlphaFold predicted model**: {{ structures.alphafold.structure_id }}{% if structures.alphafold.mean_plddt %} (mean pLDDT: {{ "%.1f" | format(structures.alphafold.mean_plddt) }}){% endif %}
{% endif %}
{% else %}
No structural data found in PDB or AlphaFold for this target. Structure-based approaches will require experimental crystallography or homology modelling.
{% endif %}

---

## Selectivity Profile

{% if selectivity.has_data %}
- **Compounds assessed for off-target activity**: {{ selectivity.compounds_assessed }}
- **Compounds with significant off-target activity** (> 3 unrelated targets at ≤ 1µM): {{ selectivity.flagged }}
{% if selectivity.flagged > 0 %}
Compounds flagged for off-target liability should be deprioritised or require selectivity optimisation before advancement.
{% else %}
No compounds showed concerning off-target activity, suggesting the current chemical matter may have reasonable selectivity profiles.
{% endif %}
{% else %}
Selectivity profiling not yet completed or no compounds available for assessment.
{% endif %}

---

## Competitive Landscape

{% if competitive.has_data %}
**Clinical trials identified** (query: {{ target.gene_name }}): {{ competitive.trial_count }}

{% if competitive.trials %}
| Status | Phase | Title | Sponsor |
|--------|-------|-------|---------|
{% for t in competitive.trials[:10] %}| {{ t.status or "—" }} | {{ t.phase or "—" }} | {{ t.title[:60] if t.title else "—" }}{% if t.title and t.title|length > 60 %}…{% endif %} | {{ t.sponsor or "—" }} |
{% endfor %}
{% endif %}

{% if competitive.approved_drugs %}
**Approved drugs targeting {{ target.gene_name }}**: {{ competitive.approved_drugs | join(", ") }}
{% else %}
No approved drugs directly targeting {{ target.gene_name }} identified — consistent with the pre-clinical status of this target.
{% endif %}
{% else %}
No competitive landscape data available.
{% endif %}

---

## Data Gaps & Limitations

{% if gaps %}
The following data gaps were identified during pipeline execution:

{% for gap in gaps %}
- {{ gap }}
{% endfor %}
{% else %}
All data sources returned results for this target. No significant data gaps identified.
{% endif %}

---

*Generated by RoxFox Bio Research Pipeline v0.1 — for internal use only. Data sourced from ChEMBL, Open Targets, RCSB PDB, AlphaFold EBI, and ClinicalTrials.gov (all public databases).*
