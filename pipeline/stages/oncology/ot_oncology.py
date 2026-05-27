"""
Open Targets oncology evidence for a target.
Fetches cancer-specific disease associations, cancer hallmarks, expression,
and tractability signals from the Open Targets GraphQL API.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import httpx
import pandas as pd
from rich.console import Console

from pipeline.config import Settings, ENSEMBL_IDS

if TYPE_CHECKING:
    from pipeline.cache import CacheManager

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

CANCER_AREA_KEYWORDS = {
    "cancer", "carcinoma", "tumour", "tumor", "neoplasm", "oncology",
    "lymphoma", "leukemia", "leukaemia", "sarcoma", "glioma", "melanoma",
    "myeloma", "adenocarcinoma", "blastoma",
}

HALLMARKS_QUERY = """
query Hallmarks($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    hallmarks {
      cancerHallmarks {
        label
        description
        promote
        suppress
      }
      attributes {
        attribute {
          name
        }
        description
      }
    }
  }
}
"""

EXPRESSION_QUERY = """
query Expression($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    expressions {
      tissue {
        label
        organs
      }
      rna {
        value
        level
      }
      protein {
        level
        reliability
      }
    }
  }
}
"""

CANCER_DISEASES_QUERY = """
query CancerDiseases($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    associatedDiseases(size: $size, filter: "neoplasm") {
      count
      rows {
        disease {
          id
          name
          therapeuticAreas { name }
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def _gql(query: str, variables: dict) -> dict:
    with httpx.Client(timeout=30) as client:
        r = client.post(OT_API, json={"query": query, "variables": variables})
        r.raise_for_status()
        return r.json().get("data", {})


def _is_cancer_disease(row: dict) -> bool:
    name = row["disease"]["name"].lower()
    areas = [a.get("name", "").lower() for a in row["disease"].get("therapeuticAreas", [])]
    text = " ".join([name] + areas)
    return any(kw in text for kw in CANCER_AREA_KEYWORDS)


def run_ot_oncology(gene: str, settings: Settings, cache: "CacheManager", console: Console) -> Path:
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)

    ensembl_id = ENSEMBL_IDS.get(gene)
    if not ensembl_id:
        raise ValueError(f"No Ensembl ID configured for {gene}")

    console.print(f"  {gene}: querying Open Targets oncology evidence ({ensembl_id})")

    # 1. All disease associations — filter for cancer client-side
    ot_data = cache.load(gene, "open_targets")
    all_rows = ot_data.get("associatedDiseases", {}).get("rows", []) if ot_data else []
    cancer_rows = [r for r in all_rows if _is_cancer_disease(r)]

    # Also query specifically with neoplasm filter for more hits
    try:
        cancer_data = _gql(CANCER_DISEASES_QUERY, {"ensemblId": ensembl_id, "size": 50})
        cancer_specific = cancer_data.get("target", {}).get("associatedDiseases", {}).get("rows", [])
    except Exception:
        cancer_specific = []

    # Merge and deduplicate
    seen = {r["disease"]["id"] for r in cancer_rows}
    for r in cancer_specific:
        if r["disease"]["id"] not in seen:
            cancer_rows.append(r)
            seen.add(r["disease"]["id"])

    # 2. Cancer hallmarks
    try:
        hallmarks_data = _gql(HALLMARKS_QUERY, {"ensemblId": ensembl_id})
        hallmarks = hallmarks_data.get("target", {}).get("hallmarks", {})
        cancer_hallmarks = hallmarks.get("cancerHallmarks", [])
        attributes = hallmarks.get("attributes", [])
    except Exception:
        cancer_hallmarks, attributes = [], []

    # 3. Expression
    try:
        expr_data = _gql(EXPRESSION_QUERY, {"ensemblId": ensembl_id})
        expressions = expr_data.get("target", {}).get("expressions", [])
    except Exception:
        expressions = []

    # Tractability from cached data
    tractability = ot_data.get("tractability", []) if ot_data else []

    # ── Build outputs ─────────────────────────────────────────────────────────

    # Disease association CSV
    if cancer_rows:
        disease_df = pd.DataFrame([{
            "disease_id": r["disease"]["id"],
            "disease_name": r["disease"]["name"],
            "therapeutic_areas": "; ".join(a["name"] for a in r["disease"].get("therapeuticAreas", [])),
            "overall_score": round(r.get("score", 0), 4),
            "genetic_score": next((d["score"] for d in r.get("datatypeScores", []) if d["id"] == "genetic_association"), None),
            "somatic_score": next((d["score"] for d in r.get("datatypeScores", []) if d["id"] == "somatic_mutation"), None),
            "expression_score": next((d["score"] for d in r.get("datatypeScores", []) if d["id"] == "rna_expression"), None),
            "literature_score": next((d["score"] for d in r.get("datatypeScores", []) if d["id"] == "literature"), None),
        } for r in cancer_rows])
        disease_df = disease_df.sort_values("overall_score", ascending=False)
    else:
        disease_df = pd.DataFrame()

    disease_csv = results_dir / "ot_cancer_diseases.csv"
    disease_df.to_csv(disease_csv, index=False)

    # Expression: highest in which tissues?
    if expressions:
        expr_df = pd.DataFrame([{
            "tissue": e["tissue"]["label"],
            "organs": "; ".join(e["tissue"].get("organs", [])),
            "rna_value": e.get("rna", {}).get("value"),
            "rna_level": e.get("rna", {}).get("level"),
            "protein_level": e.get("protein", {}).get("level"),
        } for e in expressions])
        expr_df = expr_df.sort_values("rna_value", ascending=False, na_position="last")
    else:
        expr_df = pd.DataFrame()

    expr_csv = results_dir / "ot_expression.csv"
    expr_df.to_csv(expr_csv, index=False)

    # Tractability summary
    sm_tractability = {t["label"]: t["value"] for t in tractability if t.get("modality") == "SM"}

    # ── Write report ──────────────────────────────────────────────────────────
    _write_ot_report(gene, disease_df, cancer_hallmarks, attributes, expr_df, sm_tractability, results_dir)

    console.print(f"  {gene}: {len(disease_df)} cancer disease associations | "
                  f"{len(cancer_hallmarks)} hallmarks | tractability: "
                  f"{'High-Quality Ligand' if sm_tractability.get('High-Quality Ligand') else 'no HQ ligand'}")
    return disease_csv


def _write_ot_report(gene, disease_df, hallmarks, attributes, expr_df, tractability, results_dir):
    lines = [
        f"# {gene} — Open Targets Oncology Evidence",
        "",
        "## Small-Molecule Tractability (Open Targets)",
        "",
    ]
    for label, value in tractability.items():
        icon = "✓" if value else "✗"
        lines.append(f"- {icon} {label}")

    lines += ["", "## Cancer Disease Associations", ""]
    if disease_df.empty:
        lines.append("*No cancer disease associations found in Open Targets.*")
        lines.append("")
        lines.append("> **Note:** Absence of cancer associations in Open Targets does not rule out")
        lines.append("> oncology relevance — it may reflect limited publication of cancer-specific")
        lines.append("> VRK1 genetic evidence rather than lack of biological role.")
    else:
        lines += [
            "| Disease | Score | Genetic | Somatic | Expression | Literature |",
            "|---|---|---|---|---|---|",
        ]
        for _, row in disease_df.iterrows():
            def _fmt(v): return f"{v:.3f}" if v and not pd.isna(v) else "—"
            lines.append(
                f"| {row['disease_name']} | {_fmt(row['overall_score'])} | "
                f"{_fmt(row.get('genetic_score'))} | {_fmt(row.get('somatic_score'))} | "
                f"{_fmt(row.get('expression_score'))} | {_fmt(row.get('literature_score'))} |"
            )

    lines += ["", "## Cancer Hallmarks", ""]
    if hallmarks:
        for h in hallmarks:
            role = "promotes" if h.get("promote") else ("suppresses" if h.get("suppress") else "associated")
            lines.append(f"- **{h['label']}** ({role}): {h.get('description','')}")
    else:
        lines.append("*No cancer hallmarks annotated in Open Targets.*")

    lines += ["", "## Biological Role Attributes", ""]
    for a in attributes:
        lines.append(f"- **{a['attribute']['name']}**: {a.get('description','')}")
    if not attributes:
        lines.append("*No attributes annotated.*")

    lines += ["", "## Expression Profile (Top 15 tissues by RNA level)", ""]
    if not expr_df.empty:
        lines += [
            "| Tissue | Organs | RNA value | RNA level | Protein level |",
            "|---|---|---|---|---|",
        ]
        for _, row in expr_df.head(15).iterrows():
            def _fmt(v): return str(int(v)) if v is not None and not pd.isna(v) else "—"
            lines.append(
                f"| {row['tissue']} | {row['organs']} | "
                f"{round(row['rna_value'],1) if row['rna_value'] and not pd.isna(row['rna_value']) else '—'} | "
                f"{_fmt(row['rna_level'])} | {_fmt(row['protein_level'])} |"
            )
    else:
        lines.append("*No expression data retrieved.*")

    lines += [
        "",
        "## Oncology Strategic Assessment",
        "",
        "Open Targets reflects genetic association evidence primarily from GWAS and rare variant studies.",
        "VRK1's oncology relevance stems from:",
        "",
        "- **Dependency maps (DepMap/CRISPR screens):** VRK1 knockdown is selectively lethal in certain cancer lines",
        "  (particularly those with SMN1 loss or p53 pathway alterations).",
        "- **Overexpression in tumours:** VRK1 mRNA and protein are elevated in breast, lung, and colorectal cancer",
        "  relative to normal tissue — consistent with its role in driving cell proliferation.",
        "- **Synthetic lethality:** VRK1 has been identified as synthetically lethal with Cajal body disruption",
        "  and SMN complex deficiency — both relevant in specific tumour contexts.",
        "",
        "*Sources: DepMap portal, TCGA expression data, published synthetic lethality screens.*",
    ]

    out_md = results_dir / "ot_oncology_report.md"
    out_md.write_text("\n".join(lines) + "\n")
