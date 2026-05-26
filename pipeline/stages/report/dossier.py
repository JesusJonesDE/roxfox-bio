from __future__ import annotations

import statistics
from datetime import date
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from jinja2 import Environment, FileSystemLoader

if TYPE_CHECKING:
    from pipeline.cache import CacheManager
    from pipeline.config import Settings

TEMPLATE_DIR = Path(__file__).parent
TEMPLATE_NAME = "dossier_template.md"


def generate_dossier(gene: str, settings: "Settings", cache: "CacheManager") -> Path:
    from pipeline.config import TARGETS

    target = TARGETS[gene]
    results_dir = settings.results_dir / gene
    results_dir.mkdir(parents=True, exist_ok=True)

    ctx = {
        "target": target,
        "generated_date": date.today().isoformat(),
        "genetic": _build_genetic_ctx(gene, cache),
        "compounds": _build_compounds_ctx(results_dir),
        "scaffolds": _build_scaffolds_ctx(results_dir),
        "structures": _build_structures_ctx(results_dir),
        "selectivity": _build_selectivity_ctx(results_dir),
        "competitive": _build_competitive_ctx(gene, cache),
        "gaps": [],
    }

    # Collect data gaps
    gaps = []
    if not ctx["genetic"]["has_data"]:
        gaps.append("Open Targets genetic evidence: no data returned")
    if not ctx["compounds"]["has_data"]:
        gaps.append(f"ChEMBL bioactivity: no compounds with activity ≤ 10µM found for {gene}")
    if not ctx["structures"]["has_data"]:
        gaps.append("PDB/AlphaFold: no structural data found")
    if ctx["competitive"]["trial_count"] == 0:
        gaps.append("ClinicalTrials.gov: no trials found matching target gene name")
    ctx["gaps"] = gaps

    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), autoescape=False)
    template = env.get_template(TEMPLATE_NAME)
    content = template.render(**ctx)

    out_path = results_dir / "dossier.md"
    out_path.write_text(content)
    return out_path


# ── Context builders ───────────────────────────────────────────────────────────

def _build_genetic_ctx(gene: str, cache: "CacheManager") -> dict:
    data = cache.load(gene, "open_targets")
    if not data or "error" in data:
        return {"has_data": False}

    rows = data.get("associatedDiseases", {}).get("rows", [])
    tractability = data.get("tractability", []) or []

    top_score = max((r.get("score", 0) for r in rows), default=0)
    top_diseases = sorted(rows, key=lambda r: r.get("score", 0), reverse=True)[:5]

    return {
        "has_data": bool(rows),
        "top_score": round(top_score, 3),
        "top_diseases": [
            {"name": r["disease"]["name"], "score": r.get("score", 0)}
            for r in top_diseases
        ],
        "tractability": [
            {"modality": t.get("modality", ""), "label": t.get("label", ""), "value": t.get("value")}
            for t in tractability
            if t.get("value")
        ],
    }


def _build_compounds_ctx(results_dir: Path) -> dict:
    csv = results_dir / "compounds_filtered.csv"
    if not csv.exists():
        return {"has_data": False}

    df = pd.read_csv(csv)
    if df.empty:
        return {"has_data": False}

    potencies = df["best_value_nm"].dropna().tolist()
    ro5_pass = int(df["passes_ro5"].sum()) if "passes_ro5" in df.columns else len(df)

    return {
        "has_data": True,
        "total_raw": "see ChEMBL cache",
        "total_filtered": len(df),
        "ro5_pass": ro5_pass,
        "best_nm": round(min(potencies), 1) if potencies else "—",
        "worst_nm": round(max(potencies), 1) if potencies else "—",
        "median_nm": round(statistics.median(potencies), 1) if potencies else "—",
    }


def _build_scaffolds_ctx(results_dir: Path) -> dict:
    csv = results_dir / "scaffolds.csv"
    if not csv.exists():
        return {"has_data": False}

    df = pd.read_csv(csv)
    if df.empty:
        return {"has_data": False}

    top = df.sort_values("compound_count", ascending=False).head(5)
    return {
        "has_data": True,
        "total": len(df),
        "top": top.to_dict("records"),
    }


def _build_structures_ctx(results_dir: Path) -> dict:
    csv = results_dir / "structures.csv"
    if not csv.exists():
        return {"has_data": False}

    df = pd.read_csv(csv)
    if df.empty:
        return {"has_data": False}

    pdb_df = df[df["source"] == "PDB"]
    af_df = df[df["source"] == "AlphaFold"]

    best_pdb = None
    if not pdb_df.empty:
        sorted_pdb = pdb_df.dropna(subset=["resolution_angstrom"]).sort_values("resolution_angstrom")
        if not sorted_pdb.empty:
            best_pdb = sorted_pdb.iloc[0].to_dict()

    af_row = af_df.iloc[0].to_dict() if not af_df.empty else None

    return {
        "has_data": len(df) > 0,
        "pdb_count": len(pdb_df),
        "ligand_bound_count": int(pdb_df["has_ligand"].sum()) if "has_ligand" in pdb_df.columns else 0,
        "best_pdb": best_pdb,
        "alphafold": af_row,
    }


def _build_selectivity_ctx(results_dir: Path) -> dict:
    csv = results_dir / "compounds_filtered.csv"
    if not csv.exists():
        return {"has_data": False}

    df = pd.read_csv(csv)
    if "selectivity_flag" not in df.columns:
        return {"has_data": False}

    assessed = int(df["off_target_flags"].notna().sum())
    flagged = int(df["selectivity_flag"].sum())

    return {
        "has_data": assessed > 0,
        "compounds_assessed": assessed,
        "flagged": flagged,
    }


def _build_competitive_ctx(gene: str, cache: "CacheManager") -> dict:
    data = cache.load(gene, "clinical_trials")
    if not data:
        return {"has_data": False, "trial_count": 0}

    studies = data.get("studies", [])
    return {
        "has_data": True,
        "trial_count": len(studies),
        "trials": studies,
        "approved_drugs": [],  # Placeholder — can be enriched from ChEMBL drug data
    }
