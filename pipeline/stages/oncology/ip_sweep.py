"""
IP / source sweep for lead compounds.

Looks up ChEMBL document sources for each compound in a given series
to determine whether they come from journal articles, patents, or both.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import TYPE_CHECKING

import httpx
import pandas as pd
from rich.console import Console
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings

if TYPE_CHECKING:
    from pipeline.cache import CacheManager

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_molecule_docs(mol_id: str) -> list[dict]:
    """Fetch all documents that report activity for a given molecule."""
    params = {
        "molecule_chembl_id": mol_id,
        "limit": 100,
        "format": "json",
    }
    with httpx.Client(timeout=20) as client:
        r = client.get(f"{CHEMBL_BASE}/activity", params=params,
                       headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"})
        r.raise_for_status()
    acts = r.json().get("activities", [])
    # Return unique documents for this molecule
    seen = set()
    docs = []
    for a in acts:
        doc_id = a.get("document_chembl_id")
        if doc_id and doc_id not in seen:
            seen.add(doc_id)
            docs.append({
                "document_chembl_id": doc_id,
                "document_journal": a.get("document_journal"),
                "document_year": a.get("document_year"),
                "src_id": a.get("src_id"),  # 1=literature, 2=deposited, 7=patent
            })
    return docs


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_document(doc_id: str) -> dict:
    with httpx.Client(timeout=15) as client:
        r = client.get(f"{CHEMBL_BASE}/document/{doc_id}",
                       params={"format": "json"},
                       headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"})
        r.raise_for_status()
    return r.json()


def run_ip_sweep(gene: str, series_ids: list[str], settings: Settings,
                 cache: "CacheManager", console: Console) -> Path:
    """
    Run IP sweep on compounds in specified scaffold series.
    series_ids: list of scaffold IDs to check (e.g. ['SCF-013', 'SCF-001'])
    """
    results_dir = settings.results_dir / gene
    compounds_path = results_dir / "compounds_filtered.csv"
    triage_path = results_dir / "scaffold_triage.csv"

    if not compounds_path.exists():
        raise FileNotFoundError(f"{compounds_path} not found")
    if not triage_path.exists():
        raise FileNotFoundError(f"{triage_path} not found — run triage first")

    compounds_df = pd.read_csv(compounds_path)
    triage_df = pd.read_csv(triage_path)

    # Get compound IDs from requested series
    target_compounds = compounds_df[compounds_df["scaffold_id"].isin(series_ids)][
        ["molecule_chembl_id", "scaffold_id", "best_value_nm", "passes_ro5", "selectivity_flag"]
    ].copy()

    console.print(f"  {gene}: IP sweep on {len(target_compounds)} compounds in {series_ids}")

    rows = []
    for _, comp in target_compounds.iterrows():
        mol_id = comp["molecule_chembl_id"]
        try:
            docs = _get_molecule_docs(mol_id)
            for doc in docs:
                # Enrich with full doc metadata for patents
                doc_type = "Unknown"
                title = ""
                doi = ""
                try:
                    full_doc = _get_document(doc["document_chembl_id"])
                    doc_type = full_doc.get("doc_type", "Unknown")
                    title = full_doc.get("title", "")
                    doi = full_doc.get("doi", "") or full_doc.get("patent_id", "")
                    time.sleep(0.1)
                except Exception:
                    pass

                rows.append({
                    "molecule_chembl_id": mol_id,
                    "scaffold_id": comp["scaffold_id"],
                    "best_value_nm": comp["best_value_nm"],
                    "passes_ro5": comp["passes_ro5"],
                    "selectivity_flag": comp["selectivity_flag"],
                    "document_chembl_id": doc["document_chembl_id"],
                    "doc_type": doc_type,
                    "journal": doc.get("document_journal", ""),
                    "year": doc.get("document_year"),
                    "title": title,
                    "doi_or_patent": doi,
                })
            time.sleep(0.15)
        except Exception as exc:
            rows.append({
                "molecule_chembl_id": mol_id,
                "scaffold_id": comp["scaffold_id"],
                "best_value_nm": comp["best_value_nm"],
                "passes_ro5": comp["passes_ro5"],
                "selectivity_flag": comp["selectivity_flag"],
                "document_chembl_id": "",
                "doc_type": "ERROR",
                "journal": str(exc),
                "year": None,
                "title": "",
                "doi_or_patent": "",
            })

    ip_df = pd.DataFrame(rows)
    out_csv = results_dir / "ip_sweep.csv"
    ip_df.to_csv(out_csv, index=False)

    _write_ip_report(gene, ip_df, series_ids, results_dir)

    patent_count = len(ip_df[ip_df["doc_type"].str.upper() == "PATENT"]) if not ip_df.empty else 0
    journal_count = len(ip_df[ip_df["doc_type"].str.upper() == "PUBLICATION"]) if not ip_df.empty else 0
    console.print(f"  {gene}: {journal_count} journal sources, {patent_count} patent sources found")

    return out_csv


def _write_ip_report(gene, df, series_ids, results_dir):
    lines = [
        f"# {gene} — IP Landscape Sweep",
        f"**Series analysed:** {', '.join(series_ids)}",
        "",
        "## Summary",
        "",
    ]

    if df.empty:
        lines.append("*No document data retrieved.*")
    else:
        patents = df[df["doc_type"].str.upper() == "PATENT"].drop_duplicates("document_chembl_id")
        journals = df[df["doc_type"].str.upper() == "PUBLICATION"].drop_duplicates("document_chembl_id")

        lines += [
            f"- Journal publications: **{len(journals)}**",
            f"- Patents: **{len(patents)}**",
            "",
        ]

        for series_id in series_ids:
            series_df = df[df["scaffold_id"] == series_id]
            compounds = series_df["molecule_chembl_id"].unique()
            lines += [f"### {series_id} ({len(compounds)} compounds)", ""]

            for mol_id in compounds:
                mol_docs = series_df[series_df["molecule_chembl_id"] == mol_id]
                nm = mol_docs.iloc[0]["best_value_nm"] if not mol_docs.empty else "?"
                lines.append(f"**{mol_id}** (best: {nm} nM)")
                for _, doc in mol_docs.iterrows():
                    if not doc["document_chembl_id"]:
                        continue
                    ref = doc.get("doi_or_patent") or doc["document_chembl_id"]
                    year = int(doc["year"]) if doc["year"] and not pd.isna(doc["year"]) else "?"
                    lines.append(f"  - [{doc['doc_type']}] {doc['title'][:80] if doc['title'] else doc['journal']} ({year}) — `{ref}`")
                lines.append("")

        lines += [
            "## IP Risk Assessment",
            "",
        ]

        if len(patents) == 0:
            lines += [
                "**No patents found** — compounds in this series appear to come from academic literature only.",
                "This suggests the chemical space may be free-to-operate, but a formal FTO (Freedom to Operate)",
                "analysis by a patent attorney is required before committing to this scaffold.",
            ]
        else:
            lines += [
                f"**{len(patents)} patent(s) identified** — this series is covered by patent literature.",
                "Review each patent for claim scope, filing dates, and expiry before proceeding.",
                "Consider whether claims cover the specific analogs needed for optimisation",
                "or whether design-around space exists.",
                "",
                "**Patents found:**",
            ]
            for _, p in patents.iterrows():
                year = int(p["year"]) if p["year"] and not pd.isna(p["year"]) else "?"
                lines.append(f"  - `{p['doi_or_patent'] or p['document_chembl_id']}` ({year}) — {p['title'][:100] if p['title'] else 'no title'}")

    out_md = results_dir / "ip_report.md"
    out_md.write_text("\n".join(lines) + "\n")
