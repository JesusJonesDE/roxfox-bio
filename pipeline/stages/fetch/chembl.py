import time
from typing import Any

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings
from pipeline.models import Target

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAGE_SIZE = 1000
REQUEST_DELAY = 0.3  # seconds between pages


@retry(
    stop=stop_after_attempt(4),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get(url: str, params: dict) -> dict:
    with httpx.Client(timeout=30) as client:
        r = client.get(url, params=params, headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"})
        r.raise_for_status()
        return r.json()


def _resolve_chembl_target_id(uniprot_id: str) -> str:
    data = _get(f"{CHEMBL_BASE}/target", {"search": uniprot_id, "format": "json"})
    targets = data.get("targets", [])
    for t in targets:
        # Prefer SINGLE PROTEIN entries
        if t.get("target_type") == "SINGLE PROTEIN":
            for comp in t.get("target_components", []):
                for xref in comp.get("target_component_xrefs", []):
                    if xref.get("xref_src_db") == "UniProt" and xref.get("xref_id") == uniprot_id:
                        return t["target_chembl_id"]
    # Fallback: first result
    if targets:
        return targets[0]["target_chembl_id"]
    raise ValueError(f"No ChEMBL target found for UniProt {uniprot_id}")


def fetch_chembl(target: Target, settings: Settings) -> tuple[dict, int]:
    chembl_id = _resolve_chembl_target_id(target.uniprot_id)
    target.chembl_id = chembl_id

    records = []
    offset = 0
    while True:
        params = {
            "target_chembl_id": chembl_id,
            "assay_type__in": "B,F",
            "activity_type__in": "IC50,Ki,Kd",
            "limit": PAGE_SIZE,
            "offset": offset,
            "format": "json",
        }
        data = _get(f"{CHEMBL_BASE}/activity", params)
        page = data.get("activities", [])
        records.extend(page)

        total = data.get("page_meta", {}).get("total_count", 0)
        offset += PAGE_SIZE
        if offset >= total or not page:
            break
        time.sleep(REQUEST_DELAY)

    return {"chembl_target_id": chembl_id, "activities": records}, len(records)
