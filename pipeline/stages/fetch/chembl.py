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
    # Use the component accession filter — exact UniProt match, not free-text search
    data = _get(f"{CHEMBL_BASE}/target", {
        "target_components__accession": uniprot_id,
        "target_type": "SINGLE PROTEIN",
        "format": "json",
    })
    targets = data.get("targets", [])
    if targets:
        return targets[0]["target_chembl_id"]
    raise ValueError(f"No ChEMBL SINGLE PROTEIN target found for UniProt {uniprot_id}")


def fetch_chembl(target: Target, settings: Settings) -> tuple[dict, int]:
    # Use hardcoded ChEMBL ID if available; avoids fragile API lookup
    if target.chembl_id is None:
        return {"note": f"{target.gene_name} has no ChEMBL target entry", "activities": []}, 0

    chembl_id = target.chembl_id or _resolve_chembl_target_id(target.uniprot_id)

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
