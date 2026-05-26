import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings
from pipeline.models import Target

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _search(uniprot_id: str) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "entry",
        "request_options": {"results_verbosity": "minimal", "return_all_hits": True},
    }
    with httpx.Client(timeout=30) as client:
        r = client.post(RCSB_SEARCH, json=query)
        if r.status_code == 204:
            return []
        r.raise_for_status()
        data = r.json()
        return [hit["identifier"] for hit in data.get("result_set", [])]


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_entry(pdb_id: str) -> dict:
    with httpx.Client(timeout=15) as client:
        r = client.get(f"{RCSB_DATA}/{pdb_id}")
        r.raise_for_status()
        return r.json()


def fetch_pdb(target: Target, settings: Settings) -> tuple[dict, int]:
    pdb_ids = _search(target.uniprot_id)
    entries = []
    for pdb_id in pdb_ids[:100]:  # cap at 100 structures
        try:
            entry = _get_entry(pdb_id)
            entries.append(_extract_entry_meta(pdb_id, entry))
        except Exception:
            entries.append({"pdb_id": pdb_id, "error": "metadata fetch failed"})

    return {"structures": entries}, len(entries)


def _extract_entry_meta(pdb_id: str, entry: dict) -> dict:
    struct = entry.get("struct", {})
    refine = entry.get("refine", [{}])
    exp_method = entry.get("exptl", [{}])
    pdbx_vrptx = entry.get("pdbx_vrptx_summary_geometry", {})

    resolution = None
    if refine:
        resolution = refine[0].get("ls_d_res_high")

    method = "Unknown"
    if exp_method:
        method = exp_method[0].get("method", "Unknown")

    # Ligands: check nonpolymer entities
    nonpolymer = entry.get("rcsb_entry_info", {}).get("nonpolymer_entity_count", 0)
    ligand_ids: list[str] = []
    polymer_entities = entry.get("polymer_entities", [])

    chain_ids: list[str] = []
    for pe in (entry.get("polymer_entities") or []):
        for inst in (pe.get("polymer_entity_instances") or []):
            cid = inst.get("rcsb_polymer_entity_instance_container_identifiers", {}).get("auth_asym_id")
            if cid:
                chain_ids.append(cid)

    return {
        "pdb_id": pdb_id,
        "resolution_angstrom": resolution,
        "method": method,
        "has_ligand": nonpolymer > 0,
        "ligand_ids": ligand_ids,
        "chain_ids": chain_ids[:10],
        "deposition_date": entry.get("rcsb_accession_info", {}).get("deposit_date"),
    }
