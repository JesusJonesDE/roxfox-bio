import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings
from pipeline.models import Target

AF_BASE = "https://alphafold.ebi.ac.uk/api/prediction"


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_prediction(uniprot_id: str) -> list[dict]:
    with httpx.Client(timeout=20) as client:
        r = client.get(f"{AF_BASE}/{uniprot_id}")
        if r.status_code == 404:
            return []
        r.raise_for_status()
        return r.json()


def fetch_alphafold(target: Target, settings: Settings) -> tuple[dict, int]:
    models = _get_prediction(target.uniprot_id)
    if not models:
        return {"models": [], "note": "No AlphaFold model found"}, 0

    # Take the first (canonical) model
    m = models[0]
    result = {
        "entry_id": m.get("entryId"),
        "uniprot_id": target.uniprot_id,
        "model_url": m.get("cifUrl") or m.get("pdbUrl"),
        "pae_image_url": m.get("paeImageUrl"),
        "mean_plddt": m.get("meanPlddt") or _estimate_mean_plddt(m),
        "model_created_date": m.get("modelCreatedDate"),
        "latest_version": m.get("latestVersion"),
        "raw": m,
    }
    return result, 1


def _estimate_mean_plddt(model: dict) -> float | None:
    """Some AF responses include pLDDT in a different field."""
    confidence = model.get("confidenceScore") or model.get("plddt")
    return float(confidence) if confidence is not None else None
