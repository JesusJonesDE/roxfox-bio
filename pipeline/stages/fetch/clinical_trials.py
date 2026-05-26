import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings
from pipeline.models import Target

CT_BASE = "https://clinicaltrials.gov/api/v2/studies"
PAGE_SIZE = 100


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_page(params: dict) -> dict:
    with httpx.Client(timeout=20) as client:
        r = client.get(CT_BASE, params=params)
        r.raise_for_status()
        return r.json()


def fetch_clinical_trials(target: Target, settings: Settings) -> tuple[dict, int]:
    studies = []
    next_token = None

    while True:
        params = {
            "query.term": target.gene_name,
            "pageSize": PAGE_SIZE,
            "format": "json",
            "fields": "NCTId,BriefTitle,OverallStatus,Phase,InterventionName,LeadSponsorName,StartDate",
        }
        if next_token:
            params["pageToken"] = next_token

        data = _get_page(params)
        page_studies = data.get("studies", [])
        for s in page_studies:
            proto = s.get("protocolSection", {})
            id_mod = proto.get("identificationModule", {})
            status_mod = proto.get("statusModule", {})
            design_mod = proto.get("designModule", {})
            sponsor_mod = proto.get("sponsorCollaboratorsModule", {})
            arms_mod = proto.get("armsInterventionsModule", {})

            interventions = [
                i.get("interventionName", "")
                for i in arms_mod.get("interventions", [])
            ]

            studies.append({
                "nct_id": id_mod.get("nctId"),
                "title": id_mod.get("briefTitle"),
                "status": status_mod.get("overallStatus"),
                "phase": design_mod.get("phases", [None])[0] if design_mod.get("phases") else None,
                "sponsor": sponsor_mod.get("leadSponsor", {}).get("leadSponsorName"),
                "interventions": interventions,
                "start_date": status_mod.get("startDateStruct", {}).get("date"),
            })

        next_token = data.get("nextPageToken")
        if not next_token or not page_studies:
            break

    return {"studies": studies, "gene_query": target.gene_name}, len(studies)
