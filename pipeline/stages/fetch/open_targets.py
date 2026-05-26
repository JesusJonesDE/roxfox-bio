import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from pipeline.config import Settings, ENSEMBL_IDS
from pipeline.models import Target

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

QUERY = """
query TargetInfo($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    biotype
    tractability {
      modality
      label
      value
    }
    associatedDiseases(page: {index: 0, size: 20}) {
      count
      rows {
        disease {
          id
          name
          therapeuticAreas {
            name
          }
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


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _graphql(ensembl_id: str) -> dict:
    with httpx.Client(timeout=30) as client:
        r = client.post(OT_GRAPHQL, json={"query": QUERY, "variables": {"ensemblId": ensembl_id}})
        r.raise_for_status()
        return r.json()


def fetch_open_targets(target: Target, settings: Settings) -> tuple[dict, int]:
    ensembl_id = ENSEMBL_IDS.get(target.gene_name)
    if not ensembl_id:
        return {"error": f"No Ensembl ID configured for {target.gene_name}"}, 0

    result = _graphql(ensembl_id)
    data = result.get("data", {}).get("target", {})
    if not data:
        return {"error": "No data returned from Open Targets"}, 0

    count = data.get("associatedDiseases", {}).get("count", 0)
    return data, count
