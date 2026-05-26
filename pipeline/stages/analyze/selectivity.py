from __future__ import annotations

import time
from typing import TYPE_CHECKING

import httpx
import pandas as pd
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

if TYPE_CHECKING:
    from rich.console import Console
    from pipeline.cache import CacheManager
    from pipeline.config import Settings

from pipeline.models import Compound

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
OFF_TARGET_THRESHOLD_NM = 1_000.0  # 1µM
OFF_TARGET_FLAG_COUNT = 3


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_activities(mol_id: str, target_chembl_id: str) -> list[dict]:
    params = {
        "molecule_chembl_id": mol_id,
        "standard_value__lte": OFF_TARGET_THRESHOLD_NM,
        "standard_units": "nM",
        "activity_type__in": "IC50,Ki,Kd",
        "limit": 200,
        "format": "json",
    }
    with httpx.Client(timeout=20) as client:
        r = client.get(f"{CHEMBL_BASE}/activity", params=params, headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"})
        r.raise_for_status()
        data = r.json()
    activities = data.get("activities", [])
    # Exclude the primary target itself
    return [a for a in activities if a.get("target_chembl_id") != target_chembl_id]


def run_selectivity_analysis(
    gene: str,
    compounds: list[Compound],
    settings: "Settings",
    cache: "CacheManager",
    console: "Console",
) -> list[Compound]:
    if not compounds:
        return compounds

    # Get the primary target ChEMBL ID from ChEMBL cache
    chembl_cache = cache.load(gene, "chembl")
    primary_target_id = chembl_cache.get("chembl_target_id") if chembl_cache else None

    # Query only for compounds passing Ro5 to keep request count manageable
    ro5_compounds = [c for c in compounds if c.passes_ro5]
    flagged = 0

    for c in ro5_compounds:
        try:
            off_targets = _get_activities(c.molecule_chembl_id, primary_target_id or "")
            unique_targets = len({a.get("target_chembl_id") for a in off_targets})
            c.off_target_flags = unique_targets
            c.selectivity_flag = unique_targets > OFF_TARGET_FLAG_COUNT
            if c.selectivity_flag:
                flagged += 1
            time.sleep(0.1)
        except Exception:
            pass  # Non-fatal; leave flags at defaults

    # Update CSV
    results_dir = settings.results_dir / gene
    csv_path = results_dir / "compounds_filtered.csv"
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        flags = {c.molecule_chembl_id: (c.off_target_flags, c.selectivity_flag) for c in compounds}
        df["off_target_flags"] = df["molecule_chembl_id"].map(lambda x: flags.get(x, (0, False))[0])
        df["selectivity_flag"] = df["molecule_chembl_id"].map(lambda x: flags.get(x, (0, False))[1])
        df.to_csv(csv_path, index=False)

    console.print(
        f"  [dim]{gene:10}[/dim] analyze  selectivity    [green]OK[/green]    "
        f"({len(ro5_compounds)} compounds profiled, {flagged} flagged)"
    )
    return compounds
