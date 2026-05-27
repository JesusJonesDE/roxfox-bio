from dataclasses import dataclass, field
from pathlib import Path

from pipeline.models import Target


TARGETS: dict[str, Target] = {
    "VRK1": Target(
        gene_name="VRK1",
        uniprot_id="Q99986",
        program_code="RXF-001",
        indications=["SMA", "Oncology"],
        display_name="VRK1 — Vaccinia Related Kinase 1",
        chembl_id="CHEMBL1293199",   # Serine/threonine-protein kinase VRK1
    ),
    "IGHMBP2": Target(
        gene_name="IGHMBP2",
        uniprot_id="P38935",
        program_code="RXF-002",
        indications=["SMARD1"],
        display_name="IGHMBP2 — Immunoglobulin Mu Binding Protein 2",
        chembl_id=None,   # Not present in ChEMBL — no published bioactivity data
    ),
    "VCP": Target(
        gene_name="VCP",
        uniprot_id="P55072",
        program_code="RXF-003",
        indications=["FTD"],
        display_name="VCP — Valosin-Containing Protein",
        chembl_id="CHEMBL1075145",   # Transitional endoplasmic reticulum ATPase
    ),
}

# Open Targets Ensembl IDs (resolved from prior API research)
ENSEMBL_IDS: dict[str, str] = {
    "VRK1": "ENSG00000100749",
    "IGHMBP2": "ENSG00000132471",
    "VCP": "ENSG00000197140",
}


@dataclass
class Settings:
    data_dir: Path = field(default_factory=lambda: Path("data"))
    cache_max_age_days: int = 30

    @property
    def cache_dir(self) -> Path:
        return self.data_dir / "cache"

    @property
    def results_dir(self) -> Path:
        return self.data_dir / "results"

    @property
    def shared_structures_dir(self) -> Path:
        return self.data_dir / "cache" / "shared" / "structures"

    @property
    def manifest_path(self) -> Path:
        return self.data_dir / "pipeline_manifest.json"
