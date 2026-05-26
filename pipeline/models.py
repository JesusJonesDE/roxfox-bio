from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Target:
    gene_name: str
    uniprot_id: str
    program_code: str
    indications: list[str]
    display_name: str
    chembl_id: Optional[str] = None


@dataclass
class BioactivityRecord:
    molecule_chembl_id: str
    canonical_smiles: Optional[str]
    standard_type: str        # IC50, Ki, Kd
    standard_value: float     # in reported units
    standard_units: str
    value_nm: float           # normalised to nM
    pchembl_value: Optional[float]
    target_chembl_id: str
    assay_chembl_id: str


@dataclass
class Compound:
    molecule_chembl_id: str
    smiles: str
    best_value_nm: float
    best_assay_type: str
    molecular_weight: float
    logp: float
    hbd: int
    hba: int
    rotatable_bonds: int
    ro5_violations: int
    passes_ro5: bool
    scaffold_id: Optional[str] = None
    off_target_flags: int = 0
    selectivity_flag: bool = False


@dataclass
class Scaffold:
    scaffold_smiles: str
    scaffold_id: str
    compound_count: int
    median_potency_nm: float
    best_potency_nm: float
    target_gene: str


@dataclass
class Structure:
    structure_id: str
    source: str               # "PDB" or "AlphaFold"
    method: str               # X-ray, Cryo-EM, NMR, Predicted
    has_ligand: bool
    resolution_angstrom: Optional[float] = None
    ligand_ids: list[str] = field(default_factory=list)
    chain_ids: list[str] = field(default_factory=list)
    mean_plddt: Optional[float] = None
    deposition_date: Optional[str] = None
    target_uniprot: str = ""


@dataclass
class CacheEntry:
    target_gene: str
    stage: str
    source: str
    fetched_at: str           # ISO 8601
    record_count: int
    file_path: str
    is_valid: bool = True
