from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Optional

import httpx
import pandas as pd
from rich.console import Console
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from pipeline.cache import CacheManager
from pipeline.config import Settings

RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
KLIFS_BASE = "https://klifs.net/api"
VRK2_UNIPROT = "O95551"
EGFR_PDB = "1M17"
EGFR_CHAIN = "A"
LIGAND_CUTOFF_DEFAULT = 6.0

# KLIFS position constants
GATEKEEPER_POS = 45
HINGE_POSITIONS = frozenset([46, 47, 48])
_KEY_SUBPOCKETS = frozenset(["P-loop", "Hinge", "Gatekeeper", "DFG"])

# Amino acid groups for difference_type classification
_ALIPHATIC = frozenset("GAVLIP")
_AROMATIC = frozenset("FYW")
_HYDROXYL = frozenset("ST")
_AMIDE = frozenset("NQ")
_CARBOXYLATE = frozenset("DE")
_BASIC = frozenset("RKH")
_HYDROPHOBIC = frozenset("MLIVF")
_POSITIVE = frozenset("RKH")
_NEGATIVE = frozenset("DE")
_HBOND_CAPABLE = frozenset("STNQYWHKRDCE")

_CONSERVATIVE_GROUPS = [
    _ALIPHATIC, _AROMATIC, _HYDROXYL, _AMIDE, _CARBOXYLATE, _BASIC, _HYDROPHOBIC,
]

# Solvent molecules to exclude from ligand detection
_SOLVENT = frozenset([
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "FMT",
    "ACT", "MPD", "BME", "DMS", "EOH", "IMD",
])

_PLDT_CONFIDENCE_THRESHOLD = 70.0

# KLIFS canonical 85-position subpocket assignments (1-indexed).
# Only the four subpockets in _KEY_SUBPOCKETS need to be precise.
_KLIFS_SUBPOCKET_MAP: dict[int, str] = {
    **{i: "P-loop" for i in range(1, 9)},         # I.1–I.8: glycine-rich P-loop
    **{i: "beta1-2" for i in range(9, 14)},        # β1–β2 strands
    **{i: "beta3-alphaC" for i in range(14, 25)},  # β3 + αC helix
    **{i: "Linker" for i in range(25, 45)},        # linker leading to gatekeeper
    45: "Gatekeeper",
    **{i: "Hinge" for i in (46, 47, 48)},
    **{i: "alphaD-alphaE" for i in range(49, 72)}, # αD-αE helices
    **{i: "DFG" for i in range(72, 77)},            # DFG motif
    **{i: "Other" for i in range(77, 86)},          # C-terminal extension
}


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _http_get(url: str, timeout: int = 60) -> httpx.Response:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(url)
        r.raise_for_status()
        return r


def _download_pdb(pdb_id: str, dest: Path) -> None:
    """Download a PDB file from RCSB to dest."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = _http_get(RCSB_DOWNLOAD.format(pdb_id=pdb_id))
    dest.write_bytes(r.content)


# ── PDB parsing helpers ────────────────────────────────────────────────────────

def _parse_resolution(pdb_text: str) -> Optional[float]:
    for line in pdb_text.splitlines():
        if "REMARK   2 RESOLUTION." in line:
            try:
                return float(line.split("RESOLUTION.")[1].strip().split()[0])
            except (IndexError, ValueError):
                pass
    return None


def _has_ligand(pdb_text: str) -> bool:
    for line in pdb_text.splitlines():
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname not in _SOLVENT:
                return True
    return False


def _select_best_vrk1_structure(gene_symbol: str, settings: Settings) -> tuple[str, Path]:
    """Scan the PDB cache and return (pdb_id, path) of the best structure."""
    struct_dir = settings.cache_dir / gene_symbol / "structures"
    if not struct_dir.exists() or not list(struct_dir.glob("*.pdb")):
        raise FileNotFoundError(
            f"No PDB structures found in {struct_dir}. "
            f"Run `pipeline fetch --target {gene_symbol}` first."
        )

    best: Optional[tuple[float, str, Path, bool]] = None  # (resolution, pdb_id, path, has_ligand)

    for pdb_path in struct_dir.glob("*.pdb"):
        text = pdb_path.read_text(errors="replace")
        resolution = _parse_resolution(text)
        ligand = _has_ligand(text)
        pdb_id = pdb_path.stem.upper()

        if resolution is None:
            continue

        if best is None:
            best = (resolution, pdb_id, pdb_path, ligand)
        elif ligand and not best[3]:
            best = (resolution, pdb_id, pdb_path, ligand)
        elif ligand == best[3] and resolution < best[0]:
            best = (resolution, pdb_id, pdb_path, ligand)

    if best is None:
        raise ValueError(f"No parseable PDB structures in {struct_dir}")

    return best[1], best[2]


# ── KLIFS pocket mapping ───────────────────────────────────────────────────────

def _klifs_pocket_opencadd(pdb_id: str) -> Optional[pd.DataFrame]:
    """Try opencadd KLIFS session for pocket residue mapping."""
    try:
        from opencadd.databases.klifs import setup_remote  # type: ignore
        session = setup_remote()
        structs = session.structures.by_structure_pdb_id(pdb_id)
        if structs is None or len(structs) == 0:
            return None
        klifs_id = structs.iloc[0]["structure.klifs_id"]
        pocket = session.pockets.by_structure_klifs_id(klifs_id)
        return pocket
    except Exception:
        return None


def _klifs_pocket_rest(pdb_id: str) -> Optional[pd.DataFrame]:
    """Fallback: call KLIFS REST API directly.

    The REST API returns:
      /structures_pdb_list  → [{structure_ID, pocket (85-char string), ...}]
      /interactions_match_residues?structure_ID=X → [{index, Xray_position, KLIFS_position}]

    'index' is the 1-based KLIFS pocket position (1–85).
    'Xray_position' is the PDB residue number ('_' for gaps).
    The 85-char 'pocket' string contains the amino acid letter at position index-1.
    """
    try:
        r = _http_get(f"{KLIFS_BASE}/structures_pdb_list?pdb-codes={pdb_id}")
        structs = r.json()
        if not structs:
            return None
        structure_id = structs[0]["structure_ID"]
        pocket_str = structs[0].get("pocket", "")

        r2 = _http_get(f"{KLIFS_BASE}/interactions_match_residues?structure_ID={structure_id}")
        residues = r2.json()
        if not residues:
            return None

        rows = []
        for entry in residues:
            idx = int(entry["index"])                  # 1-based KLIFS position
            xray = entry.get("Xray_position", "_")
            aa = pocket_str[idx - 1] if pocket_str and idx <= len(pocket_str) else "_"
            rows.append({
                "residue.klifs_id": idx,
                "residue.id": xray if xray != "_" else None,
                "residue.klifs_letter": aa,
                "subpocket.name": _KLIFS_SUBPOCKET_MAP.get(idx, "Other"),
            })
        return pd.DataFrame(rows) if rows else None
    except Exception:
        return None


def _get_klifs_pocket(pdb_id: str) -> pd.DataFrame:
    """Return KLIFS 85-position pocket DataFrame for a PDB structure."""
    pocket = _klifs_pocket_opencadd(pdb_id)
    if pocket is None:
        pocket = _klifs_pocket_rest(pdb_id)
    if pocket is None:
        raise ValueError(
            f"KLIFS mapping not available for {pdb_id}. "
            "This may indicate KLIFS does not have this structure or the API is unavailable. "
            "VRK1 is an atypical kinase — KLIFS coverage may be partial."
        )

    # Normalise column names from opencadd (REST path already produces clean names).
    # Check letter/amino FIRST — "id" appears inside "residue" so the ordering matters.
    col_map = {c: c for c in pocket.columns}
    for c in pocket.columns:
        cl = c.lower().replace(".", "_").replace(" ", "_")
        if "amino" in cl or ("letter" in cl and "klifs" in cl):
            col_map[c] = "residue.klifs_letter"
        elif "klifs" in cl and cl.endswith("_id") and "structure" not in cl:
            col_map[c] = "residue.klifs_id"
        elif cl in ("residue_id", "residue.id") or (
            "residue" in cl and cl.endswith("_id") and "klifs" not in cl
        ):
            col_map[c] = "residue.id"
        elif "subpocket" in cl or "region" in cl:
            col_map[c] = "subpocket.name"
    pocket = pocket.rename(columns=col_map)

    # Ensure required columns exist
    for req in ("residue.klifs_id", "residue.id", "residue.klifs_letter"):
        if req not in pocket.columns:
            raise ValueError(f"KLIFS pocket missing required column '{req}'. Got: {list(pocket.columns)}")

    if "subpocket.name" not in pocket.columns:
        pocket["subpocket.name"] = "Unknown"

    # Drop gap positions (amino acid "_" or "X") for residue ID mapping
    pocket = pocket.copy()
    pocket["is_gap"] = pocket["residue.klifs_letter"].isin(["_", "X", "-", "GAP", None])
    return pocket


# ── Binding site extraction (BioPython NeighborSearch) ────────────────────────

def _load_pdb_model(pdb_path: Path):
    from Bio.PDB import PDBParser  # type: ignore
    parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_path.stem, str(pdb_path))[0]


def _collect_ligand_atoms(model, exclude_altloc_non_A: bool = True) -> list:
    atoms = []
    for chain in model:
        for residue in chain:
            het, _, _ = residue.get_id()
            if het.startswith("H_") and residue.resname.strip() not in _SOLVENT:
                for atom in residue:
                    if exclude_altloc_non_A and atom.get_altloc() not in ("", "A"):
                        continue
                    atoms.append(atom)
    # If strict filter yields nothing, retry without altloc filter
    if not atoms and exclude_altloc_non_A:
        return _collect_ligand_atoms(model, exclude_altloc_non_A=False)
    return atoms


def _extract_binding_site(pdb_path: Path, cutoff: float) -> tuple[list, list]:
    """Return (ligand_atoms, binding_site_residues) using NeighborSearch."""
    from Bio.PDB import NeighborSearch  # type: ignore

    model = _load_pdb_model(pdb_path)
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)

    ligand_atoms = _collect_ligand_atoms(model)
    if not ligand_atoms:
        return [], []

    nearby: set = set()
    for lat in ligand_atoms:
        hits = ns.search(lat.get_vector().get_array(), cutoff, level="R")
        for res in hits:
            het, _, _ = res.get_id()
            if not het.startswith("H_") and het != "W":
                nearby.add(res)

    return ligand_atoms, list(nearby)


def _parse_resid(value) -> Optional[int]:
    """Extract integer residue number from KLIFS residue.id (may be 'I.1', '52A', etc.)."""
    if value is None or (isinstance(value, float) and value != value):
        return None
    s = str(value).strip()
    m = re.search(r"-?\d+", s)
    return int(m.group()) if m else None


def _find_residue(model, pdb_resid: int):
    """Find a residue by PDB sequence number across all chains."""
    for chain in model:
        for res in chain:
            het, seq, _ = res.get_id()
            if seq == pdb_resid and het == " ":
                return res
    return None


# ── Superimposer ───────────────────────────────────────────────────────────────

def _superimpose(
    vrk1_path: Path,
    egfr_path: Path,
    vrk1_klifs: pd.DataFrame,
    egfr_klifs: pd.DataFrame,
) -> tuple[float, int]:
    """Align EGFR onto VRK1 using shared KLIFS Cα atoms. Returns (rmsd, n_shared)."""
    from Bio.PDB import Superimposer  # type: ignore

    vrk1_model = _load_pdb_model(vrk1_path)
    egfr_model = _load_pdb_model(egfr_path)

    vrk1_pos = {
        int(row["residue.klifs_id"]): _parse_resid(row["residue.id"])
        for _, row in vrk1_klifs[~vrk1_klifs["is_gap"]].iterrows()
        if pd.notna(row["residue.id"])
    }
    vrk1_pos = {k: v for k, v in vrk1_pos.items() if v is not None}
    egfr_pos = {
        int(row["residue.klifs_id"]): _parse_resid(row["residue.id"])
        for _, row in egfr_klifs[~egfr_klifs["is_gap"]].iterrows()
        if pd.notna(row["residue.id"])
    }
    egfr_pos = {k: v for k, v in egfr_pos.items() if v is not None}

    shared = sorted(set(vrk1_pos.keys()) & set(egfr_pos.keys()))

    fixed_atoms, moving_atoms = [], []
    for pos in shared:
        vrk1_res = _find_residue(vrk1_model, vrk1_pos[pos])
        egfr_res = _find_residue(egfr_model, egfr_pos[pos])
        if vrk1_res and egfr_res:
            try:
                fixed_atoms.append(vrk1_res["CA"])
                moving_atoms.append(egfr_res["CA"])
            except KeyError:
                pass

    if len(fixed_atoms) < 3:
        return float("nan"), len(fixed_atoms)

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    return float(sup.rms), len(fixed_atoms)


# ── Difference type classification ────────────────────────────────────────────

def _classify_three_way(vrk1_aa: str, vrk2_aa: str, egfr_aa: str) -> str:
    """Classify selectivity of a KLIFS position across VRK1, VRK2, EGFR.

    Returns one of:
      "VRK1-specific"      — VRK1 differs from both VRK2 and EGFR
      "pan-VRK vs EGFR"   — VRK1 == VRK2, both differ from EGFR
      "VRK2 vs VRK1+EGFR" — VRK1 == EGFR, VRK2 differs
      "conserved"          — all three are equal
    GAP and LOW_CONF are treated as non-matching any other residue.
    """
    _non_residue = {"GAP", "_", "-", "X", "LOW_CONF", "N/A", None, ""}

    def _eq(a: str, b: str) -> bool:
        if a in _non_residue or b in _non_residue:
            return False
        return a == b

    v1_eq_v2 = _eq(vrk1_aa, vrk2_aa)
    v1_eq_eg = _eq(vrk1_aa, egfr_aa)

    if v1_eq_v2 and v1_eq_eg:
        return "conserved"
    if v1_eq_v2 and not v1_eq_eg:
        return "pan-VRK vs EGFR"
    if not v1_eq_v2 and v1_eq_eg:
        return "VRK2 vs VRK1+EGFR"
    return "VRK1-specific"


def _classify_difference(aa1: str, aa2: str) -> str:
    for aa in (aa1, aa2):
        if aa in ("GAP", "_", "-", "X", None, ""):
            return "gap"
    if aa1 == aa2:
        return "identical"
    if (aa1 in _POSITIVE and aa2 in _NEGATIVE) or (aa1 in _NEGATIVE and aa2 in _POSITIVE):
        return "charge"
    d1, d2 = aa1 in _HBOND_CAPABLE, aa2 in _HBOND_CAPABLE
    if d1 != d2:
        return "h_bond"
    for group in _CONSERVATIVE_GROUPS:
        if aa1 in group and aa2 in group:
            return "conservative"
    return "steric"


def _is_selectivity_candidate(difference_type: str, subpocket: str) -> bool:
    return difference_type not in ("identical", "gap") and subpocket in _KEY_SUBPOCKETS


# ── Comparison table ───────────────────────────────────────────────────────────

def _build_comparison(
    vrk1_klifs: pd.DataFrame,
    egfr_klifs: pd.DataFrame,
    vrk2_klifs: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Produce BindingSiteComparison table for all 85 KLIFS positions."""
    vrk1_map = {
        int(row["residue.klifs_id"]): (
            row["residue.klifs_letter"] if not row["is_gap"] else "GAP",
            str(row.get("subpocket.name", "Unknown")),
        )
        for _, row in vrk1_klifs.iterrows()
    }
    egfr_map = {
        int(row["residue.klifs_id"]): (
            row["residue.klifs_letter"] if not row["is_gap"] else "GAP",
            str(row.get("subpocket.name", "Unknown")),
        )
        for _, row in egfr_klifs.iterrows()
    }
    vrk2_map: dict[int, tuple[str, str]] = {}
    if vrk2_klifs is not None:
        for _, row in vrk2_klifs.iterrows():
            aa = "LOW_CONF" if row.get("low_confidence", False) else (
                row["residue.klifs_letter"] if not row["is_gap"] else "GAP"
            )
            vrk2_map[int(row["residue.klifs_id"])] = (aa, str(row.get("subpocket.name", "Unknown")))

    all_pos = sorted(set(vrk1_map.keys()) | set(egfr_map.keys()))
    rows = []
    for pos in all_pos:
        vrk1_aa, subpocket = vrk1_map.get(pos, ("GAP", "Unknown"))
        egfr_aa, _ = egfr_map.get(pos, ("GAP", "Unknown"))
        vrk2_aa = vrk2_map.get(pos, ("N/A", ""))[0] if vrk2_klifs is not None else "N/A"

        diff = _classify_difference(vrk1_aa, egfr_aa)
        vrk1_vrk2_diff = (
            _classify_difference(vrk1_aa, vrk2_aa)
            if vrk2_klifs is not None else None
        )
        selectivity_class = (
            _classify_three_way(vrk1_aa, vrk2_aa, egfr_aa)
            if vrk2_klifs is not None else None
        )
        row: dict = {
            "klifs_position": pos,
            "subpocket": subpocket,
            "vrk1_aa": vrk1_aa,
            "egfr_aa": egfr_aa,
            "identical_vrk1_egfr": diff == "identical",
            "difference_type": diff,
            "selectivity_candidate": _is_selectivity_candidate(diff, subpocket),
            "is_gatekeeper": pos == GATEKEEPER_POS,
            "is_hinge": pos in HINGE_POSITIONS,
            "notes": "",
        }
        if vrk2_klifs is not None:
            row["vrk2_aa"] = vrk2_aa
            row["vrk1_vrk2_diff"] = vrk1_vrk2_diff
            row["selectivity_class"] = selectivity_class
        rows.append(row)

    df = pd.DataFrame(rows)

    # Add gatekeeper notes
    gk = df[df["is_gatekeeper"]]
    if not gk.empty:
        idx = gk.index[0]
        vrk1_gk = gk.iloc[0]["vrk1_aa"]
        egfr_gk = gk.iloc[0]["egfr_aa"]
        df.at[idx, "notes"] = f"Gatekeeper: VRK1 {vrk1_gk} vs EGFR {egfr_gk}"

    return df


# ── VRK2 support ───────────────────────────────────────────────────────────────

def _fetch_vrk2_structure(settings: Settings, console: Console) -> Optional[tuple[str, Path]]:
    """Find or download a VRK2 structure. Returns (source_label, pdb_path) or None."""
    # Try KLIFS for a VRK2 crystal structure
    try:
        r = _http_get(f"{KLIFS_BASE}/kinases_list")
        kinases = r.json()
        vrk2_kinase = next(
            (k for k in kinases if k.get("kinase_name", "").upper() == "VRK2"), None
        )
        if vrk2_kinase:
            kid = vrk2_kinase.get("kinase_ID") or vrk2_kinase.get("kinaseID")
            r2 = _http_get(f"{KLIFS_BASE}/structures_list?kinase_ID={kid}")
            structs = r2.json()
            if structs:
                best = min(structs, key=lambda s: s.get("resolution", 99.0) or 99.0)
                pdb_id = best.get("pdb")
                if pdb_id:
                    dest = settings.shared_structures_dir / f"{pdb_id.upper()}.pdb"
                    if not dest.exists():
                        _download_pdb(pdb_id, dest)
                    return f"crystal/{pdb_id}", dest
    except Exception:
        pass

    # Fallback: AlphaFold
    console.print(
        f"  [dim]VRK2:[/dim] no crystal structure in KLIFS — downloading AlphaFold model..."
    )
    try:
        r = _http_get(ALPHAFOLD_API.format(uniprot_id=VRK2_UNIPROT))
        data = r.json()
        entries = data if isinstance(data, list) else [data]
        pdb_url = next(
            (e.get("pdbUrl") or e.get("pdb_url") for e in entries if e.get("pdbUrl") or e.get("pdb_url")),
            None,
        )
        if not pdb_url:
            return None
        dest = settings.shared_structures_dir / "VRK2_AF.pdb"
        if not dest.exists():
            r2 = _http_get(pdb_url, timeout=120)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(r2.content)
        return "alphafold/O95551", dest
    except Exception:
        return None


def _apply_plddt_filter(klifs_df: pd.DataFrame, pdb_path: Path) -> tuple[pd.DataFrame, int]:
    """For AlphaFold models: flag residues with pLDDT < threshold as low_confidence."""
    try:
        from Bio.PDB import PDBParser  # type: ignore
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("VRK2", str(pdb_path))[0]
        resid_to_plddt: dict[int, float] = {}
        for chain in struct:
            for res in chain:
                _, seq, _ = res.get_id()
                bfactors = [a.get_bfactor() for a in res.get_atoms()]
                if bfactors:
                    resid_to_plddt[seq] = min(bfactors)
    except Exception:
        return klifs_df, 0

    klifs_df = klifs_df.copy()
    low_conf = 0
    confidence_flags = []
    for _, row in klifs_df.iterrows():
        if row["is_gap"]:
            confidence_flags.append(False)
            continue
        resid = _parse_resid(row.get("residue.id"))
        if resid is not None:
            plddt = resid_to_plddt.get(resid, 100.0)
            is_low = plddt < _PLDT_CONFIDENCE_THRESHOLD
        else:
            is_low = False
        confidence_flags.append(is_low)
        if is_low:
            low_conf += 1
    klifs_df["low_confidence"] = confidence_flags
    return klifs_df, low_conf


# ── Report generation ─────────────────────────────────────────────────────────

def _write_report(
    gene_symbol: str,
    vrk1_pdb_id: str,
    vrk1_resolution: Optional[float],
    egfr_pdb_id: str,
    rmsd: float,
    n_shared: int,
    comparison: pd.DataFrame,
    vrk2_source: Optional[str],
    results_dir: Path,
) -> Path:
    gk = comparison[comparison["is_gatekeeper"]]
    gk_vrk1 = gk.iloc[0]["vrk1_aa"] if len(gk) > 0 else "?"
    gk_egfr = gk.iloc[0]["egfr_aa"] if len(gk) > 0 else "?"
    has_vrk2_cols = "vrk2_aa" in comparison.columns
    gk_vrk2 = (
        gk.iloc[0]["vrk2_aa"]
        if (has_vrk2_cols and len(gk) > 0 and gk.iloc[0]["vrk2_aa"] not in ("N/A", "GAP", "LOW_CONF", ""))
        else None
    )

    n_differ = int((comparison["difference_type"] != "identical").sum())
    candidates = comparison[comparison["selectivity_candidate"]]
    n_candidates = len(candidates)

    cand_rows = [
        f"| {int(r.klifs_position)} | {r.subpocket} | {r.vrk1_aa} | {r.egfr_aa} | {r.difference_type.replace('_',' ').title()} |"
        for _, r in candidates.iterrows()
    ]

    vrk2_section = ""
    if vrk2_source and has_vrk2_cols:
        vrk2_gk = gk_vrk2 or "?"
        _non = {"GAP", "LOW_CONF", "N/A", ""}
        sc = comparison["selectivity_class"]
        n_vrk1_specific = int((sc == "VRK1-specific").sum())
        n_pan_vrk = int((sc == "pan-VRK vs EGFR").sum())
        n_vrk2_specific = int((sc == "VRK2 vs VRK1+EGFR").sum())
        n_conserved = int((sc == "conserved").sum())

        vrk1_specific_cands = comparison[
            (comparison["selectivity_class"] == "VRK1-specific") &
            comparison["selectivity_candidate"]
        ]
        vrk1_cand_rows = [
            f"| {int(r.klifs_position)} | {r.subpocket} | {r.vrk1_aa} | {r.vrk2_aa} | {r.egfr_aa} |"
            for _, r in vrk1_specific_cands.iterrows()
        ]

        vrk2_section = f"""
---

## VRK2 Three-Way Comparison

**VRK2 source**: {vrk2_source}

### Selectivity Class Counts (all 85 KLIFS positions)

| Class | Count | Meaning |
|-------|-------|---------|
| VRK1-specific | {n_vrk1_specific} | VRK1 ≠ VRK2 and VRK1 ≠ EGFR — unique to VRK1 |
| pan-VRK vs EGFR | {n_pan_vrk} | VRK1 = VRK2 ≠ EGFR — shared VRK handle, avoids EGFR |
| VRK2 vs VRK1+EGFR | {n_vrk2_specific} | VRK2 ≠ VRK1 = EGFR — VRK2-specific position |
| conserved | {n_conserved} | VRK1 = VRK2 = EGFR — not useful for selectivity |

### Gatekeeper (KLIFS position 45)

VRK1 **{gk_vrk1}** | VRK2 **{vrk2_gk}** | EGFR **{gk_egfr}**

### VRK1-Specific Selectivity Candidates

Positions where VRK1 ≠ VRK2, VRK1 ≠ EGFR, **and** in a key subpocket — prime targets for VRK1-only inhibitors.

| KLIFS Position | Subpocket | VRK1 | VRK2 | EGFR |
|---------------|-----------|------|------|------|
{chr(10).join(vrk1_cand_rows) if vrk1_cand_rows else "| — | No VRK1-specific selectivity candidates | — | — | — |"}

A VRK1-selective inhibitor must exploit VRK1-specific positions (VRK1 ≠ VRK2 AND VRK1 ≠ EGFR).
Pan-VRK positions (VRK1 = VRK2 ≠ EGFR) are useful for avoiding EGFR while still hitting both VRK paralogs.
"""

    res_str = f"{vrk1_resolution:.2f} Å" if vrk1_resolution else "unknown resolution"
    rmsd_str = f"{rmsd:.2f} Å" if rmsd == rmsd else "N/A (insufficient shared positions)"

    report = f"""# {gene_symbol} Structural Selectivity Analysis — VRK1 vs. EGFR

**Generated**: {datetime.now().strftime("%Y-%m-%d")}
**Method**: KLIFS 85-position pocket alignment + BioPython Superimposer

---

## Structures Used

| Role | PDB ID | Resolution | Notes |
|------|--------|-----------|-------|
| {gene_symbol} | {vrk1_pdb_id} | {res_str} | Best available from pipeline PDB cache |
| EGFR reference | {egfr_pdb_id} | 2.60 Å | Wild-type EGFR, erlotinib complex (chain A) |

---

## Structural Alignment Quality

- Shared Cα atoms (KLIFS positions): **{n_shared}**
- Superimposition RMSD: **{rmsd_str}**

---

## Gatekeeper Residues (KLIFS Position 45)

| Kinase | Gatekeeper Residue |
|--------|--------------------|
| {gene_symbol} | **{gk_vrk1}** |
| EGFR | **{gk_egfr}** (T790M = resistance mutation) |

The gatekeeper is the primary selectivity determinant for ATP-competitive kinase inhibitors.
{"VRK1 and EGFR share the same gatekeeper — selectivity must come from other pocket positions." if gk_vrk1 == gk_egfr else f"VRK1 ({gk_vrk1}) and EGFR ({gk_egfr}) differ at the gatekeeper — this is the primary selectivity handle for SCF-013."}

---

## Binding Site Comparison

Total positions compared: **85** | Positions that differ: **{n_differ}** | Selectivity candidates: **{n_candidates}**

### Selectivity Candidates (differ AND in key subpocket)

| KLIFS Position | Subpocket | VRK1 | EGFR | Difference Type |
|---------------|-----------|------|------|----------------|
{chr(10).join(cand_rows) if cand_rows else "| — | No selectivity candidates identified | — | — | — |"}

---

## Selectivity Hypothesis

{"Based on the binding site comparison, VRK1 and EGFR differ at " + str(n_differ) + " positions. " + ("The gatekeeper difference (VRK1 " + gk_vrk1 + " vs EGFR " + gk_egfr + ") is likely the dominant selectivity driver: " + ("the larger VRK1 gatekeeper limits back-pocket access for certain inhibitor scaffolds, while EGFR's smaller gatekeeper accommodates them. SCF-013 likely exploits this asymmetry." if len(gk_vrk1) > len(gk_egfr) else "EGFR's larger gatekeeper restricts certain scaffolds that VRK1 accommodates. SCF-013 may prefer the more open VRK1 back-pocket.")) if gk_vrk1 != gk_egfr else "Selectivity over EGFR must be driven by positions outside the gatekeeper. The " + str(n_candidates) + " selectivity candidates in key subpockets (P-loop, Hinge, DFG) should be prioritised for structure-guided analog design."}
{vrk2_section}
"""
    path = results_dir / "structural_selectivity_report.md"
    path.write_text(report)
    return path


def update_research_report(
    gene_symbol: str,
    gk_vrk1: str,
    gk_egfr: str,
    n_candidates: int,
    rmsd: float,
    results_dir: Path,
) -> bool:
    """Inject structural findings into data/results/research_report.md. Returns True if updated."""
    report_path = results_dir.parent / "research_report.md"
    if not report_path.exists():
        return False

    text = report_path.read_text()
    if "Structural Selectivity (VRK1 vs EGFR)" in text:
        return False

    rmsd_str = f"{rmsd:.2f} Å" if rmsd == rmsd else "N/A"
    section = (
        f"\n\n### Structural Selectivity (VRK1 vs EGFR)\n\n"
        f"- Gatekeeper: VRK1 **{gk_vrk1}** | EGFR **{gk_egfr}**\n"
        f"- Alignment RMSD: {rmsd_str}\n"
        f"- Selectivity candidates in key subpockets: **{n_candidates}**\n\n"
        f"*Full report: [{gene_symbol}/structural_selectivity_report.md]"
        f"({gene_symbol}/structural_selectivity_report.md)*\n"
    )

    marker = f"## {gene_symbol}"
    if marker in text:
        idx = text.index(marker)
        end_of_line = text.find("\n", idx)
        text = text[: end_of_line + 1] + section + text[end_of_line + 1 :]
        report_path.write_text(text)
        return True

    report_path.write_text(text + section)
    return True


# ── Main entry point ──────────────────────────────────────────────────────────

def run_structalign(
    gene_symbol: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    include_vrk2: bool,
    cutoff: float,
    console: Console,
) -> None:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    # Cache hit check: all output files present and not forced
    if not force:
        comparison_path = results_dir / "binding_site_comparison.csv"
        report_path = results_dir / "structural_selectivity_report.md"
        if comparison_path.exists() and report_path.exists():
            console.print(
                f"  [dim]{gene_symbol:10}[/dim] structalign          "
                f"[yellow]SKIP[/yellow]  (output files present; use --force to re-run)"
            )
            return

    # Select VRK1 structure
    vrk1_pdb_id, vrk1_path = _select_best_vrk1_structure(gene_symbol, settings)
    vrk1_resolution = _parse_resolution(vrk1_path.read_text(errors="replace"))
    ligand_present = _has_ligand(vrk1_path.read_text(errors="replace"))

    # Fetch EGFR reference
    egfr_path = settings.shared_structures_dir / f"{EGFR_PDB}.pdb"
    if not egfr_path.exists():
        console.print(f"  [dim]{gene_symbol}:[/dim] downloading EGFR reference {EGFR_PDB}...")
        _download_pdb(EGFR_PDB, egfr_path)

    console.print(
        f"  [dim]{gene_symbol}:[/dim] best structure {vrk1_pdb_id} "
        f"({'ligand-bound' if ligand_present else 'apo'}"
        + (f", {vrk1_resolution:.2f} Å" if vrk1_resolution else "")
        + ")"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] EGFR reference {EGFR_PDB} (erlotinib, chain {EGFR_CHAIN})"
    )

    # KLIFS pocket mapping
    console.print(f"  [dim]{gene_symbol}:[/dim] mapping KLIFS binding site pockets...")
    vrk1_klifs = _get_klifs_pocket(vrk1_pdb_id)
    egfr_klifs = _get_klifs_pocket(EGFR_PDB)

    # NeighborSearch binding site extraction
    console.print(f"  [dim]{gene_symbol}:[/dim] extracting binding sites ({cutoff} Å cutoff)...")
    vrk1_lig, vrk1_bs = _extract_binding_site(vrk1_path, cutoff)
    egfr_lig, egfr_bs = _extract_binding_site(egfr_path, cutoff)

    n_vrk1_bs = len(vrk1_bs)
    n_egfr_bs = len(egfr_bs)

    # Superimposer alignment
    console.print(f"  [dim]{gene_symbol}:[/dim] aligning structures via shared KLIFS Cα atoms...")
    rmsd, n_shared = _superimpose(vrk1_path, egfr_path, vrk1_klifs, egfr_klifs)

    # VRK2 (optional, US3)
    vrk2_klifs: Optional[pd.DataFrame] = None
    vrk2_source: Optional[str] = None
    if include_vrk2:
        console.print(f"  [dim]{gene_symbol}:[/dim] fetching VRK2 structure...")
        result = _fetch_vrk2_structure(settings, console)
        if result:
            vrk2_source, vrk2_path = result
            try:
                vrk2_pdb_id = vrk2_path.stem.upper()
                vrk2_klifs_raw = _get_klifs_pocket(vrk2_pdb_id)
                if "alphafold" in vrk2_source:
                    vrk2_klifs, n_low_conf = _apply_plddt_filter(vrk2_klifs_raw, vrk2_path)
                    if n_low_conf > 0:
                        console.print(
                            f"  [yellow]VRK2:[/yellow] {n_low_conf} binding site residues "
                            f"excluded (pLDDT < {_PLDT_CONFIDENCE_THRESHOLD})"
                        )
                else:
                    vrk2_klifs = vrk2_klifs_raw
            except Exception as exc:
                console.print(f"  [yellow]VRK2 KLIFS mapping failed:[/yellow] {exc}")
        else:
            console.print(f"  [yellow]{gene_symbol}:[/yellow] VRK2 structure unavailable — skipping")

    # Build comparison table
    comparison = _build_comparison(vrk1_klifs, egfr_klifs, vrk2_klifs)

    # Write CSV outputs
    vrk1_klifs.to_csv(results_dir / "binding_site_vrk1.csv", index=False)
    egfr_klifs.to_csv(results_dir / "binding_site_egfr.csv", index=False)
    comparison.to_csv(results_dir / "binding_site_comparison.csv", index=False)

    # Write report
    _write_report(
        gene_symbol,
        vrk1_pdb_id,
        vrk1_resolution,
        EGFR_PDB,
        rmsd,
        n_shared,
        comparison,
        vrk2_source,
        results_dir,
    )

    # Inject into master research report (US4)
    gk = comparison[comparison["is_gatekeeper"]]
    gk_vrk1 = gk.iloc[0]["vrk1_aa"] if len(gk) > 0 else "?"
    gk_egfr = gk.iloc[0]["egfr_aa"] if len(gk) > 0 else "?"
    n_candidates = int(comparison["selectivity_candidate"].sum())
    updated = update_research_report(gene_symbol, gk_vrk1, gk_egfr, n_candidates, rmsd, results_dir)

    # Console summary
    rmsd_str = f"{rmsd:.2f} Å" if rmsd == rmsd else "N/A"
    n_differ = int((comparison["difference_type"] != "identical").sum())

    console.print(
        f"  [dim]{gene_symbol}:[/dim] binding site — {gene_symbol}: {n_vrk1_bs} residues | "
        f"EGFR: {n_egfr_bs} residues"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] KLIFS alignment RMSD: {rmsd_str} ({n_shared} shared Cα)"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] gatekeeper — {gene_symbol}: {gk_vrk1} | EGFR: {gk_egfr}"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] {n_differ} positions differ | {n_candidates} selectivity candidates"
    )
    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]structural_selectivity_report.md written[/green]"
    )
    if updated:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [dim]research_report.md updated[/dim]"
        )
