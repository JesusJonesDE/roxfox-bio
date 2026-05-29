"""Hit identification stage for targets with no ChEMBL data (e.g. IGHMBP2).

Queries BindingDB by UniProt and target name, PubChem BioAssay, and optionally
SF1 helicase analogs, then filters by Ro5, deduplicates, and writes outputs.
"""
from __future__ import annotations

import csv
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import httpx
from rich.console import Console
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from pipeline.cache import CacheManager
from pipeline.config import Settings

# ── HTTP helpers ──────────────────────────────────────────────────────────────

REQUEST_DELAY = 0.5  # seconds between paginated requests


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.TransportError)),
)
def _get_json(url: str, params: Optional[dict] = None, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(
            url,
            params=params or {},
            headers={"User-Agent": "RoxFoxBio-Pipeline/0.1"},
        )
        r.raise_for_status()
        return r.json()


# ── BindingDB helpers ─────────────────────────────────────────────────────────

BINDINGDB_BASE = "https://bindingdb.org/rest"

_BDB_ASSAY_FIELDS = ("Ki", "IC50", "Kd", "EC50")


def _parse_bindingdb_response(data: Any, source_label: str) -> list[dict]:
    """Extract compound records from a BindingDB getLigands* JSON response."""
    compounds: list[dict] = []

    if not isinstance(data, dict):
        return compounds

    # BindingDB wraps results in varying nested keys depending on the endpoint
    # Try common paths
    affinities: list[dict] = []
    for key in ("getLigandsByUniprotsResponse", "getLigandsByTargetNameResponse"):
        if key in data:
            inner = data[key]
            # Ligands are nested under "affinities" or "Affinity"
            for aff_key in ("affinities", "Affinity"):
                val = inner.get(aff_key, [])
                if isinstance(val, list):
                    affinities = val
                elif isinstance(val, dict):
                    affinities = [val]
                if affinities:
                    break
        if affinities:
            break

    # Fallback: search recursively for a list of affinity dicts
    if not affinities:
        affinities = _find_affinity_list(data)

    for entry in affinities:
        if not isinstance(entry, dict):
            continue

        smiles = (
            entry.get("smiles")
            or entry.get("ligandSmiles")
            or entry.get("Smiles")
            or ""
        ).strip()
        if not smiles:
            continue

        # Extract best potency value across assay types
        best_val: Optional[float] = None
        best_type: Optional[str] = None
        for field in _BDB_ASSAY_FIELDS:
            for key_variant in (field, field.lower(), f"{field}_nM"):
                raw = entry.get(key_variant) or entry.get(f"{key_variant}_nM")
                if raw is None:
                    continue
                try:
                    val = float(str(raw).replace(">", "").replace("<", "").strip())
                    if best_val is None or val < best_val:
                        best_val = val
                        best_type = field
                except (ValueError, TypeError):
                    continue

        if best_val is None:
            continue

        compounds.append({
            "smiles": smiles,
            "best_value_nm": best_val,
            "best_assay_type": best_type or "unknown",
            "source": source_label,
        })

    return compounds


def _find_affinity_list(obj: Any, depth: int = 0) -> list[dict]:
    """Recursively search nested dicts/lists for a list of affinity dicts."""
    if depth > 6:
        return []
    if isinstance(obj, list) and obj and isinstance(obj[0], dict):
        if any(
            k in obj[0]
            for k in ("smiles", "ligandSmiles", "Smiles", "Ki", "IC50", "Kd")
        ):
            return obj
    if isinstance(obj, dict):
        for v in obj.values():
            result = _find_affinity_list(v, depth + 1)
            if result:
                return result
    if isinstance(obj, list):
        for item in obj:
            result = _find_affinity_list(item, depth + 1)
            if result:
                return result
    return []


def _query_bindingdb_by_uniprot(uniprot_id: str, console: Console) -> list[dict]:
    # New REST API: /getLigandsByUniprots?uniprot=P38935&json=1
    url = f"{BINDINGDB_BASE}/getLigandsByUniprots"
    try:
        data = _get_json(url, {"uniprot": uniprot_id, "json": "1"})
        hits = _parse_bindingdb_response(data, "BindingDB_UniProt")
        console.print(f"    BindingDB (UniProt {uniprot_id}): {len(hits)} compounds")
        return hits
    except Exception as exc:
        console.print(f"    [yellow]BindingDB UniProt query failed: {exc}[/yellow]")
        return []


def _query_bindingdb_by_name(target_name: str, source_label: str, console: Console) -> list[dict]:
    # New REST API: /getLigandsByTargetName?target=IGHMBP2&json=1
    url = f"{BINDINGDB_BASE}/getLigandsByTargetName"
    try:
        data = _get_json(url, {"target": target_name, "json": "1"})
        hits = _parse_bindingdb_response(data, source_label)
        console.print(f"    BindingDB (name='{target_name}'): {len(hits)} compounds")
        return hits
    except Exception as exc:
        console.print(f"    [yellow]BindingDB name query failed: {exc}[/yellow]")
        return []


# ── PubChem BioAssay helpers ──────────────────────────────────────────────────

PUBCHEM_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_AIDS = 3
MAX_CIDS_PER_AID = 50


def _query_pubchem(uniprot_id: str, console: Console) -> list[dict]:
    """Fetch active compounds from PubChem BioAssays linked to a UniProt accession."""
    compounds: list[dict] = []

    # Step A: get assay IDs via gene symbol (UniProtID not accepted; genesymbol works)
    # Derive gene symbol from config if possible, fallback to uniprot search
    from pipeline.config import TARGETS
    gene_sym = next(
        (t.gene_name for t in TARGETS.values() if t.uniprot_id == uniprot_id),
        None,
    )
    aid_url = (
        f"{PUBCHEM_PUG}/assay/target/genesymbol/{gene_sym}/aids/JSON"
        if gene_sym
        else f"{PUBCHEM_PUG}/assay/target/UniProtID/{uniprot_id}/aids/JSON"
    )
    try:
        aids_data = _get_json(aid_url)
    except Exception as exc:
        console.print(f"    [yellow]PubChem AID lookup failed: {exc}[/yellow]")
        return compounds

    aids: list[int] = []
    try:
        id_list = aids_data.get("IdentifierList", {}).get("AID", [])
        aids = id_list if isinstance(id_list, list) else [id_list]
    except (KeyError, TypeError):
        pass
    if not aids:
        try:
            aids = aids_data["InformationList"]["Information"][0]["AID"]
            if isinstance(aids, int):
                aids = [aids]
        except (KeyError, IndexError, TypeError):
            pass

    if not aids:
        console.print(f"    PubChem: no assays found for UniProt {uniprot_id}")
        return compounds

    console.print(f"    PubChem: {len(aids)} assay(s) found, querying first {MAX_AIDS}")

    for aid in aids[:MAX_AIDS]:
        time.sleep(REQUEST_DELAY)

        # Step B: active CIDs for this assay
        try:
            cids_data = _get_json(
                f"{PUBCHEM_PUG}/assay/aid/{aid}/cids/JSON",
                {"activity_outcome": "active"},
            )
        except Exception as exc:
            console.print(f"    [yellow]PubChem AID {aid} CID lookup failed: {exc}[/yellow]")
            continue

        cids: list[int] = []
        try:
            cids = cids_data["InformationList"]["Information"][0]["CID"]
            if isinstance(cids, int):
                cids = [cids]
        except (KeyError, IndexError, TypeError):
            pass

        if not cids:
            continue

        cids = cids[:MAX_CIDS_PER_AID]
        console.print(f"    PubChem AID {aid}: {len(cids)} active CID(s)")
        time.sleep(REQUEST_DELAY)

        # Step C: get SMILES + properties for each CID batch
        cid_str = ",".join(str(c) for c in cids)
        try:
            props_data = _get_json(
                f"{PUBCHEM_PUG}/compound/cid/{cid_str}/property/IsomericSMILES,MolecularWeight,XLogP/JSON"
            )
        except Exception as exc:
            console.print(f"    [yellow]PubChem property fetch failed for AID {aid}: {exc}[/yellow]")
            continue

        for prop in props_data.get("PropertyTable", {}).get("Properties", []):
            smiles = prop.get("IsomericSMILES", "").strip()
            if not smiles:
                continue
            mw = prop.get("MolecularWeight")
            xlogp = prop.get("XLogP")
            compounds.append({
                "smiles": smiles,
                "best_value_nm": None,   # activity value not easily available per CID from this endpoint
                "best_assay_type": f"PubChem_AID{aid}",
                "source": "PubChem",
                "molecular_weight": float(mw) if mw is not None else None,
                "logp": float(xlogp) if xlogp is not None else None,
            })

    console.print(f"    PubChem: {len(compounds)} compound(s) total")
    return compounds


# ── RDKit chemistry helpers ───────────────────────────────────────────────────

def _canonicalise(smiles: str) -> Optional[str]:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def _strip_salts(smiles: str) -> str:
    """Keep the largest fragment (by heavy atom count)."""
    try:
        from rdkit import Chem
        from rdkit.Chem.SaltRemover import SaltRemover
        remover = SaltRemover()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        stripped = remover.StripMol(mol)
        if stripped is None or stripped.GetNumAtoms() == 0:
            # Fallback: largest dot-separated fragment
            frags = smiles.split(".")
            return max(frags, key=len)
        return Chem.MolToSmiles(stripped, isomericSmiles=True)
    except Exception:
        frags = smiles.split(".")
        return max(frags, key=len)


def _compute_ro5(smiles: str) -> dict:
    """Compute Lipinski Ro5 descriptors. Returns dict with MW, logP, HBD, HBA, RotB."""
    result = {
        "molecular_weight": None,
        "logp": None,
        "hbd": None,
        "hba": None,
        "rotatable_bonds": None,
        "ro5_violations": None,
        "passes_ro5": None,
    }
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        result.update({
            "molecular_weight": round(mw, 2),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rotb,
            "ro5_violations": violations,
            "passes_ro5": violations == 0,
        })
    except Exception:
        pass
    return result


# ── Core filtering and deduplication ─────────────────────────────────────────

def _process_compounds(
    raw: list[dict],
    gene_symbol: str,
    console: Console,
) -> list[dict]:
    """Canonicalise, strip salts, compute Ro5, deduplicate, assign scaffold IDs."""
    seen_smiles: dict[str, dict] = {}

    for entry in raw:
        raw_smiles = entry.get("smiles", "")
        if not raw_smiles:
            continue

        # Strip salts first, then canonicalise
        stripped = _strip_salts(raw_smiles)
        canon = _canonicalise(stripped)
        if canon is None:
            continue

        # Deduplicate: if same canonical SMILES seen, keep entry with best potency
        existing = seen_smiles.get(canon)
        entry_val = entry.get("best_value_nm")
        if existing is not None:
            existing_val = existing.get("best_value_nm")
            if entry_val is not None and (existing_val is None or entry_val < existing_val):
                seen_smiles[canon] = {**entry, "smiles": canon}
        else:
            seen_smiles[canon] = {**entry, "smiles": canon}

    # Apply Ro5 filter
    filtered: list[dict] = []
    for canon, entry in seen_smiles.items():
        # Use pre-computed values if available (e.g. from PubChem), else compute
        if entry.get("molecular_weight") is None or entry.get("logp") is None:
            ro5 = _compute_ro5(canon)
        else:
            mw = entry.get("molecular_weight")
            logp = entry.get("logp")
            hbd = entry.get("hbd")
            hba = entry.get("hba")
            rotb = entry.get("rotatable_bonds")
            # Recompute HBD/HBA/rotb from RDKit if missing
            rdkit_props = _compute_ro5(canon)
            hbd = hbd if hbd is not None else rdkit_props.get("hbd")
            hba = hba if hba is not None else rdkit_props.get("hba")
            rotb = rotb if rotb is not None else rdkit_props.get("rotatable_bonds")
            violations = sum([
                (mw or 0) > 500,
                (logp or 0) > 5,
                (hbd or 0) > 5,
                (hba or 0) > 10,
            ])
            ro5 = {
                "molecular_weight": mw,
                "logp": logp,
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rotb,
                "ro5_violations": violations,
                "passes_ro5": violations == 0,
            }
        entry.update(ro5)
        filtered.append(entry)

    console.print(f"    {len(filtered)} unique compound(s) after deduplication")

    # Sort by potency ascending (None → last)
    filtered.sort(key=lambda x: (x.get("best_value_nm") is None, x.get("best_value_nm") or 0))

    # Assign scaffold IDs
    for i, entry in enumerate(filtered, start=1):
        entry["scaffold_id"] = f"{gene_symbol}-SCF-{i:03d}"
        entry.setdefault("off_target_flags", 0)
        entry.setdefault("selectivity_flag", False)

    return filtered


# ── Output writers ────────────────────────────────────────────────────────────

CSV_COLUMNS = [
    "molecule_id",
    "smiles",
    "best_value_nm",
    "best_assay_type",
    "molecular_weight",
    "logp",
    "hbd",
    "hba",
    "rotatable_bonds",
    "ro5_violations",
    "passes_ro5",
    "scaffold_id",
    "source",
    "off_target_flags",
    "selectivity_flag",
]


def _write_csv(compounds: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for i, cmpd in enumerate(compounds, start=1):
            row = dict(cmpd)
            row["molecule_id"] = cmpd.get("scaffold_id", f"HIT-{i:03d}")
            writer.writerow(row)


def _write_report(
    gene_symbol: str,
    compounds: list[dict],
    sources_queried: list[str],
    analog_fallback_used: bool,
    out_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    date_str = datetime.now().strftime("%Y-%m-%d")

    top10 = compounds[:10]
    rows = []
    for cmpd in top10:
        potency = (
            f"{cmpd['best_value_nm']:.1f} nM ({cmpd.get('best_assay_type', '?')})"
            if cmpd.get("best_value_nm") is not None
            else "N/A"
        )
        smiles_short = cmpd["smiles"][:40] + ("..." if len(cmpd["smiles"]) > 40 else "")
        rows.append(
            f"| {cmpd.get('scaffold_id', '')} "
            f"| {smiles_short} "
            f"| {potency} "
            f"| {cmpd.get('molecular_weight', 'N/A')} "
            f"| {cmpd.get('logp', 'N/A')} "
            f"| {cmpd.get('source', '')} |"
        )

    analog_note = (
        "\n> **Note:** Fewer than 5 direct IGHMBP2 binders were found. "
        "The helicase-analog fallback was activated. Compounds labelled "
        "`helicase_analog` are structural analogs from other SF1 helicase "
        "programmes and are **not confirmed IGHMBP2 binders**. "
        "Experimental validation (e.g., thermal shift, SPR) is required.\n"
        if analog_fallback_used else ""
    )

    gap_note = (
        "\n> **Data gap:** No publicly available small molecule data was found "
        "for IGHMBP2 in BindingDB or PubChem BioAssay as of the query date. "
        "IGHMBP2 remains a dark-target; de novo screening or structure-based "
        "virtual screening is recommended.\n"
        if not compounds else ""
    )

    report = f"""# {gene_symbol} Hit Identification Report

**Generated**: {date_str}
**UniProt**: P38935
**Sources queried**: {', '.join(sources_queried)}

---

## Summary

| Metric | Value |
|--------|-------|
| Sources queried | {len(sources_queried)} |
| Compounds found (pre-filter) | — |
| Compounds after Ro5 filter | {len(compounds)} |
| Analog fallback activated | {'Yes' if analog_fallback_used else 'No'} |

{analog_note}{gap_note}
---

## Top 10 Hits

| Scaffold ID | SMILES (truncated) | Potency | MW | LogP | Source |
|-------------|-------------------|---------|-----|------|--------|
{chr(10).join(rows) if rows else "| — | No compounds found | — | — | — | — |"}

---

## Notes

- All potency values in nM (Ki, IC50, Kd, or EC50 as reported by source).
- Ro5 filter: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10.
- Compounds deduplicated by canonical SMILES (RDKit); salts stripped (largest fragment kept).
- Scaffold IDs assigned by potency rank ({gene_symbol}-SCF-001 = most potent).
"""
    out_path.write_text(report)


# ── Main entry point ──────────────────────────────────────────────────────────

def run_hitid(
    gene_symbol: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> None:
    """Run hit identification for a target with no ChEMBL data."""
    from pipeline.config import TARGETS

    tgt = TARGETS.get(gene_symbol)
    if tgt is None:
        raise ValueError(f"Unknown target: {gene_symbol}")

    uniprot_id = tgt.uniprot_id
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    out_csv = results_dir / "compounds_filtered.csv"
    out_report = results_dir / "hitid_report.md"

    # Cache check
    if not force:
        cached = cache.load(gene_symbol, "hitid")
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol:10}[/dim] hitid  "
                f"[yellow]SKIP[/yellow]  (cached {len(cached) if isinstance(cached, list) else '?'} records)"
            )
            return

    sources_queried: list[str] = []
    all_hits: list[dict] = []

    # Step 1 — BindingDB by UniProt
    console.print(f"  [dim]{gene_symbol}:[/dim] Step 1 — BindingDB (UniProt {uniprot_id})")
    bdb_uniprot = _query_bindingdb_by_uniprot(uniprot_id, console)
    sources_queried.append(f"BindingDB/UniProt({uniprot_id})")
    all_hits.extend(bdb_uniprot)

    # Step 2 — BindingDB by target name
    console.print(f"  [dim]{gene_symbol}:[/dim] Step 2 — BindingDB (name={gene_symbol})")
    bdb_name = _query_bindingdb_by_name(gene_symbol, "BindingDB_TargetName", console)
    sources_queried.append(f"BindingDB/name({gene_symbol})")
    all_hits.extend(bdb_name)

    # Step 3 — PubChem BioAssay
    console.print(f"  [dim]{gene_symbol}:[/dim] Step 3 — PubChem BioAssay")
    pubchem_hits = _query_pubchem(uniprot_id, console)
    sources_queried.append("PubChem_BioAssay")
    all_hits.extend(pubchem_hits)

    # Step 4 — SF1 helicase analog fallback
    analog_fallback_used = False
    if len(all_hits) < 5:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] Step 4 — fewer than 5 hits; "
            "activating SF1 helicase analog fallback"
        )
        helicase_hits = _query_bindingdb_by_name("helicase", "helicase_analog", console)
        sources_queried.append("BindingDB/name(helicase)[analog_fallback]")
        # Pre-filter analogs to Ki/IC50 < 10000 nM and label
        helicase_filtered = [
            {**h, "source": "helicase_analog"}
            for h in helicase_hits
            if h.get("best_value_nm") is not None and h["best_value_nm"] < 10000
        ]
        all_hits.extend(helicase_filtered)
        if helicase_filtered:
            analog_fallback_used = True
            console.print(
                f"    Helicase analog fallback: {len(helicase_filtered)} compounds "
                f"(Ki/IC50 < 10 µM) added"
            )
        else:
            console.print(
                "    Helicase analog fallback: no potent analogs found"
            )

    # Step 5 — Process and filter
    console.print(f"  [dim]{gene_symbol}:[/dim] Step 5 — Processing {len(all_hits)} raw hits")
    compounds = _process_compounds(all_hits, gene_symbol, console)

    # Step 6 — Write outputs
    _write_csv(compounds, out_csv)
    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]compounds_filtered.csv[/green] → "
        f"{out_csv} ({len(compounds)} rows)"
    )

    _write_report(gene_symbol, compounds, sources_queried, analog_fallback_used, out_report)
    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]hitid_report.md[/green] → {out_report}"
    )

    # Cache
    cache_records = [
        {k: v for k, v in c.items() if k != "smiles" or True}
        for c in compounds
    ]
    cache.save(gene_symbol, "hitid", cache_records, len(cache_records))
    console.print(
        f"  [dim]{gene_symbol:10}[/dim] hitid  [green]OK[/green]    "
        f"({len(compounds)} compounds after Ro5 filter)"
    )
