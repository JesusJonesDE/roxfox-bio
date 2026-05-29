"""
Fragment pipeline — ADMET screening and final output step.

Evaluates grown candidates with relaxed ADMET thresholds (BBB 0.3 instead
of 0.5) and writes compounds_filtered.csv in the VRK1 schema.
"""
from __future__ import annotations

import pandas as pd
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings

# Column order matching the downstream schema used by pipeline dock / rank
_OUTPUT_COLUMNS = [
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

# Relaxed thresholds for fragment-origin CNS candidates
_RELAXED_THRESHOLDS = {
    "BBB_Martins":        (">",  0.3),   # relaxed from 0.5 — SMARD1 is neurological
    "CYP1A2_Veith":       ("<",  0.3),
    "CYP2D6_Veith":       ("<",  0.3),
    "CYP3A4_Veith":       ("<",  0.3),
    "Solubility_AqSolDB": (">", -4.0),
    "HIA_Hou":            (">",  0.3),
}


def _check_threshold(value: float, operator: str, threshold: float) -> bool:
    if operator == ">":
        return value > threshold
    if operator == "<":
        return value < threshold
    raise ValueError(f"Unknown operator: {operator!r}")


def _run_admet_relaxed(smiles: str) -> tuple[bool, float, dict]:
    """Run ADMET-AI with relaxed BBB threshold (0.3).

    Returns (admet_pass, bbb_score, all_scores).
    """
    try:
        from admet_ai import ADMETModel
    except ImportError as exc:
        raise RuntimeError(
            "admet-ai not installed. Run: pip install admet-ai"
        ) from exc

    model = ADMETModel()
    raw_preds = model.predict(smiles=smiles)

    # Unwrap SMILES-keyed or flat dict
    if isinstance(raw_preds, dict) and smiles in raw_preds:
        preds: dict = raw_preds[smiles]
    elif isinstance(raw_preds, dict):
        first_val = next(iter(raw_preds.values()))
        if isinstance(first_val, dict):
            preds = first_val
        else:
            preds = raw_preds
    else:
        preds = raw_preds

    bbb_score = float(preds.get("BBB_Martins", 0.0))

    all_scores: dict[str, float] = {}
    admet_pass = True
    for field, (op, threshold) in _RELAXED_THRESHOLDS.items():
        value = float(preds.get(field, 0.0))
        all_scores[field] = value
        if not _check_threshold(value, op, threshold):
            admet_pass = False

    return admet_pass, bbb_score, all_scores


def run_output(
    gene_symbol: str,
    candidates_df: pd.DataFrame,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> int:
    """Screen grown candidates with relaxed ADMET and write compounds_filtered.csv.

    Returns the number of rows written.
    """
    # 1. Cache check
    cache_key = "fragment_output"
    if not force:
        cached = cache.load(gene_symbol, cache_key)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] output [yellow]SKIP[/yellow] (cached)"
            )
            return int(cached.get("n_rows", 0)) if isinstance(cached, dict) else 0

    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / "compounds_filtered.csv"

    # 2. Handle empty candidates_df
    if candidates_df is None or len(candidates_df) == 0:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [yellow]no candidates — writing empty "
            f"compounds_filtered.csv[/yellow]"
        )
        empty_df = pd.DataFrame(columns=_OUTPUT_COLUMNS)
        empty_df.to_csv(out_path, index=False)
        cache.save(gene_symbol, cache_key, {"n_rows": 0}, 0)
        return 0

    console.print(
        f"  [dim]{gene_symbol}:[/dim] running relaxed ADMET on "
        f"{len(candidates_df)} candidates (BBB threshold 0.3)..."
    )

    # 3. Evaluate each candidate
    output_rows: list[dict] = []

    for _, cand_row in candidates_df.iterrows():
        smiles = str(cand_row.get("smiles", ""))
        candidate_id = str(cand_row.get("candidate_id", ""))

        if not smiles or not candidate_id:
            continue

        bbb_score = 0.0
        admet_pass = False
        try:
            admet_pass, bbb_score, _ = _run_admet_relaxed(smiles)
        except Exception as exc:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] [yellow]ADMET failed for "
                f"{candidate_id}: {exc}[/yellow]"
            )
            # Include candidate anyway — downstream tools decide

        output_rows.append(
            {
                "molecule_id": candidate_id,
                "smiles": smiles,
                "best_value_nm": None,
                "best_assay_type": "fragment_screen_predicted",
                "molecular_weight": cand_row.get("molecular_weight"),
                "logp": cand_row.get("logp"),
                "hbd": cand_row.get("hbd"),
                "hba": cand_row.get("hba"),
                "rotatable_bonds": cand_row.get("rotatable_bonds"),
                "ro5_violations": cand_row.get("ro5_violations"),
                "passes_ro5": cand_row.get("passes_ro5"),
                "scaffold_id": candidate_id,
                "source": "fragment_virtual_screen",
                "off_target_flags": 0,
                "selectivity_flag": False,
                # Internal sort key — removed before writing
                "_bbb_score": bbb_score,
            }
        )

    # 4. Sort by BBB score descending
    output_rows.sort(key=lambda r: r["_bbb_score"], reverse=True)

    # 5. Build final DataFrame with exact column order
    for row in output_rows:
        row.pop("_bbb_score", None)

    out_df = pd.DataFrame(output_rows, columns=_OUTPUT_COLUMNS)

    # 6. Write compounds_filtered.csv
    out_df.to_csv(out_path, index=False)
    n_rows = len(out_df)

    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]output complete — "
        f"{n_rows} candidates written to compounds_filtered.csv[/green]"
    )

    # 7. Cache and return
    cache.save(gene_symbol, cache_key, {"n_rows": n_rows}, n_rows)
    return n_rows
