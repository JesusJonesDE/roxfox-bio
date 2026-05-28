from __future__ import annotations

from pipeline.cache import CacheManager
from pipeline.config import Settings
from pipeline.models import GateResult, GateStatus
from pipeline.stages.validate.validate import (
    _cache_gate_result,
    _load_cached_gate_result,
    _load_smiles,
    _write_gate_report,
)

# ── Thresholds ─────────────────────────────────────────────────────────────────

# Field names as returned by admet-ai >= 1.0 (TDC benchmark naming)
_FIELD_MAP = {
    "BBB": "BBB_Martins",
    "CYP1A2": "CYP1A2_Veith",
    "CYP2D6": "CYP2D6_Veith",
    "CYP3A4": "CYP3A4_Veith",
    "Solubility": "Solubility_AqSolDB",
    "HIA": "HIA_Hou",
    "Bioavailability": "Bioavailability_Ma",
}

_THRESHOLDS = {
    "BBB_Martins":         (">",  0.5),   # BBB penetration probability
    "CYP1A2_Veith":        ("<",  0.3),   # CYP1A2 inhibition probability
    "CYP2D6_Veith":        ("<",  0.3),   # CYP2D6 inhibition probability
    "CYP3A4_Veith":        ("<",  0.3),   # CYP3A4 inhibition probability
    "Solubility_AqSolDB":  (">", -4.0),   # logS (> -4 ≈ > 10 µg/mL at MW ~400)
    "HIA_Hou":             (">",  0.3),   # Human intestinal absorption
}


def _check_threshold(value: float, operator: str, threshold: float) -> bool:
    if operator == ">":
        return value > threshold
    if operator == "<":
        return value < threshold
    raise ValueError(f"Unknown operator: {operator!r}")


def run_admet_gate(
    gene_symbol: str,
    scaffold_id: str,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console,
) -> GateResult:
    """Run the ADMET gate for a given scaffold.

    Checks six ADMET properties from ADMET-AI and returns a GateResult with
    PASS if all thresholds are met, FAIL otherwise.
    """
    # 1. Cache check
    if not force:
        cached = _load_cached_gate_result(gene_symbol, scaffold_id, "admet", cache)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] ADMET gate — [dim]SKIP[/dim] (cached result for {scaffold_id})"
            )
            return cached

    # 2. Load SMILES (salt-stripped)
    smiles = _load_smiles(gene_symbol, scaffold_id, settings)

    # 3. Run ADMET-AI
    try:
        from admet_ai import ADMETModel
    except ImportError as exc:
        raise RuntimeError(
            "admet-ai not installed. Run: pip install admet-ai"
        ) from exc

    model = ADMETModel()
    raw_preds = model.predict(smiles=smiles)

    # raw_preds may be a dict keyed by SMILES or a flat dict depending on version
    if isinstance(raw_preds, dict) and smiles in raw_preds:
        preds: dict = raw_preds[smiles]
    elif isinstance(raw_preds, dict):
        # flat dict or single-entry dict — use as-is or unwrap first value
        first_val = next(iter(raw_preds.values()))
        if isinstance(first_val, dict):
            preds = first_val
        else:
            preds = raw_preds
    else:
        preds = raw_preds

    # 4. Extract scores — check all threshold fields are present
    missing = [f for f in _THRESHOLDS if f not in preds]
    if missing:
        raise RuntimeError(
            f"ADMET-AI returned unexpected fields. Missing: {missing}. "
            f"Available: {sorted(preds.keys())[:15]}..."
        )

    scores: dict[str, float] = {field: float(preds[field]) for field in _THRESHOLDS}

    # 5. Apply thresholds
    failures: list[str] = []
    for prop, (op, threshold) in _THRESHOLDS.items():
        value = scores[prop]
        if not _check_threshold(value, op, threshold):
            failures.append(f"{prop}={value:.3g} (threshold {op}{threshold})")

    # 6. Build GateResult
    if failures:
        status = GateStatus.FAIL
        reason = "Failed: " + ", ".join(failures)
    else:
        status = GateStatus.PASS
        reason = "All 6 ADMET properties passed"

    result = GateResult(
        gate_name="admet",
        status=status,
        score=scores["BBB_Martins"],
        reason=reason,
        details=scores,
    )

    # 7. Write report
    result.report_path = _write_gate_report(
        gene_symbol, scaffold_id, "admet", result, settings
    )

    # 8. Cache
    _cache_gate_result(gene_symbol, scaffold_id, "admet", result, cache)

    return result
