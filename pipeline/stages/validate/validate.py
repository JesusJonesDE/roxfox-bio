from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console
from rich.table import Table

from pipeline.cache import CacheManager
from pipeline.config import Settings
from pipeline.models import GateResult, GateStatus, ValidationResult

_GATE_NAMES = ("admet", "mmgbsa", "selectivity", "md")
_STATUS_COLOURS = {
    GateStatus.PASS: "green",
    GateStatus.FAIL: "red",
    GateStatus.ERROR: "yellow",
    GateStatus.PENDING: "cyan",
    GateStatus.NOT_RUN: "dim",
}


# ── SMILES helpers ─────────────────────────────────────────────────────────────

def _strip_salts(smiles: str) -> str:
    """Return the largest fragment from a multi-fragment SMILES (salt stripping)."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def _load_smiles(gene_symbol: str, scaffold_id: str, settings: Settings) -> str:
    """Load SMILES for a scaffold from compounds_filtered.csv, stripping salts."""
    csv = settings.results_dir / gene_symbol / "compounds_filtered.csv"
    if not csv.exists():
        raise FileNotFoundError(
            f"compounds_filtered.csv not found for {gene_symbol}. "
            f"Run `pipeline fetch --target {gene_symbol}` first."
        )
    df = pd.read_csv(csv)
    row = df[df.get("scaffold_id", pd.Series(dtype=str)) == scaffold_id]
    if row.empty:
        for col in ("compound_id", "molecule_chembl_id"):
            if col in df.columns:
                row = df[df[col] == scaffold_id]
                if not row.empty:
                    break
    if row.empty:
        raise ValueError(f"Scaffold '{scaffold_id}' not found in compounds_filtered.csv")
    smiles_col = next(
        (c for c in ("smiles", "canonical_smiles", "SMILES") if c in row.columns), None
    )
    raw = str(row.iloc[0][smiles_col]) if smiles_col else ""
    return _strip_salts(raw)


# ── Docking PDBQT helper ───────────────────────────────────────────────────────

def _load_docking_pdbqt(gene_symbol: str, scaffold_id: str, settings: Settings) -> Path:
    """Locate top-pose docking PDBQT. Raises FileNotFoundError with helpful message."""
    pdbqt = settings.results_dir / gene_symbol / f"docking_poses_{scaffold_id}.pdbqt"
    if not pdbqt.exists():
        raise FileNotFoundError(
            f"Docking result not found for {scaffold_id}: {pdbqt}\n"
            f"Run `pipeline dock --target {gene_symbol} --scaffold {scaffold_id}` first."
        )
    return pdbqt


# ── Cache helpers ──────────────────────────────────────────────────────────────

def _cache_key(gate_name: str, scaffold_id: str) -> str:
    return f"validate_{gate_name}_{scaffold_id}"


def _cache_gate_result(
    gene_symbol: str,
    scaffold_id: str,
    gate_name: str,
    result: GateResult,
    cache: CacheManager,
) -> None:
    payload = {
        "gate_name": result.gate_name,
        "status": result.status.value,
        "score": result.score,
        "reason": result.reason,
        "details": result.details,
        "report_path": str(result.report_path) if result.report_path else None,
        "duration_s": result.duration_s,
        "timestamp": result.timestamp,
    }
    cache.save(gene_symbol, _cache_key(gate_name, scaffold_id), payload, 1)


def _load_cached_gate_result(
    gene_symbol: str,
    scaffold_id: str,
    gate_name: str,
    cache: CacheManager,
) -> Optional[GateResult]:
    data = cache.load(gene_symbol, _cache_key(gate_name, scaffold_id))
    if data is None:
        return None
    return GateResult(
        gate_name=data["gate_name"],
        status=GateStatus(data["status"]),
        score=data["score"],
        reason=data["reason"],
        details=data.get("details", {}),
        report_path=Path(data["report_path"]) if data.get("report_path") else None,
        duration_s=data.get("duration_s", 0.0),
        timestamp=data.get("timestamp", ""),
    )


# ── Report writer ──────────────────────────────────────────────────────────────

def _write_gate_report(
    gene_symbol: str,
    scaffold_id: str,
    gate_name: str,
    result: GateResult,
    settings: Settings,
    extra_sections: str = "",
) -> Path:
    """Write a standard gate markdown report. Returns path written."""
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / f"validate_{gate_name}_{scaffold_id}.md"

    status_icon = "✅ PASS" if result.status == GateStatus.PASS else (
        "❌ FAIL" if result.status == GateStatus.FAIL else
        "⚠️ ERROR" if result.status == GateStatus.ERROR else str(result.status.value)
    )

    detail_rows = "\n".join(
        f"| {k} | {v} |" for k, v in result.details.items()
    )

    report = f"""# {gene_symbol} — {gate_name.upper()} Gate — {scaffold_id}

**Status**: {status_icon}
**Score**: {result.score:.4g}
**Reason**: {result.reason}
**Duration**: {result.duration_s:.1f}s
**Generated**: {result.timestamp}

---

## Details

| Property | Value |
|----------|-------|
{detail_rows if detail_rows else "| — | — |"}

{extra_sections}
"""
    path.write_text(report)
    return path


# ── ValidationResult persistence ──────────────────────────────────────────────

def _save_validation_result(result: ValidationResult, settings: Settings) -> None:
    path = (
        settings.results_dir / result.gene_symbol
        / f"validation_result_{result.scaffold_id}.json"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "gene_symbol": result.gene_symbol,
        "scaffold_id": result.scaffold_id,
        "smiles": result.smiles,
        "overall_pass": result.overall_pass,
        "handoff_ready": result.handoff_ready,
        "created_at": result.created_at,
        "gates": {
            name: {
                "status": gr.status.value,
                "score": gr.score,
                "reason": gr.reason,
                "details": gr.details,
                "report_path": str(gr.report_path) if gr.report_path else None,
            }
            for name, gr in result.gates.items()
        },
    }
    path.write_text(json.dumps(data, indent=2))


def _load_validation_result(
    gene_symbol: str, scaffold_id: str, settings: Settings
) -> Optional[ValidationResult]:
    path = (
        settings.results_dir / gene_symbol
        / f"validation_result_{scaffold_id}.json"
    )
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    gates = {
        name: GateResult(
            gate_name=name,
            status=GateStatus(g["status"]),
            score=g["score"],
            reason=g["reason"],
            details=g.get("details", {}),
            report_path=Path(g["report_path"]) if g.get("report_path") else None,
        )
        for name, g in data.get("gates", {}).items()
    }
    return ValidationResult(
        gene_symbol=data["gene_symbol"],
        scaffold_id=data["scaffold_id"],
        smiles=data["smiles"],
        gates=gates,
        overall_pass=data["overall_pass"],
        handoff_ready=data["handoff_ready"],
        created_at=data["created_at"],
    )


def _load_all_validation_results(
    gene_symbol: str, settings: Settings
) -> dict[str, ValidationResult]:
    results_dir = settings.results_dir / gene_symbol
    out: dict[str, ValidationResult] = {}
    for path in sorted(results_dir.glob("validation_result_*.json")):
        scaffold_id = path.stem.replace("validation_result_", "")
        r = _load_validation_result(gene_symbol, scaffold_id, settings)
        if r:
            out[scaffold_id] = r
    return out


# ── Dashboard ──────────────────────────────────────────────────────────────────

def _render_dashboard(
    gene_symbol: str,
    results: dict[str, ValidationResult],
    console: Console,
) -> None:
    table = Table(title=f"{gene_symbol} — Validation Dashboard", show_header=True)
    table.add_column("Scaffold", style="bold")
    for gate in _GATE_NAMES:
        table.add_column(gate.upper())
    table.add_column("Handoff?")

    for scaffold_id, r in sorted(results.items()):
        cells = []
        for gate in _GATE_NAMES:
            gr = r.gates.get(gate)
            if gr is None:
                cells.append("[dim]—[/dim]")
            else:
                colour = _STATUS_COLOURS.get(gr.status, "white")
                cells.append(f"[{colour}]{gr.status.value}[/{colour}]")
        handoff = "[green]✓[/green]" if r.handoff_ready else "[dim]✗[/dim]"
        table.add_row(scaffold_id, *cells, handoff)

    console.print(table)


def _write_dashboard_md(
    gene_symbol: str,
    results: dict[str, ValidationResult],
    settings: Settings,
) -> Path:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    header = f"| Scaffold | ADMET | MM-GBSA | Selectivity | MD | Handoff? |\n|---|---|---|---|---|---|\n"
    rows = []
    for scaffold_id, r in sorted(results.items()):
        cells = []
        for gate in _GATE_NAMES:
            gr = r.gates.get(gate)
            cells.append(gr.status.value if gr else "—")
        handoff = "✓" if r.handoff_ready else "✗"
        rows.append(f"| {scaffold_id} | " + " | ".join(cells) + f" | {handoff} |")

    md_path = results_dir / "validation_dashboard.md"
    md_path.write_text(
        f"# {gene_symbol} — Validation Dashboard\n\n"
        f"**Generated**: {datetime.now().strftime('%Y-%m-%d')}\n\n"
        + header + "\n".join(rows) + "\n"
    )

    json_path = results_dir / "validation_dashboard.json"
    json_path.write_text(json.dumps(
        {scaffold_id: {
            "overall_pass": r.overall_pass,
            "handoff_ready": r.handoff_ready,
            "gates": {g: r.gates[g].status.value for g in r.gates},
        } for scaffold_id, r in results.items()},
        indent=2,
    ))

    return md_path


def _write_wetlab_handoff(
    gene_symbol: str,
    scaffold_id: str,
    result: ValidationResult,
    settings: Settings,
) -> Path:
    results_dir = settings.results_dir / gene_symbol
    admet_gr = result.gates.get("admet")
    mmgbsa_gr = result.gates.get("mmgbsa")
    sel_gr = result.gates.get("selectivity")
    md_gr = result.gates.get("md")

    admet_rows = "\n".join(
        f"| {k} | {v} |" for k, v in (admet_gr.details if admet_gr else {}).items()
    )

    report = f"""# Wet-Lab Handoff Report — {gene_symbol} / {scaffold_id}

**Generated**: {datetime.now().strftime('%Y-%m-%d')}
**Target**: {gene_symbol}
**Scaffold**: {scaffold_id}
**SMILES**: `{result.smiles}`
**Status**: ✅ ALL GATES PASSED — ready for CRO submission

---

## Gate Summary

| Gate | Status | Score | Key metric |
|---|---|---|---|
| ADMET | {admet_gr.status.value if admet_gr else '—'} | {admet_gr.score:.3f if admet_gr else '—'} | BBB penetration |
| MM-GBSA | {mmgbsa_gr.status.value if mmgbsa_gr else '—'} | {mmgbsa_gr.score:.2f if mmgbsa_gr else '—'} kcal/mol | ΔG binding |
| Selectivity | {sel_gr.status.value if sel_gr else '—'} | {sel_gr.score:.1f if sel_gr else '—'}× | Selectivity index |
| MD | {md_gr.status.value if md_gr else '—'} | {md_gr.score:.2f if md_gr else '—'} Å | Mean RMSD |

---

## ADMET Properties

| Property | Score |
|---|---|
{admet_rows if admet_rows else "| — | — |"}

---

## Binding Affinity (MM-GBSA)

ΔG = **{mmgbsa_gr.score:.2f} kcal/mol** (pass threshold: ≤ −7.0 kcal/mol)

---

## Selectivity

Selectivity index = **{sel_gr.score:.1f}×** vs. best off-target (pass threshold: ≥ 10×)

Details: {sel_gr.reason if sel_gr else '—'}

---

## MD Pose Stability

Mean ligand RMSD over final 10 ns = **{md_gr.score:.2f} Å** (pass threshold: ≤ 3.0 Å)

---

## Gate Reports

{"- [ADMET report](validate_admet_" + scaffold_id + ".md)" if admet_gr and admet_gr.report_path else ""}
{"- [MM-GBSA report](validate_mmgbsa_" + scaffold_id + ".md)" if mmgbsa_gr and mmgbsa_gr.report_path else ""}
{"- [Selectivity report](validate_selectivity_" + scaffold_id + ".md)" if sel_gr and sel_gr.report_path else ""}
{"- [MD report](validate_md_" + scaffold_id + ".md)" if md_gr and md_gr.report_path else ""}
- [Docking report](docking_report_{scaffold_id}.md)

---

## Recommended next steps

1. Submit to CRO for biochemical IC50 assay against {gene_symbol} kinase
2. Run kinase selectivity panel (30–50 kinases) at 1 µM
3. Confirm BBB penetration in PAMPA-BBB assay
4. Request cellular target engagement (VRK1 phospho-BAF Western blot in GBM line)
"""
    path = results_dir / f"wetlab_handoff_{scaffold_id}.md"
    path.write_text(report)
    return path


def run_dashboard(
    gene_symbol: str,
    settings: Settings,
    console: Console,
) -> None:
    results = _load_all_validation_results(gene_symbol, settings)
    if not results:
        console.print(
            f"  [yellow]{gene_symbol}:[/yellow] no validation results found. "
            f"Run `pipeline validate --target {gene_symbol} --scaffold <ID>` first."
        )
        return
    _render_dashboard(gene_symbol, results, console)
    md_path = _write_dashboard_md(gene_symbol, results, settings)
    console.print(f"  [dim]{gene_symbol}:[/dim] [green]validation_dashboard.md written[/green]")

    for scaffold_id, r in results.items():
        if r.handoff_ready:
            _write_wetlab_handoff(gene_symbol, scaffold_id, r, settings)
            console.print(
                f"  [dim]{gene_symbol}:[/dim] [green]wetlab_handoff_{scaffold_id}.md written[/green]"
            )


# ── Orchestrator ───────────────────────────────────────────────────────────────

def run_validate(
    gene_symbol: str,
    scaffold_id: str,
    gate: Optional[str],
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
    md_max_cost: float = 5.0,
) -> int:
    """Run validation gates. Returns exit code (0=ok, 1=error)."""
    from pipeline.stages.validate.gates.admet import run_admet_gate
    from pipeline.stages.validate.gates.mmgbsa import run_mmgbsa_gate
    from pipeline.stages.validate.gates.selectivity import run_selectivity_gate
    from pipeline.stages.validate.gates.md import run_md_gate

    gate_runners = {
        "admet": lambda: run_admet_gate(gene_symbol, scaffold_id, settings, cache, force, console),
        "mmgbsa": lambda: run_mmgbsa_gate(gene_symbol, scaffold_id, settings, cache, force, console),
        "selectivity": lambda: run_selectivity_gate(gene_symbol, scaffold_id, settings, cache, force, console),
        "md": lambda: run_md_gate(gene_symbol, scaffold_id, settings, cache, force, console, md_max_cost),
    }

    smiles = _load_smiles(gene_symbol, scaffold_id, settings)
    result = _load_validation_result(gene_symbol, scaffold_id, settings) or ValidationResult(
        gene_symbol=gene_symbol,
        scaffold_id=scaffold_id,
        smiles=smiles,
        created_at=datetime.now(timezone.utc).isoformat(),
    )
    result.smiles = smiles

    exit_code = 0
    gates_to_run = [gate] if gate else list(_GATE_NAMES)

    for gate_name in gates_to_run:
        runner = gate_runners.get(gate_name)
        if runner is None:
            console.print(f"  [red]Unknown gate: {gate_name}. Choose from: {', '.join(_GATE_NAMES)}[/red]")
            return 1

        # In sequential mode, skip MD if any prior gate failed
        if gate is None and gate_name == "md":
            prior = [result.gates.get(g) for g in ("admet", "mmgbsa", "selectivity")]
            if any(gr and gr.status in (GateStatus.FAIL, GateStatus.ERROR) for gr in prior):
                console.print(
                    f"  [dim]{gene_symbol}:[/dim] [dim]MD gate skipped — prior gate did not pass[/dim]"
                )
                result.gates["md"] = GateResult(
                    gate_name="md", status=GateStatus.NOT_RUN,
                    score=float("nan"), reason="Skipped — prior gate did not pass",
                )
                continue

        console.print(f"  [dim]{gene_symbol}:[/dim] running {gate_name.upper()} gate for {scaffold_id}...")
        t0 = time.monotonic()
        try:
            gr = runner()
            gr.duration_s = time.monotonic() - t0
            gr.timestamp = datetime.now(timezone.utc).isoformat()
        except Exception as exc:
            gr = GateResult(
                gate_name=gate_name,
                status=GateStatus.ERROR,
                score=float("nan"),
                reason=str(exc),
                duration_s=time.monotonic() - t0,
                timestamp=datetime.now(timezone.utc).isoformat(),
            )
            exit_code = 1

        result.gates[gate_name] = gr
        _cache_gate_result(gene_symbol, scaffold_id, gate_name, gr, cache)

        icon = {"PASS": "✅", "FAIL": "❌", "ERROR": "⚠️"}.get(gr.status.value, "")
        console.print(
            f"  [dim]{gene_symbol}:[/dim] {gate_name.upper()} {icon} {gr.status.value} "
            f"— {gr.reason} ({gr.duration_s:.1f}s)"
        )

        if gate is None and gr.status == GateStatus.ERROR:
            console.print(f"  [red]{gene_symbol}: ERROR in {gate_name} gate — stopping sequence[/red]")
            break

    # Update overall status
    run_gates = [gr for gr in result.gates.values() if gr.status != GateStatus.NOT_RUN]
    result.overall_pass = bool(run_gates) and all(gr.status == GateStatus.PASS for gr in run_gates)
    result.handoff_ready = all(
        result.gates.get(g) and result.gates[g].status == GateStatus.PASS
        for g in _GATE_NAMES
    )
    _save_validation_result(result, settings)

    if result.handoff_ready:
        _write_wetlab_handoff(gene_symbol, scaffold_id, result, settings)
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [green bold]ALL GATES PASSED — wetlab_handoff_{scaffold_id}.md written[/green bold]"
        )

    return exit_code
