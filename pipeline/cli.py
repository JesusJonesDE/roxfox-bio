from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from pipeline.config import Settings, TARGETS

app = typer.Typer(help="RoxFox Bio drug discovery research pipeline.", add_completion=False)
console = Console()

# ── Shared option defaults ─────────────────────────────────────────────────────

def _resolve_targets(target: Optional[str], all_targets: bool) -> list[str]:
    if all_targets:
        return list(TARGETS.keys())
    if target:
        name = target.upper()
        if name not in TARGETS:
            console.print(f"[red]Unknown target: {target}. Available: {', '.join(TARGETS.keys())}[/red]")
            raise typer.Exit(1)
        return [name]
    console.print("[red]Specify --target <NAME> or --all[/red]")
    raise typer.Exit(1)


def _make_settings(data_dir: Optional[Path], max_age: int) -> Settings:
    s = Settings()
    if data_dir:
        s.data_dir = data_dir
    s.cache_max_age_days = max_age
    return s


# ── fetch ──────────────────────────────────────────────────────────────────────

@app.command()
def fetch(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    force: Optional[str] = typer.Option(None, "--force", help="Force stage: fetch | all"),
    max_age: int = typer.Option(30, "--max-age", help="Cache freshness threshold in days"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Fetch raw data from all external sources for the specified target(s)."""
    from pipeline.cache import CacheManager
    from pipeline.stages.fetch.chembl import fetch_chembl
    from pipeline.stages.fetch.open_targets import fetch_open_targets
    from pipeline.stages.fetch.pdb import fetch_pdb
    from pipeline.stages.fetch.alphafold import fetch_alphafold
    from pipeline.stages.fetch.clinical_trials import fetch_clinical_trials

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)
    force_fetch = force in ("fetch", "all")

    exit_code = 0
    for gene in targets:
        tgt = TARGETS[gene]
        console.rule(f"[bold cyan]{gene}[/bold cyan] — fetch")
        sources = [
            ("chembl", lambda t=tgt: fetch_chembl(t, settings)),
            ("open_targets", lambda t=tgt: fetch_open_targets(t, settings)),
            ("pdb", lambda t=tgt: fetch_pdb(t, settings)),
            ("alphafold", lambda t=tgt: fetch_alphafold(t, settings)),
            ("clinical_trials", lambda t=tgt: fetch_clinical_trials(t, settings)),
        ]
        for source_name, fetcher in sources:
            if not force_fetch and cache.is_source_fresh(gene, source_name):
                meta = cache.get_source_meta(gene, source_name)
                console.print(f"  [dim]{gene:10}[/dim] fetch  {source_name:20} [yellow]SKIP[/yellow]  (cached {meta['record_count']} records)")
                continue
            try:
                import time
                t0 = time.monotonic()
                data, count = fetcher()
                cache.save(gene, source_name, data, count)
                elapsed = time.monotonic() - t0
                console.print(f"  [dim]{gene:10}[/dim] fetch  {source_name:20} [green]OK[/green]    ({count} records, {elapsed:.1f}s)")
            except Exception as exc:
                console.print(f"  [dim]{gene:10}[/dim] fetch  {source_name:20} [red]FAIL[/red]  {exc}")
                exit_code = 1

        if exit_code == 0:
            cache.mark_stage_complete("fetch", gene)

    raise typer.Exit(exit_code)


# ── analyze ────────────────────────────────────────────────────────────────────

@app.command()
def analyze(
    target: Optional[str] = typer.Option(None, "--target", "-t"),
    all_targets: bool = typer.Option(False, "--all"),
    force: Optional[str] = typer.Option(None, "--force"),
    max_age: int = typer.Option(30, "--max-age"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
):
    """Run analysis on cached data: filter compounds, scaffold clustering, structural inventory."""
    from pipeline.cache import CacheManager
    from pipeline.stages.analyze.bioactivity import run_bioactivity_analysis
    from pipeline.stages.analyze.scaffolds import run_scaffold_analysis
    from pipeline.stages.analyze.selectivity import run_selectivity_analysis
    from pipeline.stages.analyze.structures import run_structures_analysis

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)
    force_analyze = force in ("analyze", "all")

    exit_code = 0
    for gene in targets:
        if not cache.is_stage_complete("fetch", gene):
            console.print(f"[red]{gene}: fetch stage not complete. Run `pipeline fetch --target {gene}` first.[/red]")
            exit_code = 1
            continue

        if not force_analyze and cache.is_stage_complete("analyze", gene):
            console.print(f"  [dim]{gene:10}[/dim] analyze [yellow]SKIP[/yellow]  (cached)")
            continue

        console.rule(f"[bold cyan]{gene}[/bold cyan] — analyze")
        try:
            results_dir = settings.results_dir / gene
            results_dir.mkdir(parents=True, exist_ok=True)

            compounds = run_bioactivity_analysis(gene, settings, cache, console)
            compounds = run_scaffold_analysis(gene, compounds, settings, console)
            compounds = run_selectivity_analysis(gene, compounds, settings, cache, console)
            run_structures_analysis(gene, settings, cache, console)
            cache.mark_stage_complete("analyze", gene)
        except Exception as exc:
            console.print(f"  [red]{gene} analyze FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── report ─────────────────────────────────────────────────────────────────────

@app.command()
def report(
    target: Optional[str] = typer.Option(None, "--target", "-t"),
    all_targets: bool = typer.Option(False, "--all"),
    force: Optional[str] = typer.Option(None, "--force"),
    max_age: int = typer.Option(30, "--max-age"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
):
    """Generate markdown research dossier from analysis results."""
    from pipeline.cache import CacheManager
    from pipeline.stages.report.dossier import generate_dossier

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)
    force_report = force in ("report", "all")

    exit_code = 0
    for gene in targets:
        if not cache.is_stage_complete("analyze", gene):
            console.print(f"[red]{gene}: analyze stage not complete. Run `pipeline analyze --target {gene}` first.[/red]")
            exit_code = 1
            continue

        if not force_report and cache.is_stage_complete("report", gene):
            console.print(f"  [dim]{gene:10}[/dim] report  [yellow]SKIP[/yellow]  (cached)")
            continue

        console.rule(f"[bold cyan]{gene}[/bold cyan] — report")
        try:
            out_path = generate_dossier(gene, settings, cache)
            console.print(f"  [dim]{gene:10}[/dim] report  [green]OK[/green]    → {out_path}")
            cache.mark_stage_complete("report", gene)
        except Exception as exc:
            console.print(f"  [red]{gene} report FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── run ────────────────────────────────────────────────────────────────────────

@app.command()
def run(
    target: Optional[str] = typer.Option(None, "--target", "-t"),
    all_targets: bool = typer.Option(False, "--all"),
    force: Optional[str] = typer.Option(None, "--force", help="Force stage: fetch | analyze | report | all"),
    max_age: int = typer.Option(30, "--max-age"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
):
    """Execute all stages in sequence: fetch → analyze → report."""
    from pipeline.cache import CacheManager
    from pipeline.stages.fetch.chembl import fetch_chembl
    from pipeline.stages.fetch.open_targets import fetch_open_targets
    from pipeline.stages.fetch.pdb import fetch_pdb
    from pipeline.stages.fetch.alphafold import fetch_alphafold
    from pipeline.stages.fetch.clinical_trials import fetch_clinical_trials
    from pipeline.stages.analyze.bioactivity import run_bioactivity_analysis
    from pipeline.stages.analyze.scaffolds import run_scaffold_analysis
    from pipeline.stages.analyze.selectivity import run_selectivity_analysis
    from pipeline.stages.analyze.structures import run_structures_analysis
    from pipeline.stages.report.dossier import generate_dossier
    import time

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)
    force_fetch = force in ("fetch", "all")
    force_analyze = force in ("analyze", "all")
    force_report = force in ("report", "all")

    exit_code = 0
    for gene in targets:
        tgt = TARGETS[gene]

        # ── Fetch ──
        console.rule(f"[bold cyan]{gene}[/bold cyan] — fetch")
        sources = [
            ("chembl", lambda t=tgt: fetch_chembl(t, settings)),
            ("open_targets", lambda t=tgt: fetch_open_targets(t, settings)),
            ("pdb", lambda t=tgt: fetch_pdb(t, settings)),
            ("alphafold", lambda t=tgt: fetch_alphafold(t, settings)),
            ("clinical_trials", lambda t=tgt: fetch_clinical_trials(t, settings)),
        ]
        fetch_ok = True
        for source_name, fetcher in sources:
            if not force_fetch and cache.is_source_fresh(gene, source_name):
                meta = cache.get_source_meta(gene, source_name)
                console.print(f"  [dim]{gene:10}[/dim] fetch  {source_name:20} [yellow]SKIP[/yellow]  (cached {meta['record_count']} records)")
                continue
            try:
                t0 = time.monotonic()
                data, count = fetcher()
                cache.save(gene, source_name, data, count)
                elapsed = time.monotonic() - t0
                console.print(f"  [dim]{gene:10}[/dim] fetch  {source_name:20} [green]OK[/green]    ({count} records, {elapsed:.1f}s)")
            except Exception as exc:
                console.print(f"  [dim]{gene:10}[/dim] fetch  {source_name:20} [red]FAIL[/red]  {exc}")
                fetch_ok = False
                exit_code = 1
        if fetch_ok:
            cache.mark_stage_complete("fetch", gene)

        if not cache.is_stage_complete("fetch", gene):
            console.print(f"[red]{gene}: fetch incomplete, skipping analyze + report[/red]")
            continue

        # ── Analyze ──
        console.rule(f"[bold cyan]{gene}[/bold cyan] — analyze")
        if not force_analyze and cache.is_stage_complete("analyze", gene):
            console.print(f"  [dim]{gene:10}[/dim] analyze [yellow]SKIP[/yellow]  (cached)")
        else:
            try:
                results_dir = settings.results_dir / gene
                results_dir.mkdir(parents=True, exist_ok=True)
                compounds = run_bioactivity_analysis(gene, settings, cache, console)
                compounds = run_scaffold_analysis(gene, compounds, settings, console)
                compounds = run_selectivity_analysis(gene, compounds, settings, cache, console)
                run_structures_analysis(gene, settings, cache, console)
                cache.mark_stage_complete("analyze", gene)
            except Exception as exc:
                console.print(f"  [red]{gene} analyze FAIL: {exc}[/red]")
                exit_code = 1
                continue

        # ── Report ──
        console.rule(f"[bold cyan]{gene}[/bold cyan] — report")
        if not force_report and cache.is_stage_complete("report", gene):
            console.print(f"  [dim]{gene:10}[/dim] report  [yellow]SKIP[/yellow]  (cached)")
        else:
            try:
                out_path = generate_dossier(gene, settings, cache)
                console.print(f"  [dim]{gene:10}[/dim] report  [green]OK[/green]    → {out_path}")
                cache.mark_stage_complete("report", gene)
            except Exception as exc:
                console.print(f"  [red]{gene} report FAIL: {exc}[/red]")
                exit_code = 1

    raise typer.Exit(exit_code)


# ── oncology ───────────────────────────────────────────────────────────────────

@app.command()
def oncology(
    target: Optional[str] = typer.Option(None, "--target", "-t"),
    all_targets: bool = typer.Option(False, "--all"),
    series: str = typer.Option("SCF-013,SCF-001", "--series", help="Comma-separated scaffold IDs for SAR + IP sweep"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
    max_age: int = typer.Option(30, "--max-age"),
):
    """Run oncology analysis: OT evidence, kinome selectivity, IP sweep, and SAR."""
    from pipeline.cache import CacheManager
    from pipeline.stages.oncology.ot_oncology import run_ot_oncology
    from pipeline.stages.oncology.selectivity_kinome import run_selectivity_kinome
    from pipeline.stages.oncology.ip_sweep import run_ip_sweep
    from pipeline.stages.oncology.sar_analysis import run_sar_analysis

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)
    series_ids = [s.strip() for s in series.split(",") if s.strip()]

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — oncology")

        steps = [
            ("OT oncology evidence", lambda g=gene: run_ot_oncology(g, settings, cache, console)),
            ("kinome selectivity",   lambda g=gene: run_selectivity_kinome(g, settings, cache, console)),
            ("IP sweep",             lambda g=gene: run_ip_sweep(g, series_ids, settings, cache, console)),
        ]
        for label, fn in steps:
            try:
                fn()
            except Exception as exc:
                console.print(f"  [red]{gene} {label} FAIL: {exc}[/red]")
                exit_code = 1

        for sid in series_ids:
            try:
                run_sar_analysis(gene, sid, settings, console)
            except Exception as exc:
                console.print(f"  [red]{gene} SAR {sid} FAIL: {exc}[/red]")
                exit_code = 1

    raise typer.Exit(exit_code)


# ── pockets ────────────────────────────────────────────────────────────────────

@app.command()
def pockets(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. IGHMBP2)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
    max_age: int = typer.Option(30, "--max-age"),
):
    """Run fpocket druggability analysis on X-ray structures for the specified target(s)."""
    from pipeline.cache import CacheManager
    from pipeline.stages.pockets.pocket_analysis import run_pocket_analysis

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — pockets")
        try:
            run_pocket_analysis(gene, settings, cache, console)
        except Exception as exc:
            console.print(f"  [red]{gene} pockets FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── triage ─────────────────────────────────────────────────────────────────────

@app.command()
def triage(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
    max_age: int = typer.Option(30, "--max-age"),
):
    """Score and rank scaffolds from clean (Ro5-passing, selective) compounds."""
    from pipeline.stages.triage.scaffold_triage import run_scaffold_triage

    settings = _make_settings(data_dir, max_age)
    targets = _resolve_targets(target, all_targets)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — triage")
        try:
            run_scaffold_triage(gene, settings, console)
        except Exception as exc:
            console.print(f"  [red]{gene} triage FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── depmap ─────────────────────────────────────────────────────────────────────

@app.command()
def depmap(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    force: bool = typer.Option(False, "--force", help="Re-fetch even if cache is fresh"),
    max_age: int = typer.Option(30, "--max-age", help="Cache freshness threshold in days"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Fetch DepMap CRISPR gene effect data and produce a ranked cancer lineage report."""
    from pipeline.cache import CacheManager
    from pipeline.stages.depmap.depmap import run_depmap

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — depmap")
        try:
            run_depmap(gene, settings, cache, force, console)
        except Exception as exc:
            console.print(f"  [red]{gene} depmap FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── structalign ────────────────────────────────────────────────────────────────

@app.command()
def structalign(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    include_vrk2: bool = typer.Option(False, "--include-vrk2", help="Include VRK2 in three-way comparison"),
    force: bool = typer.Option(False, "--force", help="Re-run even if output files exist"),
    cutoff: float = typer.Option(6.0, "--cutoff", help="Binding site distance cutoff in Å"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Align VRK1 vs. EGFR ATP binding sites and produce a residue-level selectivity report."""
    from pipeline.cache import CacheManager
    from pipeline.stages.structalign.structalign import run_structalign

    settings = _make_settings(data_dir, 30)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — structalign")
        try:
            run_structalign(gene, settings, cache, force, include_vrk2, cutoff, console)
        except Exception as exc:
            console.print(f"  [red]{gene} structalign FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── biomarker ──────────────────────────────────────────────────────────────────

@app.command()
def biomarker(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    lineage: str = typer.Option(..., "--lineage", help='OncotreeLineage to analyse (e.g. "CNS/Brain")'),
    force: bool = typer.Option(False, "--force", help="Re-run even if cached output exists"),
    min_lines: int = typer.Option(3, "--min-lines", help="Min dependent lines with mutation for inclusion"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Compute co-mutation enrichment biomarker for a target's DepMap lineage."""
    from pipeline.cache import CacheManager
    from pipeline.stages.biomarker.biomarker import run_biomarker

    settings = _make_settings(data_dir, 30)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — biomarker")
        try:
            run_biomarker(gene, lineage, settings, cache, force, min_lines, console)
        except ValueError as exc:
            console.print(f"  [red]{gene} biomarker FAIL: {exc}[/red]")
            exit_code = 1
        except Exception as exc:
            console.print(f"  [red]{gene} biomarker FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── dock ───────────────────────────────────────────────────────────────────────

def _resolve_scaffolds(
    scaffold: Optional[str],
    all_scaffolds: bool,
    gene: str,
    settings: "Settings",
    top_n: Optional[int] = None,
) -> list[str]:
    """Return deduplicated scaffold ID list for a target.

    With --all-scaffolds: deduplicate by scaffold_id, keep only Ro5-passing +
    non-flagged (selectivity_flag==False), sort by best_value_nm descending
    (higher thermal shift first), then apply --top-n if set.
    """
    import pandas as _pd

    if scaffold:
        return [scaffold]

    results_dir = settings.results_dir / gene
    compounds_csv = results_dir / "compounds_filtered.csv"
    if not compounds_csv.exists():
        raise FileNotFoundError(f"compounds_filtered.csv not found for {gene}")

    df = _pd.read_csv(compounds_csv)

    # Deduplicate: one row per scaffold, keep the best (highest best_value_nm)
    id_col = next(
        (c for c in ("scaffold_id", "compound_id", "molecule_chembl_id") if c in df.columns),
        None,
    )
    if id_col is None:
        return []

    df = df.dropna(subset=[id_col])

    # Filter: Ro5-passing and not flagged for selectivity problems
    if "passes_ro5" in df.columns:
        df = df[df["passes_ro5"] == True]  # noqa: E712
    if "selectivity_flag" in df.columns:
        df = df[df["selectivity_flag"] == False]  # noqa: E712

    # One row per scaffold: keep highest thermal shift / best potency value
    if "best_value_nm" in df.columns:
        df = df.sort_values("best_value_nm", ascending=False)
    df = df.drop_duplicates(subset=[id_col])

    scaffolds = df[id_col].tolist()
    if top_n is not None:
        scaffolds = scaffolds[:top_n]
    return scaffolds


@app.command()
def dock(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    scaffold: Optional[str] = typer.Option(None, "--scaffold", help="Scaffold ID from compound library (e.g. SCF-013)"),
    all_scaffolds: bool = typer.Option(False, "--all-scaffolds", help="Dock all Ro5-passing, non-flagged scaffolds"),
    top_n: Optional[int] = typer.Option(None, "--top-n", help="With --all-scaffolds: dock only the top N scaffolds by potency"),
    force: bool = typer.Option(False, "--force", help="Re-run even if cached output exists"),
    exhaustiveness: int = typer.Option(32, "--exhaustiveness", help="Vina exhaustiveness parameter"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Dock a scaffold into the target binding site using AutoDock Vina."""
    from pipeline.cache import CacheManager
    from pipeline.stages.dock.dock import run_dock

    settings = _make_settings(data_dir, 30)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)

    if not scaffold and not all_scaffolds:
        console.print("[red]Specify --scaffold <ID> or --all-scaffolds[/red]")
        raise typer.Exit(1)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — dock")
        try:
            scaffolds_to_run = _resolve_scaffolds(scaffold, all_scaffolds, gene, settings, top_n)
        except FileNotFoundError as exc:
            console.print(f"  [red]{exc}[/red]")
            exit_code = 1
            continue

        if all_scaffolds:
            console.print(f"  [dim]{gene}:[/dim] {len(scaffolds_to_run)} scaffolds queued")

        for scf in scaffolds_to_run:
            try:
                run_dock(gene, scf, settings, cache, force, exhaustiveness, console)
            except RuntimeError as exc:
                console.print(f"  [red]{exc}[/red]")
                exit_code = 1
            except Exception as exc:
                console.print(f"  [red]{gene} dock {scf} FAIL: {exc}[/red]")
                exit_code = 1

    raise typer.Exit(exit_code)


# ── cocrystal ──────────────────────────────────────────────────────────────────

@app.command()
def cocrystal(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Run for all configured targets"),
    scaffold: Optional[str] = typer.Option(None, "--scaffold", help="Scaffold ID (e.g. SCF-013)"),
    all_scaffolds: bool = typer.Option(False, "--all-scaffolds", help="Generate briefs for all Ro5-passing scaffolds"),
    top_n: Optional[int] = typer.Option(None, "--top-n", help="With --all-scaffolds: generate only top N by potency"),
    force: bool = typer.Option(False, "--force", help="Re-generate even if brief already exists"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Generate a co-crystallisation experimental brief for a target/scaffold pair."""
    from pipeline.cache import CacheManager
    from pipeline.stages.cocrystal.cocrystal import run_cocrystal

    settings = _make_settings(data_dir, 30)
    cache = CacheManager(settings)
    targets = _resolve_targets(target, all_targets)

    if not scaffold and not all_scaffolds:
        console.print("[red]Specify --scaffold <ID> or --all-scaffolds[/red]")
        raise typer.Exit(1)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — cocrystal")
        try:
            scaffolds_to_run = _resolve_scaffolds(scaffold, all_scaffolds, gene, settings, top_n)
        except FileNotFoundError as exc:
            console.print(f"  [red]{exc}[/red]")
            exit_code = 1
            continue

        for scf in scaffolds_to_run:
            try:
                run_cocrystal(gene, scf, settings, cache, force, console)
            except Exception as exc:
                console.print(f"  [red]{gene} cocrystal {scf} FAIL: {exc}[/red]")
                exit_code = 1

    raise typer.Exit(exit_code)


# ── rank ───────────────────────────────────────────────────────────────────────

@app.command()
def rank(
    target: Optional[str] = typer.Option(None, "--target", "-t", help="Target gene name (e.g. VRK1)"),
    all_targets: bool = typer.Option(False, "--all", help="Rank scaffolds for all configured targets"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir", help="Override data directory"),
):
    """Rank all docked scaffolds by affinity and produce a cross-scaffold summary table."""
    from pipeline.stages.rank.rank import run_rank

    settings = _make_settings(data_dir, 30)
    targets = _resolve_targets(target, all_targets)

    exit_code = 0
    for gene in targets:
        console.rule(f"[bold cyan]{gene}[/bold cyan] — rank")
        try:
            run_rank(gene, settings, console)
        except Exception as exc:
            console.print(f"  [red]{gene} rank FAIL: {exc}[/red]")
            exit_code = 1

    raise typer.Exit(exit_code)


# ── status ─────────────────────────────────────────────────────────────────────

@app.command()
def status(
    max_age: int = typer.Option(30, "--max-age"),
    data_dir: Optional[Path] = typer.Option(None, "--data-dir"),
):
    """Show cache and manifest state for all configured targets."""
    from pipeline.cache import CacheManager

    settings = _make_settings(data_dir, max_age)
    cache = CacheManager(settings)

    table = Table(title="Pipeline Status", show_header=True)
    table.add_column("Target", style="bold")
    table.add_column("Fetch")
    table.add_column("Analyze")
    table.add_column("Report")
    table.add_column("Status")

    for gene in TARGETS:
        def fmt(stage: str) -> str:
            date = cache.get_stage_date(stage, gene)
            if date is None:
                return "[dim]—[/dim]"
            complete = cache.is_stage_complete(stage, gene)
            label = date[:10]
            return f"[green]{label}[/green]" if complete else f"[yellow]{label} STALE[/yellow]"

        fetch_ok = cache.is_stage_complete("fetch", gene)
        analyze_ok = cache.is_stage_complete("analyze", gene)
        report_ok = cache.is_stage_complete("report", gene)

        if report_ok:
            overall = "[green]✓ current[/green]"
        elif fetch_ok or analyze_ok:
            overall = "[yellow]⚠ partial[/yellow]"
        else:
            overall = "[red]✗ not run[/red]"

        table.add_row(gene, fmt("fetch"), fmt("analyze"), fmt("report"), overall)

    console.print(table)
