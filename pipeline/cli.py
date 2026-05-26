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
