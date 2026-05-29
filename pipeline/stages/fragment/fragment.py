from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings

_SASCORER_PATH = (
    Path(sys.executable).parent.parent
    / "lib"
    / f"python{sys.version_info.major}.{sys.version_info.minor}"
    / "site-packages/rdkit/Contrib/SA_Score"
)

_STEPS = ("pocket", "library", "dock", "cluster", "grow", "admet")


@dataclass
class FragmentState:
    gene_symbol: str
    step_pocket: bool = False
    step_library: bool = False
    step_dock: bool = False
    step_cluster: bool = False
    step_grow: bool = False
    step_admet: bool = False
    n_fragments_docked: int = 0
    n_hits: int = 0
    n_clusters: int = 0
    n_grown: int = 0
    n_candidates_final: int = 0
    library_fallback_used: bool = False
    completed_at: str = ""


def _state_path(gene_symbol: str, settings: Settings) -> Path:
    return settings.cache_dir / gene_symbol / "fragment_state.json"


def _load_state(gene_symbol: str, settings: Settings) -> FragmentState:
    path = _state_path(gene_symbol, settings)
    if path.exists():
        try:
            data = json.loads(path.read_text())
            return FragmentState(**{k: v for k, v in data.items() if k in FragmentState.__dataclass_fields__})
        except Exception:
            pass
    return FragmentState(gene_symbol=gene_symbol)


def _save_state(state: FragmentState, settings: Settings) -> None:
    path = _state_path(state.gene_symbol, settings)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(state), indent=2))


def _write_fragment_report(
    gene_symbol: str,
    state: FragmentState,
    pocket: Optional[dict],
    hits_df: Optional[pd.DataFrame],
    clusters_df: Optional[pd.DataFrame],
    candidates_df: Optional[pd.DataFrame],
    settings: Settings,
) -> Path:
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    pocket_section = "No pocket identified yet." if not pocket else (
        f"| Volume | {pocket.get('volume_A3', '?'):.1f} Å³ |\n"
        f"| Druggability score | {pocket.get('druggability_score', '?'):.3f} |\n"
        f"| Centroid | ({pocket.get('centroid_x', 0):.1f}, {pocket.get('centroid_y', 0):.1f}, {pocket.get('centroid_z', 0):.1f}) |\n"
        f"| Box size | {pocket.get('box_size_A', 20):.0f} Å |\n"
        f"| Mean pLDDT | {pocket.get('plddt_mean', 'N/A')} |"
    )

    hits_section = "Fragment screen not yet complete."
    if hits_df is not None and len(hits_df) > 0:
        top10 = hits_df.head(10)
        rows = [
            f"| {r.get('fragment_id','?')} | {str(r.get('smiles',''))[:40]} | {r.get('affinity_kcal_mol', 0):.2f} | {r.get('cluster_id', '—')} |"
            for _, r in top10.iterrows()
        ]
        hits_section = "| Fragment ID | SMILES | Affinity (kcal/mol) | Cluster |\n|---|---|---|---|\n" + "\n".join(rows)

    clusters_section = "Clustering not yet complete."
    if clusters_df is not None and len(clusters_df) > 0:
        n_clusters = clusters_df["cluster_id"].nunique() if "cluster_id" in clusters_df.columns else 0
        clusters_section = f"{n_clusters} distinct chemotype clusters identified from {len(clusters_df)} top hits."

    candidates_section = "Fragment growing not yet complete."
    if candidates_df is not None and len(candidates_df) > 0:
        top10c = candidates_df.head(10)
        rows = [
            f"| {r.get('candidate_id','?')} | {str(r.get('smiles',''))[:40]} | {r.get('molecular_weight', 0):.0f} | {r.get('sa_score', 0):.1f} | {r.get('grow_method','?')} |"
            for _, r in top10c.iterrows()
        ]
        candidates_section = "| Candidate | SMILES | MW | SA score | Method |\n|---|---|---|---|---|\n" + "\n".join(rows)

    report = f"""# IGHMBP2 Fragment Virtual Screening Report

**Target**: {gene_symbol} (UniProt P38935, SMARD1)
**Generated**: {datetime.now().strftime('%Y-%m-%d')}
**Pipeline**: Fragment-based virtual screen against AlphaFold2 structure

---

## Pipeline Status

| Step | Status |
|---|---|
| Pocket identification | {'✅ Complete' if state.step_pocket else '⏳ Pending'} |
| Fragment library | {'✅ Complete' if state.step_library else '⏳ Pending'} |
| Fragment docking | {'✅ ' + str(state.n_fragments_docked) + ' fragments docked' if state.step_dock else '⏳ Pending'} |
| Clustering | {'✅ ' + str(state.n_clusters) + ' clusters' if state.step_cluster else '⏳ Pending'} |
| Fragment growing | {'✅ ' + str(state.n_grown) + ' candidates grown' if state.step_grow else '⏳ Pending'} |
| ADMET screening | {'✅ ' + str(state.n_candidates_final) + ' candidates' if state.step_admet else '⏳ Pending'} |

{'> ⚠️ **Note**: Fragment library fallback was used (ZINC download failed). Results based on 500-compound bundled library.' if state.library_fallback_used else ''}

---

## Binding Pocket

| Property | Value |
|---|---|
{pocket_section}

---

## Top Fragment Hits

{hits_section}

---

## Chemical Diversity (Clusters)

{clusters_section}

---

## Grown Drug-Like Candidates

{candidates_section}

---

## Interpretation

IGHMBP2 is an SF1 helicase with no published small molecule inhibitors. This virtual screen
provides **first-in-class computational hits** predicted to bind the ATP/helicase pocket.

**BBB note**: ADMET threshold relaxed to 0.3 for BBB penetration (standard 0.5) given
the neurological context of SMARD1 and expected CNS involvement.

---

## Next Steps

```bash
# Run full docking on top candidates
pipeline dock --target {gene_symbol} --all-scaffolds

# Run validation gates
pipeline validate --target {gene_symbol} --all-scaffolds --gate admet
pipeline validate --target {gene_symbol} --dashboard
```

**Experimental validation**: Top candidates with SA score < 3 are synthetically accessible
and can be ordered from Enamine REAL (check ZINC ID for commercial availability).
Recommended first assay: thermal shift (ΔTm) against purified IGHMBP2 helicase domain.
"""
    path = results_dir / "fragment_screen_report.md"
    path.write_text(report)
    return path


def run_fragment(
    gene_symbol: str,
    step: Optional[str],
    settings: Settings,
    cache: CacheManager,
    force: bool,
    top_n: int,
    exhaustiveness: int,
    library_size: int,
    console: Console,
) -> int:
    from pipeline.stages.fragment.pocket import run_pocket
    from pipeline.stages.fragment.library import run_library
    from pipeline.stages.fragment.screen import run_screen
    from pipeline.stages.fragment.cluster import run_cluster
    from pipeline.stages.fragment.grow import run_grow
    from pipeline.stages.fragment.output import run_output

    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)

    state = _load_state(gene_symbol, settings)
    if force:
        state = FragmentState(gene_symbol=gene_symbol)
        _save_state(state, settings)

    steps_to_run = [step] if step else list(_STEPS)

    pocket = None
    hits_df = None
    clusters_df = None
    candidates_df = None

    # Load existing intermediate results for report writing
    pocket_path = results_dir / "pocket_analysis.json"
    if pocket_path.exists():
        try:
            pocket = json.loads(pocket_path.read_text())
        except Exception:
            pass
    hits_path = results_dir / "fragment_hits.csv"
    if hits_path.exists():
        try:
            hits_df = pd.read_csv(hits_path)
        except Exception:
            pass
    clusters_path = results_dir / "fragment_clusters.csv"
    if clusters_path.exists():
        try:
            clusters_df = pd.read_csv(clusters_path)
        except Exception:
            pass
    candidates_path = results_dir / "grown_candidates.csv"
    if candidates_path.exists():
        try:
            candidates_df = pd.read_csv(candidates_path)
        except Exception:
            pass

    for step_name in steps_to_run:
        step_attr = f"step_{step_name}"

        if not force and getattr(state, step_attr, False):
            console.print(f"  [dim]{gene_symbol}:[/dim] {step_name} [yellow]SKIP[/yellow] (cached)")
            continue

        console.print(f"  [dim]{gene_symbol}:[/dim] running step: [bold]{step_name}[/bold]...")

        try:
            if step_name == "pocket":
                pocket = run_pocket(gene_symbol, settings, cache, force, console)
                state.step_pocket = True

            elif step_name == "library":
                lib_path = run_library(library_size, settings, cache, force, console)
                state.step_library = True

            elif step_name == "dock":
                if pocket is None:
                    raise RuntimeError("Pocket not identified. Run `--step pocket` first.")
                lib_path = settings.cache_dir / "shared" / "fragment_library" / "fragments_ro3.smi"
                if not lib_path.exists():
                    raise RuntimeError("Fragment library not ready. Run `--step library` first.")
                hits_df = run_screen(gene_symbol, lib_path, pocket, top_n, exhaustiveness, settings, cache, force, console)
                state.step_dock = True
                state.n_fragments_docked = len(hits_df) if hits_df is not None else 0
                state.n_hits = min(top_n, state.n_fragments_docked)

            elif step_name == "cluster":
                if hits_df is None:
                    raise RuntimeError("No fragment hits found. Run `--step dock` first.")
                clusters_df = run_cluster(gene_symbol, hits_df, settings, cache, force, console)
                state.step_cluster = True
                state.n_clusters = clusters_df["cluster_id"].nunique() if clusters_df is not None and "cluster_id" in clusters_df.columns else 0

            elif step_name == "grow":
                if clusters_df is None:
                    raise RuntimeError("No clusters found. Run `--step cluster` first.")
                candidates_df = run_grow(gene_symbol, clusters_df, settings, cache, force, console)
                state.step_grow = True
                state.n_grown = len(candidates_df) if candidates_df is not None else 0

            elif step_name == "admet":
                if candidates_df is None:
                    raise RuntimeError("No grown candidates found. Run `--step grow` first.")
                n_final = run_output(gene_symbol, candidates_df, settings, cache, force, console)
                state.step_admet = True
                state.n_candidates_final = n_final

            _save_state(state, settings)

        except Exception as exc:
            console.print(f"  [red]{gene_symbol} fragment {step_name} FAIL: {exc}[/red]")
            _save_state(state, settings)
            return 1

    # Write report
    if not step or step == "admet":
        report_path = _write_fragment_report(gene_symbol, state, pocket, hits_df, clusters_df, candidates_df, settings)
        console.print(f"  [dim]{gene_symbol}:[/dim] [green]fragment_screen_report.md written[/green]")

    if state.step_admet:
        state.completed_at = datetime.now(timezone.utc).isoformat()
        _save_state(state, settings)
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [green bold]Fragment screen complete — "
            f"{state.n_candidates_final} candidates in compounds_filtered.csv[/green bold]"
        )
        console.print(f"  [dim]{gene_symbol}:[/dim] Next: `pipeline dock --target {gene_symbol} --all-scaffolds`")

    return 0
