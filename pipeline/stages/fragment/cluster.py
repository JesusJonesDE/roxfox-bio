from __future__ import annotations

import pandas as pd
from rich.console import Console

from pipeline.cache import CacheManager
from pipeline.config import Settings


def run_cluster(
    gene_symbol: str,
    hits_df: pd.DataFrame,
    settings: Settings,
    cache: CacheManager,
    force: bool,
    console: Console,
) -> pd.DataFrame:
    """Cluster fragment hits by Morgan ECFP4 / Butina; return annotated DataFrame."""
    # 1. Cache check (before any heavy imports)
    cache_key = "fragment_cluster"
    if not force:
        cached = cache.load(gene_symbol, cache_key)
        if cached is not None:
            console.print(
                f"  [dim]{gene_symbol}:[/dim] cluster [yellow]SKIP[/yellow] (cached)"
            )
            return pd.DataFrame(cached)

    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.ML.Cluster import Butina

    # 2. Compute Morgan ECFP4 fingerprints
    fps = []
    valid_idx = []
    for i, row in hits_df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol:
            fps.append(
                AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
            )
            valid_idx.append(i)

    n = len(fps)
    if n == 0:
        console.print(
            f"  [dim]{gene_symbol}:[/dim] [yellow]no valid SMILES in hits_df — "
            f"returning empty clusters[/yellow]"
        )
        return hits_df.copy()

    # 3. Compute distance matrix (lower-triangular)
    dists: list[float] = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1.0 - s for s in sims])

    # 4. Butina clustering
    clusters = Butina.ClusterData(dists, n, 0.4, isDistData=True)
    # clusters is a tuple of tuples of indices into fps/valid_idx

    # 5. Build cluster_id assignment and attach to a working copy of hits_df
    cluster_id_by_pos: dict[int, int] = {}
    for cluster_idx, member_positions in enumerate(clusters):
        for pos in member_positions:
            cluster_id_by_pos[pos] = cluster_idx

    # 6. Attach cluster_id to a copy of hits_df using valid_idx alignment
    out_df = hits_df.copy()
    # valid_idx[pos] is the DataFrame index for fps[pos]
    cluster_id_col: dict = {
        valid_idx[pos]: cid for pos, cid in cluster_id_by_pos.items()
    }
    out_df["cluster_id"] = out_df.index.map(cluster_id_col)

    # Mark is_representative: the fragment with the most negative affinity per cluster
    # (hits_df arrives sorted ascending — most negative first — so the minimum affinity
    # is the best hit in each cluster regardless of Butina's centroid selection)
    out_df["is_representative"] = False
    if "affinity_kcal_mol" in out_df.columns:
        rep_idx = out_df.groupby("cluster_id")["affinity_kcal_mol"].idxmin()
        out_df.loc[rep_idx.values, "is_representative"] = True
    else:
        # Fallback: first member per cluster
        first_per_cluster = out_df.groupby("cluster_id").apply(lambda g: g.index[0])
        out_df.loc[first_per_cluster.values, "is_representative"] = True

    # 7. Write fragment_clusters.csv
    results_dir = settings.results_dir / gene_symbol
    results_dir.mkdir(parents=True, exist_ok=True)
    cols = ["fragment_id", "smiles", "affinity_kcal_mol", "cluster_id", "is_representative"]
    out_cols = [c for c in cols if c in out_df.columns]
    out_df[out_cols].to_csv(results_dir / "fragment_clusters.csv", index=False)

    n_clusters = out_df["cluster_id"].nunique()
    console.print(
        f"  [dim]{gene_symbol}:[/dim] [green]clustering complete — "
        f"{n_clusters} clusters from {n} valid fragments[/green]"
    )

    # 8. Cache result and return
    cache.save(gene_symbol, cache_key, out_df.to_dict(orient="records"), n_clusters)
    return out_df
