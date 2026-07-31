"""
Robust GRN inference via GRNBoost2 + bootstrap.

Principle:
  1. Subsample the samples (rows) N times (e.g. 80% of samples, without replacement)
  2. Run GRNBoost2 on each subsample
  3. For each edge (TF -> target gene), compute:
       - its frequency of appearance across the N runs (stability)
       - its mean importance score
  4. Filter: keep stable edges (frequency >= threshold), sorted by
     decreasing mean score.

Input format:
    --expr : plain text/TSV file WITH a header row (sample names).
             Rows = genes, columns = samples.
             One column is named "gene_name" and holds the gene identifiers.

Usage:
    python grnboost2_bootstrap.py \
        --expr matrix_log2cpm.txt \
        --tf TF_list.txt \
        --out grn_final.tsv \
        --n_boot 100 \
        --frac 0.8 \
        --stability_min 0.7
"""

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from distributed import LocalCluster, Client


def parse_args():
    p = argparse.ArgumentParser(description="GRNBoost2 with bootstrap for a robust network")
    p.add_argument("--expr", required=True,
                   help="Expression matrix: genes as rows (column 'gene_name'), samples as columns")
    p.add_argument("--gene_col", default="gene_name",
                   help="Name of the column holding gene identifiers (default: gene_name)")
    p.add_argument("--tf", required=True, help="TF list file (one per line)")
    p.add_argument("--out", default="grn_final.tsv", help="Output file for the final network")
    p.add_argument("--sep", default="\t", help="Separator for the expression file (default: tab)")
    p.add_argument("--n_boot", type=int, default=100, help="Number of bootstrap iterations")
    p.add_argument("--frac", type=float, default=0.8,
                   help="Fraction of samples kept at each subsampling step (0-1)")
    p.add_argument("--stability_min", type=float, default=0.7,
                   help="Minimum frequency of appearance of an edge across the N bootstraps (0-1)")
    p.add_argument("--seed", type=int, default=42, help="Seed for reproducibility")
    p.add_argument("--n_workers", type=int, default=4, help="Number of Dask workers (local computation)")
    return p.parse_args()


def run_bootstrap(ex_matrix, tf_names, n_boot, frac, seed, client):
    """Run GRNBoost2 on n_boot subsamples and aggregate the results."""
    rng = np.random.default_rng(seed)
    n_samples = ex_matrix.shape[0]
    n_sub = int(n_samples * frac)

    # edge_key -> list of importance scores observed across the runs where it appears
    edge_scores = defaultdict(list)

    for i in range(n_boot):
        idx = rng.choice(n_samples, size=n_sub, replace=False)
        sub_matrix = ex_matrix.iloc[idx, :]

        print(f"[Bootstrap {i + 1}/{n_boot}] {n_sub} samples out of {n_samples}")

        net = grnboost2(
            expression_data=sub_matrix,
            tf_names=tf_names,
            client_or_address=client,
            seed=int(rng.integers(0, 1_000_000)),
        )
        # net has columns: TF, target, importance
        for tf, target, importance in net.itertuples(index=False):
            edge_scores[(tf, target)].append(importance)

    return edge_scores


def build_final_network(edge_scores, n_boot, stability_min):
    """Filter by stability, sorted by mean score."""
    rows = []
    for (tf, target), scores in edge_scores.items():
        frequency = len(scores) / n_boot
        if frequency >= stability_min:
            rows.append({
                "TF": tf,
                "target": target,
                "mean_importance": np.mean(scores),
                "std_importance": np.std(scores),
                "frequency": frequency,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError(
            "No edge passes the stability threshold. "
            "Lower --stability_min or increase --n_boot."
        )

    df = df.sort_values("mean_importance", ascending=False).reset_index(drop=True)
    return df


def main():
    args = parse_args()

    print("Loading expression matrix...")
    raw_matrix = pd.read_csv(args.expr, sep=args.sep, header=None, index_col=None)

    with open(args.genes) as f:
        gene_names = [line.strip() for line in f if line.strip()]

    if len(gene_names) != raw_matrix.shape[0]:
        raise ValueError(...)

    raw_matrix.index = gene_names

    ex_matrix = raw_matrix.T
    ex_matrix.columns = gene_names

    print(f"Matrix: {ex_matrix.shape[0]} samples x {ex_matrix.shape[1]} genes")

    tf_names = load_tf_names(args.tf)
    # Keep only TFs present in the matrix
    tf_names = [tf for tf in tf_names if tf in ex_matrix.columns]
    print(f"{len(tf_names)} TFs found in the expression matrix")

    print(f"Starting local Dask cluster ({args.n_workers} workers)...")
    cluster = LocalCluster(n_workers=args.n_workers, threads_per_worker=1)
    client = Client(cluster)

    try:
        edge_scores = run_bootstrap(
            ex_matrix, tf_names,
            n_boot=args.n_boot,
            frac=args.frac,
            seed=args.seed,
            client=client,
        )
    finally:
        client.close()
        cluster.close()

    print(f"\n{len(edge_scores)} unique edges observed across {args.n_boot} bootstraps")

    final_net = build_final_network(
        edge_scores,
        n_boot=args.n_boot,
        stability_min=args.stability_min,
    )

    print(f"Final network: {len(final_net)} edges "
          f"(stability >= {args.stability_min})")

    final_net.to_csv(args.out, sep="\t", index=False)
    print(f"Saved to {args.out}")


if __name__ == "__main__":
    main()