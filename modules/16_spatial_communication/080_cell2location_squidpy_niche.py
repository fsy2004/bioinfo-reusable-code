#!/usr/bin/env python
"""Spatial niche wrapper for cell2location outputs and Squidpy neighborhood analysis.

This wrapper intentionally keeps cell2location training explicit because model setup
depends on tissue, reference cell labels, and GPU availability. It can:
1. read spatial h5ad;
2. optionally attach a precomputed cell abundance table;
3. run Squidpy spatial graph and neighborhood enrichment;
4. save a standardized h5ad and CSV outputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--spatial_h5ad", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--cell_abundance_csv", help="Rows are spots, columns are cell types.")
    p.add_argument("--cluster_key", default="cluster")
    p.add_argument("--library_key", default=None)
    p.add_argument("--coord_type", default="visium", choices=["visium", "generic"])
    p.add_argument("--n_neighs", type=int, default=6)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import pandas as pd
    import scanpy as sc
    import squidpy as sq

    adata = sc.read_h5ad(args.spatial_h5ad)

    if args.cell_abundance_csv:
      abundance = pd.read_csv(args.cell_abundance_csv, index_col=0)
      abundance = abundance.reindex(adata.obs_names)
      for col in abundance.columns:
          adata.obs[f"abundance_{col}"] = abundance[col].values
      abundance.to_csv(outdir / "cell_abundance_reindexed.csv")

    sq.gr.spatial_neighbors(
        adata,
        coord_type=args.coord_type,
        n_neighs=args.n_neighs,
        library_key=args.library_key,
    )

    if args.cluster_key in adata.obs:
        sq.gr.nhood_enrichment(adata, cluster_key=args.cluster_key)
        z = adata.uns[f"{args.cluster_key}_nhood_enrichment"]["zscore"]
        clusters = adata.obs[args.cluster_key].astype("category").cat.categories
        pd.DataFrame(z, index=clusters, columns=clusters).to_csv(outdir / "squidpy_neighborhood_enrichment_zscore.csv")

    adata.write_h5ad(outdir / "spatial_niche_squidpy_result.h5ad")

    note = outdir / "cell2location_next_steps.md"
    note.write_text(
        "cell2location integration note\n\n"
        "Use the external source or installed package to train the regression model on scRNA reference, "
        "then export spot-level abundance as --cell_abundance_csv and rerun this wrapper for Squidpy "
        "neighborhood enrichment and downstream spatial niche analysis.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
