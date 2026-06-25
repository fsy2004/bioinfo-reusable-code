"""CellRank fate probability and driver gene wrapper."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run scVelo + CellRank fate mapping on AnnData.")
    parser.add_argument("--adata", required=True, help="Input h5ad with spliced/unspliced layers for RNA velocity.")
    parser.add_argument("--output-dir", default="results/cellrank")
    parser.add_argument("--n-states", type=int, default=6)
    parser.add_argument("--n-top-drivers", type=int, default=200)
    args = parser.parse_args()

    try:
        import scanpy as sc
        import scvelo as scv
        import cellrank as cr
    except ImportError as exc:
        raise SystemExit("Install dependencies first: pip install scanpy scvelo cellrank") from exc

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise SystemExit("CellRank velocity mode requires adata.layers['spliced'] and ['unspliced'].")

    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)

    kernel = cr.kernels.VelocityKernel(adata)
    kernel.compute_transition_matrix()
    estimator = cr.estimators.GPCCA(kernel)
    estimator.compute_macrostates(n_states=args.n_states)
    estimator.compute_terminal_states()
    estimator.compute_fate_probabilities()
    drivers = estimator.compute_lineage_drivers(return_drivers=True)

    adata.write_h5ad(out / "cellrank_fate_result.h5ad")
    if drivers is not None:
        drivers.head(args.n_top_drivers).to_csv(out / "cellrank_top_driver_genes.csv")

    print(f"Done. CellRank results saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

