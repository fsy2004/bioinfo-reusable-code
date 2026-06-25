#!/usr/bin/env python
"""Palantir trajectory wrapper for h5ad input."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--root_cell", required=True)
    p.add_argument("--n_components", type=int, default=10)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import palantir
    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(args.h5ad)
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata)

    dm_res = palantir.utils.run_diffusion_maps(adata, n_components=args.n_components)
    ms_data = palantir.utils.determine_multiscale_space(adata)
    pr_res = palantir.core.run_palantir(ms_data, args.root_cell)

    pd.DataFrame(pr_res.pseudotime).to_csv(outdir / "palantir_pseudotime.csv")
    pd.DataFrame(pr_res.branch_probs).to_csv(outdir / "palantir_branch_probabilities.csv")
    pd.DataFrame(pr_res.entropy).to_csv(outdir / "palantir_entropy.csv")
    adata.write_h5ad(outdir / "palantir_result.h5ad")


if __name__ == "__main__":
    main()
