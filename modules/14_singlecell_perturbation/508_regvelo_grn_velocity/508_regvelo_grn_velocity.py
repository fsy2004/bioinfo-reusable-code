"""508 · RegVelo — gene-regulatory-informed RNA velocity & regulon perturbation.

Baseline (always runnable, CPU): scVelo velocity + CellRank fate probabilities, no GRN.
RegVelo path (--run-regvelo, GPU): couples splicing dynamics to a TF->target skeleton.

The baseline exists so the GRN-informed model is never reported without a comparator.
RegVelo: Wang et al., Cell 2026, doi:10.1016/j.cell.2026.04.022 (PMID 42119563).
Repo https://github.com/theislab/regvelo · Docs https://regvelo.readthedocs.io
"""
from __future__ import annotations
import argparse, os, sys, warnings
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")


def _synthetic(n_cells: int = 600, n_genes: int = 120, seed: int = 0):
    """Small spliced/unspliced AnnData along a synthetic differentiation axis."""
    import numpy as np, anndata as ad, pandas as pd
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 1, n_cells))                    # latent progression
    base = rng.uniform(0.5, 3.0, n_genes)
    # genes ramp up or down along t; unspliced leads spliced (that is what velocity reads)
    direction = rng.choice([1, -1], n_genes)
    prof = base[None, :] * (1 + direction[None, :] * t[:, None])
    spliced = rng.poisson(np.clip(prof, 0.05, None))
    unspliced = rng.poisson(np.clip(prof * (1 + 0.35 * direction[None, :]), 0.05, None))
    A = ad.AnnData(spliced.astype("float32"))
    A.layers["spliced"] = spliced.astype("float32")
    A.layers["unspliced"] = unspliced.astype("float32")
    A.var_names = [f"G{i:03d}" for i in range(n_genes)]
    A.obs_names = [f"C{i:04d}" for i in range(n_cells)]
    A.obs["true_time"] = t
    return A


def baseline(adata, outdir: str):
    """No-GRN floor: scVelo velocity -> CellRank fates. Reported alongside RegVelo."""
    import scanpy as sc, scvelo as scv, numpy as np, pandas as pd
    # n_top_genes is deliberately not passed: some scvelo versions forward it into
    # normalize_per_cell, which does not accept it. Keep all genes for the demo.
    scv.pp.filter_and_normalize(adata, min_shared_counts=0)
    scv.pp.moments(adata, n_pcs=20, n_neighbors=15)
    # Mode fallback: stochastic needs well-behaved second moments and fails on
    # near-degenerate data (leastsq_generalized raises). Deterministic is the safe floor.
    mode_used = None
    for mode in ("stochastic", "deterministic"):
        try:
            scv.tl.velocity(adata, mode=mode)
            mode_used = mode
            break
        except Exception as e:
            print(f"       velocity mode='{mode}' failed ({type(e).__name__}), trying next")
    if mode_used is None:
        return {"error": "scVelo could not fit a velocity field on this input"}
    scv.tl.velocity_graph(adata)
    sc.tl.umap(adata)

    out = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    try:                                                # CellRank is optional for the floor
        import cellrank as cr
        vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
        g = cr.estimators.GPCCA(vk)
        g.compute_schur(n_components=4)
        g.compute_macrostates(n_states=2)
        g.set_terminal_states_from_macrostates()
        g.compute_fate_probabilities()
        fp = g.fate_probabilities
        pd.DataFrame(np.asarray(fp), index=adata.obs_names,
                     columns=[str(c) for c in fp.names]).to_csv(
            os.path.join(outdir, "baseline_fate_probabilities.csv"))
        out["fates"] = [str(c) for c in fp.names]
    except Exception as e:                              # never fail the floor on an optional dep
        out["cellrank"] = f"skipped: {type(e).__name__}"

    # velocity-vs-latent-time consistency: does the no-GRN field track the known axis?
    if "velocity_pseudotime" not in adata.obs:
        try:
            scv.tl.velocity_pseudotime(adata)
        except Exception:
            pass
    if "velocity_pseudotime" in adata.obs and "true_time" in adata.obs:
        from scipy.stats import spearmanr
        r = spearmanr(adata.obs["velocity_pseudotime"], adata.obs["true_time"]).correlation
        out["baseline_pseudotime_vs_truth_rho"] = round(float(r), 3)
    return out


def run_regvelo(adata, skeleton=None):
    """RegVelo path. Guarded: needs `pip install regvelo` and ideally a GPU.

    Exported API (verified from the package): REGVELOVI, VELOVAE, ModelComparison,
    and pp / tl / pl / mt submodules. scvi-tools pattern applies
    (setup_anndata -> construct -> train). Confirm exact signatures against the
    official tutorial before a production run; they are deliberately not pinned here.
    """
    try:
        import regvelo as rgv
        import torch
    except ImportError as e:
        return {"status": "skipped", "reason": f"regvelo not installed ({e.name}); pip install regvelo"}
    gpu = bool(getattr(__import__("torch"), "cuda", None) and __import__("torch").cuda.is_available())
    if not gpu:
        return {"status": "skipped", "reason": "no CUDA GPU; RegVelo training is impractical on CPU"}

    exported = [n for n in ("REGVELOVI", "VELOVAE", "ModelComparison") if hasattr(rgv, n)]
    return {
        "status": "ready",
        "regvelo_version": getattr(rgv, "__version__", "?"),
        "exported_api": exported,
        "next": ("follow https://regvelo.readthedocs.io tutorials: "
                 "REGVELOVI.setup_anndata(adata, spliced_layer=..., unspliced_layer=...) -> "
                 "REGVELOVI(adata, skeleton=<TF x target prior>) -> .train() -> "
                 "regulon perturbation via the CellRank framework"),
        "skeleton_supplied": skeleton is not None,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", help="AnnData with spliced/unspliced layers; omit for synthetic demo")
    ap.add_argument("--run-regvelo", action="store_true", help="attempt the RegVelo path (needs install + GPU)")
    ap.add_argument("--outdir", default=RESULTS)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True); os.makedirs(ASSETS, exist_ok=True)

    if a.h5ad:
        import anndata as ad
        adata = ad.read_h5ad(a.h5ad)
        for L in ("spliced", "unspliced"):
            if L not in adata.layers:
                sys.exit(f"input lacks layers['{L}'] — RegVelo and scVelo both need spliced/unspliced")
    else:
        print("[508] no --h5ad given, using synthetic spliced/unspliced demo")
        adata = _synthetic()

    print("[508] baseline (no GRN): scVelo velocity + CellRank fates")
    b = baseline(adata, a.outdir)
    for k, v in b.items():
        print(f"       {k}: {v}")

    if a.run_regvelo:
        print("[508] RegVelo path")
        r = run_regvelo(adata)
        for k, v in r.items():
            print(f"       {k}: {v}")
    else:
        print("[508] RegVelo path not requested (--run-regvelo); baseline only")

    import json
    with open(os.path.join(a.outdir, "508_summary.json"), "w") as fh:
        json.dump({"baseline": b}, fh, indent=1, default=str)
    print(f"[508] wrote {os.path.join(a.outdir, '508_summary.json')}")


if __name__ == "__main__":
    main()
