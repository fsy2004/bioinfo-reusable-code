# 14 · Single-cell virtual perturbation and perturbation databases

Methods for moving single-cell analysis from describing cell states to predicting how cells respond to perturbation.

## Status

This category covers single-cell virtual perturbation methods. Most depend on deep learning or graph neural networks (GEARS, CellOracle, Squidiff, GenKI), which require a GPU and trained models, or on dedicated perturbation databases. No example figures are rendered locally; the original scripts are kept for reference, with dependencies and reproducible commands in each script header. `067_scPerturb` (R) and `495_bulkVGK` (scTenifoldKnk) are lighter and can be run when suitable data is available. See [unified framework](../_framework/CONVENTIONS.md) for figure conventions.

## Scripts

| Script | Purpose |
|---|---|
| `067_scperturb_etest.R` | Compute perturbation distance and E-test between perturbation and control groups. |
| `068_gears_combo_perturbation.py` | Run GEARS to predict expression response to single-gene or multi-gene combinatorial perturbation. |
| `069_celloracle_grn_perturbation.py` | GRN virtual knockdown/knockout using a CellOracle Oracle object. |
| `085_squidiff_diffusion_perturbation.py` | Run external Squidiff/PerturbDiff training or sampling scripts, recording reproducible commands. |

## Input

- Seurat RDS or AnnData h5ad.
- Perturbation/control grouping column.
- Candidate gene list.
- Trained model or GRN/Oracle object.

## Outputs

- Perturbation distance and E-test results.
- GEARS predicted expression matrix and run summary.
- CellOracle post-perturbation state transition results.
- Squidiff/PerturbDiff prediction results and run logs.
