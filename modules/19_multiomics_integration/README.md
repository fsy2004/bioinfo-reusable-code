# 19 · Multi-omics Integration / Subtyping Templates

Multi-omics latent-variable integration, supervised multi-omics classification, and unsupervised molecular subtyping.

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [084 NMF consensus subtyping](084_nmf_consensus_clustering/) | NMF + consensus clustering molecular subtyping | R | consensus matrix, rank curve, subtype heatmap | Ready |
| 083 MOFA / DIABLO multi-omics integration | Multi-omics latent-variable integration | R | factor plot, heatmap | Requires MOFA2 (python) |

## Notes

- 084: produces standard subtyping figures directly from a feature matrix (follows the [unified framework conventions](../_framework/CONVENTIONS.md)).
- 083: MOFA2 depends on the python package `mofapy2` (via reticulate); DIABLO depends on mixOmics. Not rendered locally; the original scripts are kept for reference.

## Input

- Matrix: rows = features, columns = samples (expression / immune scores / pathway activity / spatial niches).
- metadata: sample groups, disease status, clinical phenotypes.
