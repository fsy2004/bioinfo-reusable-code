# 08 · Single-cell / Spatial transcriptomics / Cell trajectory

Largest category in the collection. Module [046](046_scrna_publication_figures/) runs a standard Seurat workflow and produces a standard set of figures; the remaining modules are heavier or specialized workflows (CellChat, monocle, Scissor, spatial, deep learning), kept as original scripts for reference.

## Standard workflow

| Module | Purpose | Language | Output figures |
|------|------|------|--------|
| [046 Single-cell publication figures](046_scrna_publication_figures/) | Full Seurat workflow + standard figures | R | UMAP, marker dot plot, marker heatmap, FeaturePlot, violin |

## Engine / data preprocessing

| Module | Function |
|------|------|
| 023 RDS object structure inspection · 024/025 single-cell QC tidying | Data loading / QC / tidying (upstream for 046 and others) |

## Heavy / specialized (kept for reference)

| Module | Method | Not rendered locally because |
|------|------|----------------|
| 026 Seurat workflow + scTenifoldKnk knockout | Virtual knockout | scTenifoldKnk is heavy |
| 049 Manual annotation + CellChat + monocle | Cell communication + trajectory | CellChat/monocle are heavy |
| 051 CellChat cell communication | Circle / chord / bubble plots | CellChat is heavy |
| 058 Scissor disease-associated cells | Phenotype-associated cells | Scissor + cohort |
| 044 AD single-cell + monocle | Pseudotime | monocle is heavy |
| 027 / 050 Spatial transcriptomics | Visium spatial analysis | Requires spatial data |
| 062 scTour · 082 Palantir/Slingshot | Trajectory / vector field | Python deep learning / heavy |
| 061 scFOCAL/CellOracle GUI preparation · 491 scTour environment | Perturbation entry / environment | External GUI / environment |

Module 046 follows the [unified framework conventions](../_framework/CONVENTIONS.md). To run a heavy module, see the dependency notes at the top of each script.
