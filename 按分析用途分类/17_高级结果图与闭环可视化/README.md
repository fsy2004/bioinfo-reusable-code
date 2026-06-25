# 17 · Advanced result figures and closed-loop visualization

Reusable tool entry points and template concepts for organizing multi-omics evidence into a closed loop rather than restyling bar charts.

## Modules

| Module | Purpose | Language | Output |
|------|------|------|--------|
| [498 ggalluvial Sankey/alluvial](498_ggalluvial桑基冲积图/) | Multi-layer alluvial flow (drug, hub, pathway, etc.) | R | Alluvial/Sankey plot |

498 follows the [unified framework conventions](../_framework/CONVENTIONS.md) and runs without modification. The rest of this module documents the design concepts and tool index for advanced figures in this category.

## Recommended closed loop

`phenotype -> cell state -> communication / spatial niche -> perturbable target -> validation / prediction`

## Priority figure types

1. Patient-level risk landscape
   Places machine learning predicted probability, true labels, dataset source, cell-state proportions, and SHAP contributions in a single patient-level heatmap.

2. Response circuit circos
   Uses circular tracks to show cell-state shift, SHAP, and Scissor enrichment, with chord links for communication rewiring.

3. Target evidence wheel
   Connects virtual perturbation target, cell-state specificity, drug support, and external validation into a therapeutic evidence figure.

4. Ligand-target-response circuit
   Links sender cell, ligand, receiver cell, target gene program, and response phenotype.

5. Spatial niche response map
   Maps single-cell states onto spatial sections and overlays niche, communication vector, or pathway activity.

## Files

- `advanced_figure_tools.csv`: recommended tools, purpose, paper, and GitHub address.
- `download_advanced_figure_tools.ps1`: clone/update the relevant tools to the local `external_tools` directory.
- `literature_download_links_for_fdm.txt`: batch download links for FDM.
- `templates/`: portable template notes.

## Usage notes

- Do not commit third-party tool source code into this repository; pull it with `download_advanced_figure_tools.ps1`.
- For a new project, copy only the template scripts you need and keep this module as a tool index.
- Advanced figures should serve the evidence closed loop; replacing a bar chart with a circular chart is usually not worth a main figure.
