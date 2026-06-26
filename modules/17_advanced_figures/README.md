# 17 · Advanced result figures and closed-loop visualization

Reusable tool entry points and template concepts for organizing multi-omics evidence into a closed loop rather than restyling bar charts.

## Modules

| Module | Purpose | Language | Output |
|------|------|------|--------|
| [498 ggalluvial Sankey/alluvial](498_ggalluvial_sankey/) | Multi-layer alluvial flow (drug, hub, pathway, etc.) | R | Alluvial/Sankey plot |
| [512 raincloud](512_raincloud_plot/) | Distribution comparison that replaces grouped bar charts | R | Raincloud (half-violin + box + jitter) |
| [513 ridgeline](513_ridgeline_plot/) | Distribution shift across an ordered factor (time/stage/pseudotime) | R | Ridgeline / joyplot |
| [514 dumbbell + slope](514_dumbbell_slope_plot/) | Paired change across two conditions | R | Dumbbell, slopegraph |
| [515 chord](515_chord_diagram/) | Directed relations / flows (e.g. communication strength) | R | Chord diagram |
| [516 composite multi-panel](516_composite_multipanel/) | "Figure 1" template (UMAP + volcano + heatmap + forest) | R | Composite multi-panel figure |

All six modules follow the [unified framework conventions](../_framework/CONVENTIONS.md) and run on bundled/synthetic data without modification. The sections below document the closed-loop design concepts and the external-tool index for this category. These figure types favour lollipop/dot/violin/raincloud/chord over plain bar charts.

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
