# 021 · Immune infiltration visualization

Generates grouped difference boxplots, stacked composition plots, and a correlation heatmap from an immune cell proportion matrix (e.g. CIBERSORT).

## Input

CSV file (`example_data/CIBERSORT_results.csv`): the first column `Sample` holds sample names (group encoded by suffix `*_con`/`*_tre`); the remaining columns are immune cell proportions (row sums approximately 1). This is the standard output of deconvolution tools such as CIBERSORT or quanTIseq (produced by modules 017/018).

## Method

Wilcoxon test per group with significance annotation (`ggpubr::stat_compare_means`); stacked composition plot of per-sample cell proportions; Spearman correlation matrix heatmap of immune cell co-infiltration relationships.

## Outputs

| File | Plot type | Description |
|------|------|------|
| `assets/Immune_boxplot.png` | Grouped boxplot | Two-group comparison per cell type with significance |
| `assets/Immune_stackbar.png` | Stacked bar | Per-sample immune composition |
| `assets/Immune_correlation.png` | Correlation heatmap | Immune cell co-infiltration |

Tables and figures are written to `results/`; example figures are in `assets/`.

![boxplot](assets/Immune_boxplot.png)
![stack](assets/Immune_stackbar.png)

## Usage

```bash
Rscript 021_immune_visualization.R                              # 示例
Rscript 021_immune_visualization.R --input data/CIBERSORT_results.csv
```

## Dependencies

R, with `ggpubr`, `ComplexHeatmap`, `ggplot2`.

```r
install.packages(c("ggpubr","ggplot2","reshape2","circlize"))
BiocManager::install("ComplexHeatmap")
```
