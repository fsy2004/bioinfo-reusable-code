# 03 · GEO transcriptome preparation and differential analysis

A pipeline for GEO microarray/transcriptome data covering raw download, tidying and normalization, and differential analysis. Each module ships with example data and runs with a single `Rscript <script>` command.

## Module chain

```
008 probe->gene matrix ──▶ 009 normalization + grouping ──▶ 010 differential analysis (volcano/heatmap/PCA)
056 multi-cohort merge + batch correction (standalone module, can be used before 010 for multi-cohort integration)
```

| Module | Purpose | Language | Output figures |
|------|------|------|--------|
| [008 expression matrix tidying](008_GEO表达矩阵整理/) | GSE+GPL probe to gene-level matrix | R | none (produces matrix) |
| [009 sample grouping](009_GEO样本分组整理/) | normalization + writing group suffix | R | none (produces matrix) |
| [010 differential analysis volcano/heatmap/PCA](010_GEO差异分析_火山热图PCA/) | limma DEG + three-part visualization | R | gradient volcano plot, PCA, clustering heatmap |
| [056 multi-cohort merge + batch correction](056_GEO多队列合并_批次校正/) | merge multiple GEO datasets + remove batch effects | R | PCA before/after correction, boxplots |

## Usage (full workflow)

```bash
# 1) probe->gene
Rscript 008_GEO表达矩阵整理/008_GEO_expr_matrix_tidy.R --gse GSExxx_series_matrix.txt --gpl GPLxxx.txt --symcol 11
# 2) normalization + grouping (produces Sample_Type_Matrix.csv)
Rscript 009_GEO样本分组整理/009_GEO_sample_grouping.R --expr results/geneMatrix.csv --group group.csv
# 3) differential analysis + figures
Rscript 010_GEO差异分析_火山热图PCA/010_GEO_DEG_volcano_heatmap_PCA.R --input results/Sample_Type_Matrix.csv
```

All modules follow the [shared framework conventions](../_framework/CONVENTIONS.md): no hardcoded paths, `--input/--outdir`, a shared plotting theme, and standalone vector figures.
