# 010 · GEO differential expression analysis — volcano / heatmap / PCA

Two-group transcriptome differential expression with limma, plus volcano plot, PCA, and DEG clustering heatmap.

| | |
|---|---|
| **Language / main dependencies** | R · `limma` `ComplexHeatmap` `ggplot2` `ggrepel` |
| **Purpose** | Two-group transcriptome differential expression with the standard three-figure visualization set |
| **Input** | `example_data/expr_matrix.csv` (gene × sample) |
| **Output** | Tables and figures in `results/`; display figures in `assets/` |

## Input

**File**: `expr_matrix.csv` (CSV; rows = genes, columns = samples)

| Column | Type | Required | Example | Notes |
|------|------|:---:|------|------|
| Column 1 (gene name) | str | Yes | `TP53` | Set as row names |
| Sample columns ×N | num | Yes | `8.21` | Expression values, log2-normalized recommended |

**Naming / format convention**: sample column-name suffixes distinguish groups — control `*_con`, treatment `*_tre` (suffixes can be changed with `--ctrl/--case`). At least 3 replicates per group recommended.

**Example**:
```
Gene,S01_con,...,S01_tre,...
NUMA1,8.10,...,10.32,...
```

## Method

`limma`: `lmFit` linear model, `makeContrasts(Disease−Control)`, `eBayes` empirical Bayes, `topTable` (BH/FDR correction). Significant DEGs selected by `|log2FC| > threshold & FDR < threshold`. PCA uses `prcomp` (standardized). The heatmap uses `ComplexHeatmap` (row z-score standardization and two-way clustering).

Method citations: Ritchie *et al.*, *NAR* 2015 (limma); Gu *et al.*, *Bioinformatics* 2016 (ComplexHeatmap).

## Use case

Standard differential expression workflow for two-group comparisons in GEO/RNA-seq data (disease vs control, treated vs untreated). Produces the DEG lists needed for downstream enrichment (007), machine-learning features (04 modules), and diagnostic models (05 modules).

## Features

- Runs the example without edits; switch data with `--input`; group suffixes, thresholds, and label counts are configurable.
- Volcano plot with logFC gradient coloring, point size mapped to significance, italic labels for top genes, and up/down-regulated counts.
- Automatic delimiter detection; explicit messages when group suffixes are missing or no significant genes are found.
- Each figure exported as PDF and 300 dpi PNG.

## Outputs

Each figure is a separate file (PDF and PNG).

| File | Figure type | Notes |
|------|------|------|
| `assets/DEG_volcano.png` | Gradient volcano plot | Italic labels for top up/down-regulated genes |
| `assets/DEG_PCA.png` | PCA scatter | 95% confidence ellipses, between-group separation |
| `assets/DEG_heatmap.png` | Clustering heatmap | Top DEGs × samples, group annotation |
| `results/DE_results.csv` · `DE_significant_genes.csv` | Table | All / significant DEGs |

![Volcano plot](assets/DEG_volcano.png)
![Heatmap](assets/DEG_heatmap.png)

## Usage

```bash
Rscript 010_GEO_DEG_volcano_heatmap_PCA.R                          # 跑示例
Rscript 010_GEO_DEG_volcano_heatmap_PCA.R --input data/expr.csv --logfc 1 --padj 0.05 --topn 20
```

## Dependencies

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("limma","ComplexHeatmap"))
install.packages(c("ggplot2","ggrepel","circlize"))
```
