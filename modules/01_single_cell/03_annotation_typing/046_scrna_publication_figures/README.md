# 046 · Single-cell publication figures (standard Seurat workflow)

Runs the standard Seurat workflow on a count matrix and produces publication figures: UMAP, marker dot plot, marker heatmap, and target-gene plots.

| | |
|---|---|
| Language / main dependencies | R · `Seurat` `ggplot2` |
| Purpose | Standard single-cell analysis with publication figures |
| Input | `example_data/counts.csv` |
| Output | `results/` (object + markers + figures); preview figures in `assets/` |

## Input

Count matrix CSV: first column is gene names, remaining columns are cells (raw counts). Target genes to plot can be specified with `--genes "CD3D,MS4A1"`.

## Method

Standard Seurat workflow: QC, `NormalizeData`, HVG, `ScaleData`, `RunPCA`, `FindClusters`, `RunUMAP`, `FindAllMarkers`, then plot the top marker per cluster.

Method citation: Hao *et al.*, *Cell* 2021 (Seurat v4/v5).

## Use

Clustering, marker identification, and figure generation for any single-cell or single-nucleus RNA-seq dataset.

## Notes

- Runs the full workflow directly from counts; falls back to tSNE if UMAP fails.
- Figures: UMAP (discrete journal palette), marker dot plot, marker heatmap (viridis), target-gene FeaturePlot (viridis), and violin plot. The default Set3 and grey-red palettes have been replaced.

## Outputs

| File | Type | Description |
|------|------|------|
| `assets/UMAP_clusters.png` | UMAP | Clusters |
| `assets/Marker_dotplot.png` | Dot plot | Top marker per cluster (block-diagonal) |
| `assets/Marker_heatmap.png` | Heatmap | Top marker expression |
| `assets/<gene>_FeaturePlot.png` · `<gene>_violin.png` | FeaturePlot / violin | Target-gene distribution |

![dotplot](assets/Marker_dotplot.png)
![umap](assets/UMAP_clusters.png)

## Usage

```bash
Rscript 046_scRNA_publication_figures.R                                   # 示例
Rscript 046_scRNA_publication_figures.R --input data/counts.csv --genes "CD3D,MS4A1" --resolution 0.6
```

## Dependencies

```r
install.packages(c("Seurat","ggplot2","dplyr"))
```
