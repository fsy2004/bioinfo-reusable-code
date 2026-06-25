# 014 · Random forest feature gene selection

Ranks candidate genes by random forest Gini importance and produces an OOB error rate plot and an importance lollipop plot.

| | |
|---|---|
| **Language / main dependencies** | R · `randomForest` `ggplot2` |
| **Purpose** | Rank and select candidate genes by Gini importance |
| **Input** | `example_data/Sample_Type_Matrix.csv` + `candidate_genes.csv` |
| **Output** | Importance table and plots in `results/`; display figures in `assets/` |

## Input

Same as [012](../012_LASSO特征基因筛选/): `--input` expression matrix (group encoded in the sample name suffix) plus `--genes` candidate genes (optional).

## Method

`randomForest` (500 trees by default) is fit with the group as the response, then `importance()` returns `MeanDecreaseGini` for ranking, and features are selected by threshold or topN. The OOB error rate curve helps assess whether the number of trees is sufficient.

Method citation: Breiman, *Machine Learning* 2001 (Random Forests).

## Use

Nonlinear, collinearity-robust feature importance assessment, often intersected with LASSO/SVM-RFE (see 015) to obtain a robust feature set.

## Notes

- Runs on the example data without modification; `--ntree/--top/--threshold` are configurable.
- OOB error rate curves for multiple classes and a viridis importance lollipop (gene names in italics).

## Outputs

| File | Plot type | Description |
|------|------|------|
| `assets/RF_importance_lollipop.png` | Lollipop plot | Gini importance of top genes |
| `assets/RF_OOB_error.png` | Line plot | OOB and per-class error rate vs number of trees |
| `results/RF_gene_importance.csv` | Table | Importance for all genes |

![importance](assets/RF_importance_lollipop.png)

## Usage

```bash
Rscript 014_RandomForest_feature_selection.R                                  # 示例
Rscript 014_RandomForest_feature_selection.R --input data/expr.csv --top 20 --ntree 1000
```

## Dependencies

```r
install.packages(c("randomForest","ggplot2","reshape2","viridisLite"))
```
