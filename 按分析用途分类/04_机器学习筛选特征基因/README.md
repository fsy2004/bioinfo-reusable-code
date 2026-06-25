# 04 · Machine learning feature gene selection

Select robust feature/marker genes from candidate genes (differential genes, network pharmacology targets, etc.) using machine learning, for diagnostic and prognostic modeling.

## Modules

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [012 LASSO](012_LASSO特征基因筛选/) | L1-regularized shrinkage selection | R | CV curve, coefficient path | Ready |
| [014 RandomForest](014_RandomForest特征基因筛选/) | Gini importance selection | R | OOB error rate, importance lollipop | Ready |
| [034 Multiple ML comparison](034_12种机器学习特征基因筛选/) | 10 ML methods + feature intersection | R | ROC overlay, AUC leaderboard, UpSet | Ready |
| [013 SVM-RFE](013_SVM_RFE特征基因筛选/) | Recursive feature elimination | R | CV accuracy curve, ranking | Ready |
| [015 ML feature intersection](015_机器学习特征交集_Venn_UpSet/) | Multi-method Venn / UpSet | R | Venn, UpSet | Ready |
| [035 Multi-method combination intersection](035_5种机器学习组合_交集选择/) | Combination intersection selection | R | Combination ranking, UpSet | Ready |
| [052 SHAP interpretation](052_SHAP机器学习解释分析/) | SHAP multi-plot interpretation | R | SHAP beeswarm/waterfall/force, ROC | Ready |
| 045 Multi-ML integrated signature | RSF/BART/... integration | R | ROC, heatmap | Heavy environment |
| 059 Dual-disease 15ML×175 combinations | Large-scale combination comparison | R | AUC heatmap, ROC | Heavy environment |
| 496 Mime 101-combination prognostic signature | Mime framework | R | AUC heatmap, KM | Heavy environment |

> Heavy environment: 045/059/496 are Mime-style 100+ combination prognostic signatures, depending on `randomForestSRC/BART/mboost/plsRglm/Mime` plus external helpers. They require a large ML environment, are not rendered locally, and the original scripts are kept for reference.

> Shared ML example data: `Sample_Type_Matrix.csv` (expression matrix, sample names `*_con`/`*_tre`) plus `candidate_genes.csv` (candidate genes).
> All modules follow the [unified framework conventions](../_framework/CONVENTIONS.md). Ready modules run without modification; modules pending refactoring may require additional packages.

## Typical workflow

```
03 differential genes -> 012/013/014 per-method selection -> 015/035 intersection -> 016/063 diagnostic model
                034 multi-method comparison (ROC/AUC ranking)         052 SHAP interpretation
```
