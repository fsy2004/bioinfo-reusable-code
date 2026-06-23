# 03 · GEO 转录组整理与差异分析

GEO 微阵列/转录组从**原始下载 → 整理归一化 → 差异分析**的完整 turnkey 链路。每个模块自带示例数据,`Rscript <脚本>` 一条命令即可跑通。

## 模块链路

```
008 探针→基因矩阵 ──▶ 009 归一化+分组 ──▶ 010 差异分析(火山/热图/PCA)
056 多队列合并+批次校正(独立模块,可在 010 前用于多队列整合)
```

| 模块 | 用途 | 语言 | 输出图 |
|------|------|------|--------|
| [008 表达矩阵整理](008_GEO表达矩阵整理/) | GSE+GPL 探针→基因级矩阵 | R | —(产矩阵) |
| [009 样本分组整理](009_GEO样本分组整理/) | 归一化 + 写入分组后缀 | R | —(产矩阵) |
| [010 差异分析 火山/热图/PCA](010_GEO差异分析_火山热图PCA/) | limma DEG + 三件套可视化 | R | 渐变火山图 · PCA · 聚类热图 |
| [056 多队列合并+批次校正](056_GEO多队列合并_批次校正/) | 合并多 GEO + 去批次 | R | 校正前后 PCA · 箱线图 |

## 典型用法(完整流程)

```bash
# 1) 探针→基因
Rscript 008_GEO表达矩阵整理/008_GEO_expr_matrix_tidy.R --gse GSExxx_series_matrix.txt --gpl GPLxxx.txt --symcol 11
# 2) 归一化+分组(产出 Sample_Type_Matrix.csv)
Rscript 009_GEO样本分组整理/009_GEO_sample_grouping.R --expr results/geneMatrix.csv --group group.csv
# 3) 差异分析+出图
Rscript 010_GEO差异分析_火山热图PCA/010_GEO_DEG_volcano_heatmap_PCA.R --input results/Sample_Type_Matrix.csv
```

> 全部遵循 [统一框架规范](../_framework/CONVENTIONS.md):去硬编码路径、`--input/--outdir`、共享顶刊主题、独立矢量图。
