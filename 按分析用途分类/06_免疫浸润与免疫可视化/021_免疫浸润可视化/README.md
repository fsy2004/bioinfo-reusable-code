# 021 · 免疫浸润可视化

> 免疫细胞比例矩阵(CIBERSORT 等) → 一条命令 → 分组差异箱线图 + 堆叠组成图 + 相关性热图。

| | |
|---|---|
| **语言 / 主依赖** | R · `ggpubr` `ComplexHeatmap` `ggplot2` |
| **一句话用途** | 免疫去卷积结果的三大标配可视化 |
| **输入** | `example_data/CIBERSORT_results.csv` |
| **输出** | `results/` 表+图 · 展示图见 `assets/` |

---

## ① 输入数据

CSV:首列 `Sample`(样本名,后缀分组 `*_con`/`*_tre`),其余列 = 各免疫细胞比例(行和≈1)。即 CIBERSORT/quanTIseq 等去卷积的标准输出(由 017-020 产生)。

## ② 方法 / 原理

按分组做 Wilcoxon 差异检验(`ggpubr::stat_compare_means`,标显著性)→ 堆叠组成图(各样本细胞比例)→ Spearman 相关矩阵热图(免疫细胞共浸润关系)。

## ③ 用途

展示疾病/处理组与对照组的免疫微环境差异,定位显著改变的免疫细胞类型与共浸润模式。

## ④ 特点 / 亮点

- **Turnkey**:一条命令出三图;自动分组 + 显著性标注。
- **顶刊图**:期刊配色箱线(带星号)+ 堆叠组成 + ComplexHeatmap 相关热图。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/Immune_boxplot.png` | 分组箱线 | 各细胞两组对比 + 显著性 |
| `assets/Immune_stackbar.png` | 堆叠柱状 | 样本免疫组成 |
| `assets/Immune_correlation.png` | 相关热图 | 免疫细胞共浸润 |

![boxplot](assets/Immune_boxplot.png)
![stack](assets/Immune_stackbar.png)

---

## 运行

```bash
Rscript 021_immune_visualization.R                              # 示例
Rscript 021_immune_visualization.R --input data/CIBERSORT_results.csv
```

## 依赖安装

```r
install.packages(c("ggpubr","ggplot2","reshape2","circlize"))
BiocManager::install("ComplexHeatmap")
```
