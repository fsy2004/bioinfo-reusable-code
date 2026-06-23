# 20 · 突变 / CNV / 甲基化 / 蛋白组 / 代谢组模板

本类别为**命令行 turnkey 模板**(已是 `--maf/--outdir` 等 argparse 风格,无硬编码路径),覆盖各组学标准分析。各模板依赖其专属的重型 Bioconductor 包,需在装好对应包的环境运行。

| 模板 | 用途 | 主依赖 | 典型输出图 |
|------|------|--------|-----------|
| `mutation_maftools_pipeline.R` | 体细胞突变汇总 | maftools | oncoplot · summary |
| `methylation_minfi_champ_pipeline.R` | 甲基化差异分析 | minfi · ChAMP · limma | M 值分布 · 热图 |
| `proteomics_limma_msstats_pipeline.R` | 蛋白组差异分析 | limma · MSstats | 火山 · 热图 |
| `metabolomics_metaboanalystR_pipeline.R` | 代谢组差异分析 | MetaboAnalystR | 火山 · 热图 |
| `cnv_gistic_or_cnvkit_pipeline.md` | CNV 分析流程说明 | GISTIC2 / CNVKit / inferCNV | — |

## 运行示例

```bash
Rscript mutation_maftools_pipeline.R --maf cohort.maf --outdir results/mutation
```

## 推荐链条
多组学矩阵 → 差异筛选 → [083 MOFA/DIABLO 整合](../19_多组学整合_分型模板/) → [084 NMF 分型](../19_多组学整合_分型模板/084_NMF共识聚类分型/) → 免疫/空间/预后解释。

> ⏭️ **依赖说明**:这些 Bioconductor 包(maftools/minfi/ChAMP/MSstats/MetaboAnalystR)体积大、依赖多,**本地未渲染示例图**;模板结构本身已 turnkey,装好对应包即可直接运行出图。
> 如需期刊级火山/热图,可复用 03 类([010](../03_GEO转录组整理与差异分析/))与 ComplexHeatmap 写法。
