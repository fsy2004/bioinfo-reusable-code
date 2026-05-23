# 20_突变_CNV_甲基化_蛋白组_代谢组模板

本模块用于补齐当前代码库中突变、CNV、甲基化、蛋白组和代谢组的标准入口。

## 脚本

| 脚本 | 作用 |
|---|---|
| `mutation_maftools_pipeline.R` | 读取 MAF，输出突变样本/基因汇总和 oncoplot。 |
| `cnv_gistic_or_cnvkit_pipeline.md` | 记录 GISTIC2、CNVkit、inferCNV 的标准输入输出和封装建议。 |
| `methylation_minfi_champ_pipeline.R` | 对 beta-value 矩阵做 limma 差异甲基化位点分析。 |
| `proteomics_limma_msstats_pipeline.R` | 对蛋白矩阵做 limma 差异蛋白分析。 |
| `metabolomics_metaboanalystR_pipeline.R` | 对代谢物矩阵做两组差异代谢物分析。 |

## 推荐链条

多组学矩阵 -> 差异层面筛选 -> MOFA/DIABLO 潜变量整合 -> NMF/分型 -> 免疫、空间、药敏、MR 或预后解释。
