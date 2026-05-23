# workflow_mr_coloc_twas

用途：MR、coloc、TWAS、pQTL/MVMR 组成的遗传因果和药物靶点证据链。

## 适用场景

- 候选暴露、免疫细胞、蛋白、基因表达或药物靶点与疾病之间需要因果证据。
- 需要把 bulk/scRNA 候选基因升级为遗传支持的治疗靶点。
- 需要连接 eQTL、pQTL、TWAS、coloc 和单细胞定位。

## 输入规范

| 输入 | 推荐格式 | 关键字段 |
|---|---|---|
| 暴露 GWAS | TSV/CSV/VCF | SNP、beta、se、effect_allele、other_allele、eaf、pval |
| 结局 GWAS | TSV/CSV/VCF | SNP、beta、se、alleles、pval |
| QTL | eQTL/pQTL summary | SNP、gene/protein、beta、se、pval |
| LD reference | PLINK 或内置 OpenGWAS | ancestry 一致 |
| TWAS 权重 | FUSION .RDat/weights | cell type、gene |

## 推荐顺序

| 步骤 | 脚本 | 核心输出 | 说明 |
|---:|---|---|---|
| 1 | [028_MR_GWAS暴露_VCF显著SNP筛选.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/028_MR_GWAS暴露_VCF显著SNP筛选.R>) | 显著 SNP | 暴露工具变量入口 |
| 2 | [029_MR_GWAS暴露_LD去除.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/029_MR_GWAS暴露_LD去除.R>) | clumped SNP | 去 LD |
| 3 | [030_MR_GWAS暴露_添加EAF.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/030_MR_GWAS暴露_添加EAF.R>) | 补全 EAF | harmonise 前检查 |
| 4 | [031_MR_GWAS暴露_弱工具变量过滤.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/031_MR_GWAS暴露_弱工具变量过滤.R>) | 强工具变量 | F 统计量 |
| 5 | [075_TwoSampleMR_coloc_药物靶点因果证据链.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/075_TwoSampleMR_coloc_药物靶点因果证据链.R>) | MR、异质性、多效性、coloc | 主证据链 |
| 6 | [055_免疫细胞_疾病MR方向性检验.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/055_免疫细胞_疾病MR方向性检验.R>) | 批量 MR 和 Steiger | 免疫细胞暴露 |
| 7 | [036-042 OneK1K/FUSION TWAS](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/10_TWAS_单细胞eQTL权重>) | TWAS 权重和结果 | 细胞类型表达遗传证据 |
| 8 | [079_pQTL_MVMR_蛋白中介MR.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/09_孟德尔随机化_GWAS处理/079_pQTL_MVMR_蛋白中介MR.R>) | MVMR、两步 MR、中介比例 | 蛋白和药物靶点强化 |
| 9 | 单细胞/空间定位 | cell type expression、spatial map | 回到机制位置 |

## 结果链条写法

候选暴露/基因/蛋白 -> 工具变量筛选 -> 两样本 MR -> 敏感性分析 -> coloc 共定位 -> TWAS 或 pQTL/MVMR -> 单细胞/空间定位 -> 药物重定位或分子对接。

## 质量控制检查

- 暴露和结局人群 ancestry 尽量一致。
- 工具变量 F 统计量应报告。
- 方向性用 Steiger 检查。
- MR 阳性后优先做 coloc，避免 LD 混杂。
- 多暴露高度相关时优先补 MVMR。
- 药物靶点 MR 要尽量加入 pQTL、蛋白定位和药物数据库证据。
