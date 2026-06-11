# Reusable Bioinformatics Code Library

Curated script library organized by analysis purpose. Covers disease-target networks, enrichment, transcriptomics, single-cell, spatial, Mendelian randomization, machine learning, molecular docking / MD, virtual perturbation, drug repurposing, multi-omics integration, and disease-burden analysis.

All scripts are reusable templates. Adapt input paths, phenotype labels, and package versions before running in a new project. Keep large data files outside this repository.

---

## Module Catalog

### 01 · 网络药理学与靶点数据库

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 001_CTD_化合物靶点提取.R | CTD 化合物-靶点关联提取 | R | — |
| 002_SwissTargetPrediction_化合物靶点提取.R | SwissTarget 靶点预测提取 | R | — |
| 003_CTD_SwissTargetPrediction_靶点并集_Venn.R | 两库靶点合并 Venn | R | Venn |
| 004_GeneCards_疾病靶点提取.R | GeneCards 疾病靶点提取 | R | — |
| 005_OMIM_GeneCards_疾病靶点并集_Venn.R | 两库疾病靶点 Venn | R | Venn |
| 006_疾病_化合物靶点交集_Venn.R | 疾病×化合物靶点交集 Venn | R | Venn |
| 011_差异基因_药物靶点交集_Venn_UpSet.R | DEG ∩ 药物靶点 Venn + UpSet | R | Venn / UpSet |
| 493_OpenTargets_DGIdb_ChEMBL_可成药性评分.py | 三库靶点可成药性综合评分 | Python | — |

### 02 · GO / KEGG 富集分析

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 007_GO_KEGG_富集分析.R | GO / KEGG 通路富集 | R | 柱状 / 气泡 / 点阵 |

### 03 · GEO 转录组整理与差异分析

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 008_GEO_表达矩阵整理.R | GEO 表达矩阵读入整理 | R | — |
| 009_GEO_样本分组整理.R | 样本分组注释整理 | R | — |
| 010_GEO_差异分析_火山热图PCA.R | limma DEG + 火山 / 热图 / PCA | R | 火山 / 热图 / PCA |
| 056_GEO_多队列合并_批次校正.R | SVA 批次校正多队列合并 | R | PCA |

### 04 · 机器学习特征基因筛选

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 012_LASSO_特征基因筛选.R | LASSO 回归特征筛选 | R | 系数图 / λ 折线 |
| 013_SVM_RFE_特征基因筛选.R | SVM-RFE 递归特征排序 | R | Forest |
| 014_RandomForest_特征基因筛选.R | 随机森林特征重要性 | R | 柱状 / Forest |
| 015_机器学习特征交集_Venn_UpSet.R | 多 ML 特征交集 Venn / UpSet | R | Venn / UpSet |
| 034_12种机器学习_特征基因筛选.R | 12-ML 组合筛选 + AUC 热图 | R | AUC 热图 / ROC |
| 035_5种机器学习组合_交集选择.R | 5-ML 组合交集特征 | R | Venn |
| 045_MultimodalAD_机器学习模型_EN.R | 弹性网 (Elastic Net) 模型 | R | ROC |
| 052_SHAP_机器学习解释分析.R | SHAP 11 种解释图 | Python | SHAP / waterfall / force / summary |
| 059_双疾病_15种机器学习175组合.R | 15-ML × 175-combo 双疾病 AUC | R | AUC 热图 / ROC |
| 496_Mime_101组合机器学习预后签名.R | Mime 框架 101 组合预后签名 | R | AUC 热图 / ROC / KM |

### 05 · 诊断模型与验证

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 016_诊断模型_ROC_校准_DCA.R | ROC + 校准曲线 + DCA + nomogram | R | ROC / DCA / 校准 / nomogram |
| 063_GEO_诊断模型验证.R | 独立队列完整诊断验证 | R | ROC / DCA / 校准 / nomogram / Forest |

### 06 · 免疫浸润与免疫可视化

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 017_免疫浸润_来源函数.R | CIBERSORT 来源函数 | R | — |
| 018_免疫浸润_评分分析.R | 免疫评分 + 相关矩阵 | R | 热图 / 相关矩阵 |
| 019_免疫浸润_通用来源函数.R | 通用去卷积来源函数 | R | — |
| 020_免疫浸润_通用评分分析.R | 通用评分分析 | R | 相关矩阵 |
| 021_免疫浸润_可视化.R | 免疫浸润柱状 / 箱图可视化 | R | 柱状 / 箱图 |
| 492_IOBR_多算法免疫去卷积.R | IOBR 多方法联合去卷积 | R | 热图 / 相关矩阵 |

### 07 · 分子对接与结合能可视化

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 022_分子对接_结合能气泡图.R | 对接结合能气泡图 | R | 气泡 / 散点 |
| 086_Vina_GROMACS_MMPBSA_MDAnalysis_对接动力学自动化.py | AutoDock Vina + GROMACS + gmx_MMPBSA 全流程 | Python | RMSD / RMSF / Rg / SASA / MM-PBSA 折线 |

### 08 · 单细胞 / 空间转录组 / 细胞轨迹

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 023_RDS对象_结构检查.R | scRNA RDS 对象结构检查 | R | — |
| 024_单细胞RDS_数据整理.R | 单细胞 QC 整理 | R | — |
| 025_单细胞_数据整理.R | 单细胞数据预处理 | R | — |
| 026_单细胞_Seurat全流程_scTenifoldKnk敲除.R | Seurat 全套 + scTenifoldKnk 虚拟敲除 | R | UMAP / violin / 热图 |
| 027_空间转录组_Seurat自动分析.R | Visium 空间转录组 Seurat 全流程 | R | 空间特征图 / UMAP |
| 044_MultimodalAD_单细胞分析_EN.R | AD scRNA-seq + monocle 拟时序 | R | UMAP / 轨迹 |
| 046_单细胞_单基因发表级图.R | 单基因出版级 UMAP / violin / alluvial / marker 热图 | R | UMAP / violin / alluvial / 热图 |
| 049_单细胞_手工注释_CellChat_轨迹.R | 手工注释 + CellChat + monocle 轨迹 | R | UMAP / CellChat 圆圈 / chord / 轨迹 |
| 050_空间转录组_聚类注释_轨迹.R | 空间转录组聚类注释 + 轨迹 | R | 空间图 / UMAP / 轨迹 |
| 051_单细胞_CellChat细胞通讯.R | CellChat 圆圈 / chord / LR-bubble 通讯 | R | 圆圈 / chord / 气泡 |
| 058_单细胞_Scissor疾病相关细胞.R | Scissor 疾病表型关联细胞识别 | R | UMAP |
| 061_scFOCAL_GUI输入准备_虚拟扰动入口.R | scFOCAL / CellOracle GUI 输入准备 | R | — |
| 062_scTour_拟时序向量场教程.py | scTour 速率场轨迹 | Python | 向量场 / UMAP |
| 082_轨迹多算法一致性_Slingshot_tradeSeq_CytoTRACE2.R | Slingshot + tradeSeq + CytoTRACE2 轨迹一致性 | R | 轨迹 / pseudotime |
| 082_Palantir_轨迹分支概率.py | Palantir 分支概率轨迹 | Python | 轨迹 / pseudotime |
| 491_SCTOUR_数据外文件/ | scTour 环境配置 + 官方教程运行脚本 | PS1/Py | UMAP / 向量场 |

### 09 · 孟德尔随机化 / GWAS 处理

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 028_MR_GWAS暴露_VCF显著SNP筛选.R | VCF 格式 GWAS 显著 SNP 筛选 | R | — |
| 029_MR_GWAS暴露_LD去除.R | LD clumping 工具变量精简 | R | — |
| 030_MR_GWAS暴露_添加EAF.R | 无 EAF GWAS 添加效应等位基因频率 | R | — |
| 031_MR_GWAS暴露_弱工具变量过滤.R | F 统计过滤弱工具变量 | R | — |
| 032_MR_GWAS暴露_FinnGen结局分析.R | TwoSampleMR 全套 MR 分析 | R | scatter / forest / funnel / leave-one-out |
| 033_MR_GWAS暴露_FinnGen结局分析_备用.R | MR 备用分析模板 | R | scatter / forest / funnel |
| 043_MultimodalAD_MendelianRandomization_EN.R | MR + 方向性检验 | R | scatter / forest |
| 055_免疫细胞_疾病MR方向性检验.R | 免疫细胞-疾病双向 MR 方向性检验 | R | Forest |
| 075_TwoSampleMR_coloc_药物靶点因果证据链.R | MR + colocalization 药物靶点因果链 | R | scatter / forest / funnel / LocusZoom |
| 079_pQTL_MVMR_蛋白中介MR.R | pQTL 驱动 MVMR 蛋白中介 MR | R | Forest |
| 497_lavaan_SEM_转录中介路径模型.R | lavaan SEM 路径图 | R | SEM 路径图 |

### 10 · TWAS / 单细胞 eQTL 权重

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 036_OneK1K_TWAS_同异质成分拟合.R | OneK1K sc-eQTL 同/异质 TWAS 权重拟合 | R | 表格 |
| 037_OneK1K_TWAS_成分模型拟合.R | TWAS 成分模型拟合 | R | 表格 |
| 038_OneK1K_TWAS_权重预处理.R | TWAS 权重预处理 | R | — |
| 039_OneK1K_TWAS_FUSION权重生成.R | FUSION TWAS 权重文件生成 | R | — |
| 040_FUSION_TWAS_targetC.R | FUSION TWAS targetC 模型 | R | 表格 |
| 041_FUSION_TWAS_S_targetC.R | FUSION TWAS S+targetC 模型 | R | 表格 |
| 042_FUSION_TWAS_S_allC.R | FUSION TWAS S+allC 模型 | R | 表格 |

### 11 · WGCNA 共表达网络

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 054_WGCNA_无符号共表达网络.R | 软阈值 / 树状图 / 模块-特征热图 / hub 散点 | R | 模块热图 / hub 散点 |

### 12 · TCGA 肿瘤预后生存 (仅参考)

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 048_TCGA_单基因生存曲线.R | 单基因 KM 生存曲线 | R | KM |
| 057_TCGA_预后风险模型可视化.R | 风险评分 riskplot + timeROC + KM | R | riskplot / timeROC / KM |
| 060_TCGA_免疫相关双蝴蝶图.R | 免疫相关双蝴蝶相关图 | R | butterfly |
| 497_scSurvival/ | 单细胞-队列联合生存分析包 (外部) | Python | 队列生存曲线 |

### 13 · 转录因子调控 / 基因组圈图

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 047_RcisTarget_转录因子Motif网络.R | RcisTarget motif-TF 网络 (igraph / visNetwork / Sankey) | R | 网络 / Sankey |
| 053_circlize_基因染色体圈图.R | circlize 基因组染色体圈图 | R | circos |
| 081_pySCENIC_Regulon_TF活性.py | pySCENIC regulon + TF 活性 UMAP | Python | UMAP / 热图 |

### 14 · 单细胞虚拟扰动 / 扰动数据库

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 067_scPerturb_扰动数据Etest.R | scPerturb 数据库 E-test 统计检验 | R | 统计表 |
| 068_GEARS_单细胞组合扰动预测.py | GEARS 深度学习组合扰动预测 | Python | 散点 / 热图 |
| 069_CellOracle_GRN虚拟扰动.py | CellOracle GRN + in-silico TF KO + 流场 | Python | UMAP / 流场 / delta 排名 |
| 085_Squidiff_扩散模型单细胞扰动预测.py | Squidiff 扩散模型扰动预测 | Python | UMAP |
| 494_GenKI_VGAE虚拟敲除.py | GenKI VGAE 图神经网络虚拟敲除 | Python | 热图 / 散点 |
| 495_bulkVGK_scTenifoldKnk网络DI.R | 批量 scTenifoldKnk 扰动 DI 指数 | R | 排名柱状 |

### 15 · 药物扰动 / 药物重定位

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 070_chemCPA_药物扰动表达预测.py | chemCPA 药物扰动表达谱预测 | Python | 散点 / 热图 |
| 071_scDrug_单细胞药物响应预测.py | scDrug 单细胞药物响应预测 | Python | UMAP / 热图 |
| 078_FAERS_ROR_PRR_BCPNN_EBGM药物警戒.R | FAERS 四种信号挖掘算法药物警戒 | R | Forest / ROR 图 |

### 16 · 空间通讯 / 细胞命运

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 072_CellRank_命运概率与驱动基因.py | CellRank 细胞命运概率 + 驱动基因 | Python | UMAP / 轨迹 |
| 073_COMMOT_空间细胞通讯.py | COMMOT 空间分解细胞通讯 | Python | 空间通讯图 |
| 074_Tangram_单细胞到空间映射.py | Tangram sc → 空间映射 | Python | 空间特征图 |
| 076_decoupler_TF通路活性评分.py | decoupleR TF / 通路活性评分 | Python | violin / UMAP / 热图 |
| 077_NicheNet_配体靶基因通信推断.R | NicheNet 配体-靶基因通信推断 | R | 热图 / 网络 |
| 080_cell2location_Squidpy_空间生态位.py | cell2location + Squidpy 空间生态位 | Python | 空间图 |

### 17 · 高级结果图与闭环可视化

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 498_ggalluvial_桑基冲积图.R | ggalluvial 桑基冲积图 | R | Sankey / alluvial |

### 18 · 外部方法源码 (待整合)

`14_AI科学示意图生成/AutoFigure-Edit/` — AI 驱动 SVG 科学示意图编辑工具 (外部项目)。`autofigure2.py` 为核心脚本。

### 19 · 多组学整合 / 分型模板

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| 083_MOFA_DIABLO_多组学潜变量整合.R | MOFA + DIABLO 多组学潜变量整合 | R | 因子图 / 热图 |
| 084_NMF_ConsensusClusterPlus_共浸润分型.R | NMF + ConsensusClusterPlus 免疫分型 | R | 热图 / KM |

### 20 · 突变 / CNV / 甲基化 / 蛋白组 / 代谢组模板

| Script | Purpose | Lang | Output figures |
|--------|---------|------|----------------|
| mutation_maftools_pipeline.R | maftools 体细胞突变分析 | R | lollipop / oncoprint |
| methylation_minfi_champ_pipeline.R | minfi + ChAMP 甲基化分析 | R | M 值分布 / 热图 |
| metabolomics_metaboanalystR_pipeline.R | MetaboAnalystR 代谢组分析 | R | 火山 / 热图 |
| proteomics_limma_msstats_pipeline.R | limma + MSstats 蛋白组差异分析 | R | 火山 / 热图 |
| cnv_gistic_or_cnvkit_pipeline.md | GISTIC2 / CNVKit CNV 分析流程说明 | — | — |

### 21 · 疾病负担 / 共病 / GBD-NHANES-CHARLS

`01_GBD/` `02_NHANES/` `03_CHARLS/` `04_共病网络与模式/` 各子目录含 README 和数据处理规范。  
`99_外部源码_待整合/` 含 GBD2021、R-script-for-GBD、comorbidity_networks、CTBN-Multimorbidity-Paper、nhanes 包等外部参考实现。

---

## Figure-type → Module Quick Lookup

| Figure type | Modules |
|-------------|---------|
| UMAP / tSNE | 026 027 044 046 050 051 058 062 072 074 076 081 082 |
| Violin / split violin | 026 046 049 050 076 |
| Heatmap / ComplexHeatmap | 026 051 054 059 060 063 083 084 492 496 |
| ROC | 016 034 045 052 057 059 063 |
| DCA / calibration / nomogram | 016 063 |
| Forest / lollipop | 013 014 016 032 033 034 052 055 063 075 079 078 |
| KM / survival | 048 057 063 083 084 496 |
| Volcano | 010 083 084 |
| Bubble / dot | 007 010 022 026 051 052 |
| Venn / UpSet | 003 005 006 011 015 034 035 |
| PCA | 010 026 027 056 |
| Scatter (MR) | 032 033 043 054 |
| Circos / chord | 051 053 |
| Network (igraph / visNetwork) | 047 051 077 081 |
| Correlation matrix | 018 020 021 060 492 |
| MR funnel / Manhattan / QQ / radial | 032 033 043 |
| Docking / binding energy | 022 086 |
| MD (RMSD / RMSF / MM-PBSA) | 086 |
| Spatial feature / niche maps | 027 050 073 074 080 |
| Alluvial / Sankey | 046 047 049 498 |
| SHAP / waterfall / force | 052 |
| Pseudotime / trajectory | 044 050 062 072 082 |
| SEM path diagram | 497 |
| Pharmacovigilance forest | 078 |
| Virtual-KO delta ranking | 069 494 495 |
| AUC heatmap (ML leaderboard) | 034 059 496 |
