# 01 · 单细胞分析 — Single-cell analysis

本域共 33 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 上游与质控 — Pipeline & QC

- [023_rds_structure_check.R](01_pipeline_qc/023_rds_structure_check.R) — Inspect an RDS object's structure
- [024_scrna_rds_prep.R](01_pipeline_qc/024_scrna_rds_prep.R) — Build a Seurat/RDS object from raw data
- [025_scrna_data_prep.R](01_pipeline_qc/025_scrna_data_prep.R) — Read 10x data into a Seurat object
- [061_scfocal_gui_input_prep.R](01_pipeline_qc/061_scfocal_gui_input_prep.R) — Prepare scFOCAL GUI input + launch Shiny
- [562_mixhvg_hvg_selection](01_pipeline_qc/562_mixhvg_hvg_selection) — 混合多种 HVG 打分方法(按秩取max)选高变基因,自带 Seurat vst 基线 + ground-truth recall / ARI / silhouette 评测

## 整合与批次校正 — Integration & batch correction

- [506_scvi_scanvi_integration](02_integration_batch/506_scvi_scanvi_integration) — scVI/scANVI integration + label transfer (vs PCA baseline)
- [563_concord_contrastive_integration](02_integration_batch/563_concord_contrastive_integration) — 多批次单细胞整合模块:3 个本机可跑基线(PCA/批次中心化/ComBat+PCA)+ 守卫式 CONCORD 接口,用「批次混合熵 / 生物保真(kNN纯度+ARI/NMI) / 全局几何(trustworthiness+成对距离 Spearman)」三类指标同台对照评估。
- [564_scextract_prior_integration](02_integration_batch/564_scextract_prior_integration) — 多批次 scRNA 整合评测:以「批次混合熵 × 细胞类型 kNN 保真度 × 稀有类型保真度」双轴对比未校正 PCA / ComBat(/Harmony),并守卫式封装 scExtract 的 scanorama_prior / cellhint_prior 先验整合
- [565_scmultibench_integration_benchmark](02_integration_batch/565_scmultibench_integration_benchmark) — 把 scMultiBench 的 scIB 评测层封装成模块:给任意整合 embedding 打生物保留/批次校正/综合分,强制与朴素 PCA 基线对比,出热图+权衡散点+lollipop 排名。

## 注释与细胞分型 — Annotation & cell typing

- [044_multimodalad_scrna.R](03_annotation_typing/044_multimodalad_scrna.R) — AD brain scRNA pipeline + Monocle pseudotime
- [046_scrna_publication_figures](03_annotation_typing/046_scrna_publication_figures) — Standard Seurat flow → publication figures
- [049_scrna_manual_annot_cellchat_trajectory.R](03_annotation_typing/049_scrna_manual_annot_cellchat_trajectory.R) — Manual annotation + CellChat + trajectory
- [566_phispace_soft_annotation](03_annotation_typing/566_phispace_soft_annotation) — PhiSpace 连续表型软注释：把 query 细胞投影到参考类型张成的表型空间，给出每细胞×每类型的连续得分而非硬标签；自带质心相关 + PCA 回归两条本机可跑基线与已知混合比例真值评估，PhiSpace 本体为守卫式封装。

## 组成与丰度差异 — Composition / differential abundance

- [557_sccomp_composition_da](04_composition_da/557_sccomp_composition_da) — sccomp Bayesian beta-binomial cell-composition DA vs 3 baselines
- [558_milo_neighborhood_da](04_composition_da/558_milo_neighborhood_da) — Milo KNN-neighborhood differential abundance vs discrete-cluster baseline

## 差异表达(含 pseudobulk) — Differential expression

- [559_muscat_pseudobulk_ds](05_differential_expression/559_muscat_pseudobulk_ds) — muscat multi-sample pseudobulk differential-state vs cell-level baseline
- [567_glimes_mixed_effect_de](05_differential_expression/567_glimes_mixed_effect_de) — 多供体单细胞原始 UMI 计数上的 Poisson-GLMM(供体随机截距)差异表达,与朴素细胞级 t 检验、pseudobulk 两条基线同台对比,量化供体伪重复造成的一类错误膨胀。

## 轨迹与 RNA 速率 — Trajectory & RNA velocity

- [062_sctour_pseudotime_vectorfield.py](06_trajectory_velocity/062_sctour_pseudotime_vectorfield.py) — scTour pseudotime + latent space + vector field
- [072_cellrank_fate_drivers.py](06_trajectory_velocity/072_cellrank_fate_drivers.py) — scVelo + CellRank fate probabilities & drivers
- [082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R](06_trajectory_velocity/082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R) — Slingshot / tradeSeq / CytoTRACE2 trajectory consensus
- [087_palantir_branch_probability.py](06_trajectory_velocity/087_palantir_branch_probability.py) — Palantir pseudotime + branch probability + entropy
- [491_sctour_extra_files](06_trajectory_velocity/491_sctour_extra_files) — scTour 官方教程的复现脚本与环境记录(062 的配套材料)
- [517_vector_trajectory_direction](06_trajectory_velocity/517_vector_trajectory_direction) — VECTOR expression-potential differentiation direction

## 拷贝数与克隆 — CNV & clonality

- [560_copykat_scrna_cnv](07_cnv_clonality/560_copykat_scrna_cnv) — copyKAT scRNA CNV inference + aneuploid/diploid calling

## 通路/转录因子活性打分 — Pathway & TF activity

- [076_decoupler_tf_pathway_activity.py](08_activity_scoring/076_decoupler_tf_pathway_activity.py) — decoupler TF / pathway activity inference
- [510_scmetabolism_pathway_activity](08_activity_scoring/510_scmetabolism_pathway_activity) — Single-cell metabolic pathway activity (AUCell/UCell-style)

## 单细胞↔bulk 表型关联 — Single-cell to bulk phenotype

- [058_scrna_scissor.R](09_bulk_phenotype_link/058_scrna_scissor.R) — Scissor — link bulk phenotype to disease-relevant cells
- [497_scsurvival_cohort](09_bulk_phenotype_link/497_scsurvival_cohort) — scSurvival — single-cell cohort survival (vendored pkg)

## 单细胞基础模型 — Foundation models

- [568_scprint_foundation_grn](10_foundation_models/568_scprint_foundation_grn) — scPRINT 单细胞基础模型三任务(GRN 推断 / 去噪 / 细胞嵌入与标签预测)的本机可跑朴素基线 + 守卫式官方 API 封装,评估口径逐条对齐上游源码
- [569_nicheformer_sc_spatial_fm](10_foundation_models/569_nicheformer_sc_spatial_fm) — Nicheformer(单细胞+空间联合预训练基础模型)守卫式封装 + 本机可跑的线性对照基线:比较 intrinsic 表达 PCA 与 niche-aware(⊕空间 kNN 邻域均值 PCA)在 niche / cell-type 标签上的同折 CV macro-F1,并做解离参考→空间 query 的跨模态标签迁移地板值。
- [570_epiagent_scatac_fm](10_foundation_models/570_epiagent_scatac_fm) — scATAC 细胞×cCRE 矩阵 → TF-IDF+SVD(LSI) 基线做嵌入/聚类/细胞类型预测/填补/批次混合评估,并守卫式封装 EpiAgent 基础模型路径(仅环境探测,不臆造调用)
- [571_captain_rna_protein_fm](10_foundation_models/571_captain_rna_protein_fm) — 配对 CITE-seq「RNA→表面蛋白」填补基准台：matched-gene 与 PCA+Ridge 两条防泄漏基线 + CAPTAIN 本体守卫式探测（含 Drive 占位符识别）
- [572_cellvq_discrete_cell_fm](10_foundation_models/572_cellvq_discrete_cell_fm) — 离散「细胞词表」(VQ codebook)单细胞表征模块：PCA+k-means 码本基线量化离散化的信息损失/码本坍缩/重构 R²，并对官方 CellVQ 提供逐符号核实过的守卫式封装。
