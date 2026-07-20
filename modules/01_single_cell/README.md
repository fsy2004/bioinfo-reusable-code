# 01 · 单细胞分析 — Single-cell analysis

本域共 32 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 上游与质控 — Pipeline & QC

- [023_rds_structure_check.R](01_pipeline_qc/023_rds_structure_check.R) — Inspect an RDS object's structure
- [024_scrna_rds_prep.R](01_pipeline_qc/024_scrna_rds_prep.R) — Build a Seurat/RDS object from raw data
- [025_scrna_data_prep.R](01_pipeline_qc/025_scrna_data_prep.R) — Read 10x data into a Seurat object
- [061_scfocal_gui_input_prep.R](01_pipeline_qc/061_scfocal_gui_input_prep.R) — Prepare scFOCAL GUI input + launch Shiny
- [562_mixhvg_hvg_selection](01_pipeline_qc/562_mixhvg_hvg_selection) — (用途待补)

## 整合与批次校正 — Integration & batch correction

- [506_scvi_scanvi_integration](02_integration_batch/506_scvi_scanvi_integration) — scVI/scANVI integration + label transfer (vs PCA baseline)
- [563_concord_contrastive_integration](02_integration_batch/563_concord_contrastive_integration) — (用途待补)
- [564_scextract_prior_integration](02_integration_batch/564_scextract_prior_integration) — (用途待补)
- [565_scmultibench_integration_benchmark](02_integration_batch/565_scmultibench_integration_benchmark) — (用途待补)

## 注释与细胞分型 — Annotation & cell typing

- [044_multimodalad_scrna.R](03_annotation_typing/044_multimodalad_scrna.R) — AD brain scRNA pipeline + Monocle pseudotime
- [046_scrna_publication_figures](03_annotation_typing/046_scrna_publication_figures) — Standard Seurat flow → publication figures
- [049_scrna_manual_annot_cellchat_trajectory.R](03_annotation_typing/049_scrna_manual_annot_cellchat_trajectory.R) — Manual annotation + CellChat + trajectory
- [566_phispace_soft_annotation](03_annotation_typing/566_phispace_soft_annotation) — (用途待补)

## 组成与丰度差异 — Composition / differential abundance

- [557_sccomp_composition_da](04_composition_da/557_sccomp_composition_da) — sccomp Bayesian beta-binomial cell-composition DA vs 3 baselines
- [558_milo_neighborhood_da](04_composition_da/558_milo_neighborhood_da) — Milo KNN-neighborhood differential abundance vs discrete-cluster baseline

## 差异表达(含 pseudobulk) — Differential expression

- [559_muscat_pseudobulk_ds](05_differential_expression/559_muscat_pseudobulk_ds) — muscat multi-sample pseudobulk differential-state vs cell-level baseline
- [567_glimes_mixed_effect_de](05_differential_expression/567_glimes_mixed_effect_de) — (用途待补)

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

- [569_nicheformer_sc_spatial_fm](10_foundation_models/569_nicheformer_sc_spatial_fm) — (用途待补)
- [570_epiagent_scatac_fm](10_foundation_models/570_epiagent_scatac_fm) — (用途待补)
- [571_captain_rna_protein_fm](10_foundation_models/571_captain_rna_protein_fm) — (用途待补)
- [572_cellvq_discrete_cell_fm](10_foundation_models/572_cellvq_discrete_cell_fm) — (用途待补)
