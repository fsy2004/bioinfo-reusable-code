# Module catalog

按 **域 → 子类** 两级组织。想做某类分析时,先定位域,再在子类里挑模块;
每行给出用途、输入→输出、依赖、语言与产出图型。新项目脚手架与统一出图样式见
[`_framework/`](_framework/)。

> **复用优先,绝不从头写。** 先在本目录挑模块(或挑一个真实已发表工具)再适配,
> 不要凭记忆手写分析代码 —— 那会带来假 API 和错参数。见
> [`_framework/CONVENTIONS.md` §0](_framework/CONVENTIONS.md)。

## 状态图例

| 标记 | 含义 |
|------|------|
| ✅ | 开箱即跑 —— 用自带合成示例数据本机跑通,零改动 |
| 🟡 | 核心复现或诚实基线本机可跑;完整方法需在分析服务器装包(见 [`_framework/SERVER_DEPENDENCIES.md`](_framework/SERVER_DEPENDENCIES.md)) |
| 🔴 | 重型 / GPU / 外部工具链 —— 守卫式引用封装,不在本机渲染 |
| 📄 | 模板或上游脚本 —— 自带数据 + 自行安装,无捆绑示例 |
| 📦 | Vendored 第三方包 —— 仅保留清单 / 本地 |
| 🗃️ | 仅本地,git 忽略 —— 不在公开仓库中 |

---


## 01 · 单细胞分析 — Single-cell analysis  (32)

### 01.01 · 上游与质控 — Pipeline & QC

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 023 | [023_rds_structure_check.R](01_single_cell/01_pipeline_qc/023_rds_structure_check.R) | Inspect an RDS object's structure | RDS → structure print + csv | R base | R | — |
| 📄 | 024 | [024_scrna_rds_prep.R](01_single_cell/01_pipeline_qc/024_scrna_rds_prep.R) | Build a Seurat/RDS object from raw data | raw → Seurat object | R · Seurat, SingleR, celldex | R | — |
| 📄 | 025 | [025_scrna_data_prep.R](01_single_cell/01_pipeline_qc/025_scrna_data_prep.R) | Read 10x data into a Seurat object | 10x raw → Seurat object | R · Seurat, Matrix | R | — |
| 📄 | 061 | [061_scfocal_gui_input_prep.R](01_single_cell/01_pipeline_qc/061_scfocal_gui_input_prep.R) | Prepare scFOCAL GUI input + launch Shiny | RData + map csv → RDS + GUI | R · Seurat, scFOCAL, shiny | R | — (interactive) |
| — | 562 | [562_mixhvg_hvg_selection](01_single_cell/01_pipeline_qc/562_mixhvg_hvg_selection) | — | — | — | — | — |

### 01.02 · 整合与批次校正 — Integration & batch correction

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 506 | [506_scvi_scanvi_integration](01_single_cell/02_integration_batch/506_scvi_scanvi_integration) | scVI/scANVI integration + label transfer (vs PCA baseline) | h5ad (batch/label) → integration + labels | Py · scvi-tools, scanpy, sklearn | Python | UMAP, scatter, heatmap (confusion) |
| — | 563 | [563_concord_contrastive_integration](01_single_cell/02_integration_batch/563_concord_contrastive_integration) | — | — | — | — | — |
| — | 564 | [564_scextract_prior_integration](01_single_cell/02_integration_batch/564_scextract_prior_integration) | — | — | — | — | — |
| — | 565 | [565_scmultibench_integration_benchmark](01_single_cell/02_integration_batch/565_scmultibench_integration_benchmark) | — | — | — | — | — |

### 01.03 · 注释与细胞分型 — Annotation & cell typing

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 044 | [044_multimodalad_scrna.R](01_single_cell/03_annotation_typing/044_multimodalad_scrna.R) | AD brain scRNA pipeline + Monocle pseudotime | GSE157827 → object + trajectory | R · Seurat, SingleR, monocle | R | violin, UMAP, trajectory, feature-map |
| ✅ | 046 | [046_scrna_publication_figures](01_single_cell/03_annotation_typing/046_scrna_publication_figures) | Standard Seurat flow → publication figures | counts.csv → object + figures | R · Seurat, ggplot2 | R | UMAP, dotplot, heatmap, feature-map, violin |
| 📄 | 049 | [049_scrna_manual_annot_cellchat_trajectory.R](01_single_cell/03_annotation_typing/049_scrna_manual_annot_cellchat_trajectory.R) | Manual annotation + CellChat + trajectory | Seurat + markers → annotation + figures | R · Seurat, CellChat, monocle3 | R | violin, UMAP, heatmap, dotplot, trajectory |
| — | 566 | [566_phispace_soft_annotation](01_single_cell/03_annotation_typing/566_phispace_soft_annotation) | — | — | — | — | — |

### 01.04 · 组成与丰度差异 — Composition / differential abundance

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 557 | [557_sccomp_composition_da](01_single_cell/04_composition_da/557_sccomp_composition_da) | sccomp Bayesian beta-binomial cell-composition DA vs 3 baselines | composition counts → DA table + plots | R · sccomp, voomCLR, limma, ggbeeswarm | R | boxplot, lollipop, raincloud, dot-matrix |
| 🟡 | 558 | [558_milo_neighborhood_da](01_single_cell/04_composition_da/558_milo_neighborhood_da) | Milo KNN-neighborhood differential abundance vs discrete-cluster baseline | SCE (reducedDim + condition) → DA table | R · miloR (or BiocNeighbors baseline), igraph, ggbeeswarm | R | beeswarm, network, volcano, violin |

### 01.05 · 差异表达(含 pseudobulk) — Differential expression

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 559 | [559_muscat_pseudobulk_ds](01_single_cell/05_differential_expression/559_muscat_pseudobulk_ds) | muscat multi-sample pseudobulk differential-state vs cell-level baseline | SCE (.rds) → DS table + plots | R · muscat, SingleCellExperiment, edgeR, limma | R | MDS, volcano, heatmap, lollipop, dumbbell, raincloud |
| — | 567 | [567_glimes_mixed_effect_de](01_single_cell/05_differential_expression/567_glimes_mixed_effect_de) | — | — | — | — | — |

### 01.06 · 轨迹与 RNA 速率 — Trajectory & RNA velocity

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 062 | [062_sctour_pseudotime_vectorfield.py](01_single_cell/06_trajectory_velocity/062_sctour_pseudotime_vectorfield.py) | scTour pseudotime + latent space + vector field | h5ad → pseudotime + vector field | Py · scanpy, sctour | Python | UMAP, pseudotime, vector-field |
| 🔴 | 072 | [072_cellrank_fate_drivers.py](01_single_cell/06_trajectory_velocity/072_cellrank_fate_drivers.py) | scVelo + CellRank fate probabilities & drivers | h5ad (spliced) → driver genes | Py · scanpy, scvelo, cellrank | Python | — |
| 📄 | 082 | [082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R](01_single_cell/06_trajectory_velocity/082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R) | Slingshot / tradeSeq / CytoTRACE2 trajectory consensus | Seurat RDS → pseudotime table | R · slingshot, tradeSeq, CytoTRACE2 | R | pseudotime |
| 📄 | 087 | [087_palantir_branch_probability.py](01_single_cell/06_trajectory_velocity/087_palantir_branch_probability.py) | Palantir pseudotime + branch probability + entropy | h5ad + root → pseudotime/branch csv | Py · palantir, scanpy | Python | — (tables) |
| 📄 | 491 | [491_sctour_extra_files](01_single_cell/06_trajectory_velocity/491_sctour_extra_files) | scTour 官方教程的复现脚本与环境记录(062 的配套材料) | 教程数据 → 复现结果 + 环境说明 | Py · sctour | Python/PS | — |
| ✅ | 517 | [517_vector_trajectory_direction](01_single_cell/06_trajectory_velocity/517_vector_trajectory_direction) | VECTOR expression-potential differentiation direction | embedding + expr → potential + field | R · ggplot2 | R | vector-field |

### 01.07 · 拷贝数与克隆 — CNV & clonality

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 560 | [560_copykat_scrna_cnv](01_single_cell/07_cnv_clonality/560_copykat_scrna_cnv) | copyKAT scRNA CNV inference + aneuploid/diploid calling | gene×cell counts → prediction + CNAmat | R · copykat, ggplot2, uwot | R | heatmap (CNV), scatter (embedding), lollipop, heatmap (confusion) |

### 01.08 · 通路/转录因子活性打分 — Pathway & TF activity

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 076 | [076_decoupler_tf_pathway_activity.py](01_single_cell/08_activity_scoring/076_decoupler_tf_pathway_activity.py) | decoupler TF / pathway activity inference | h5ad → activity scores | Py · decoupler, scanpy | Python | — (downstream heatmap) |
| ✅ | 510 | [510_scmetabolism_pathway_activity](01_single_cell/08_activity_scoring/510_scmetabolism_pathway_activity) | Single-cell metabolic pathway activity (AUCell/UCell-style) | expr + meta (+gmt) → activity + figures | R · ggplot2 (self-contained scorer) | R | dotplot, row-z heatmap, distribution |

### 01.09 · 单细胞↔bulk 表型关联 — Single-cell to bulk phenotype

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 058 | [058_scrna_scissor.R](01_single_cell/09_bulk_phenotype_link/058_scrna_scissor.R) | Scissor — link bulk phenotype to disease-relevant cells | bulk pheno + Seurat → Scissor cells | R · Seurat, Scissor | R | UMAP, lollipop |
| 📦 | 497 | [497_scsurvival_cohort](01_single_cell/09_bulk_phenotype_link/497_scsurvival_cohort) | scSurvival — single-cell cohort survival (vendored pkg) | sc cohort + survival → risk model | Py · PyTorch, scanpy, lifelines | Python | cohort-survival |

### 01.10 · 单细胞基础模型 — Foundation models

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| — | 569 | [569_nicheformer_sc_spatial_fm](01_single_cell/10_foundation_models/569_nicheformer_sc_spatial_fm) | — | — | — | — | — |
| — | 570 | [570_epiagent_scatac_fm](01_single_cell/10_foundation_models/570_epiagent_scatac_fm) | — | — | — | — | — |
| — | 571 | [571_captain_rna_protein_fm](01_single_cell/10_foundation_models/571_captain_rna_protein_fm) | — | — | — | — | — |
| — | 572 | [572_cellvq_discrete_cell_fm](01_single_cell/10_foundation_models/572_cellvq_discrete_cell_fm) | — | — | — | — | — |

## 02 · 空间转录组 — Spatial transcriptomics  (21)

### 02.01 · 上游与细胞分割 — Pipeline & segmentation

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 027 | [027_spatial_seurat_auto.R](02_spatial_transcriptomics/01_pipeline_segmentation/027_spatial_seurat_auto.R) | Spatial transcriptomics read/cluster/visualize | Visium h5 → spatial figures | R · Seurat, SingleR, glmGamPoi | R | PCA, violin, UMAP, feature-map, niche-map |
| 📄 | 050 | [050_spatial_cluster_annot_trajectory.R](02_spatial_transcriptomics/01_pipeline_segmentation/050_spatial_cluster_annot_trajectory.R) | Spatial cluster annotation + monocle3 pseudotime | Visium → spatial annotation + trajectory | R · Seurat, monocle3, patchwork | R | niche-map, violin, UMAP, feature-map, pseudotime |
| — | 573 | [573_proseg_cell_segmentation](02_spatial_transcriptomics/01_pipeline_segmentation/573_proseg_cell_segmentation) | — | — | — | — | — |

### 02.02 · 空间域、空间可变基因与空间统计 — Domains, SVG & spatial statistics

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 505 | [505_spatial_advanced](02_spatial_transcriptomics/02_domains_svg_stats/505_spatial_advanced) | Spatial advanced: RCTD deconv + NMF niche + interface degree | sc ref + spatial rds → fractions/niche | R · spacexr, RcppML, mistyR | R | niche-map (×3) |
| ✅ | 541 | [541_banksy_spatial_domains](02_spatial_transcriptomics/02_domains_svg_stats/541_banksy_spatial_domains) | BANKSY neighbor-augmented spatial-domain segmentation vs non-spatial baseline | spatial csv → domains + ARI | R · Banksy, SpatialExperiment, aricode | R | spatial-scatter, lollipop, UMAP |
| ✅ | 542 | [542_nnsvg_spatial_svg](02_spatial_transcriptomics/02_domains_svg_stats/542_nnsvg_spatial_svg) | nnSVG spatially-variable genes (NNGP) vs non-spatial HVG baseline | counts + coords → SVG ranking | R · nnSVG, SpatialExperiment, scran | R | spatial-scatter, lollipop, scatter, violin |
| ✅ | 543 | [543_squidpy_spatial_statistics](02_spatial_transcriptomics/02_domains_svg_stats/543_squidpy_spatial_statistics) | squidpy spatial stats (Moran / nhood enrichment / co-occurrence / Ripley) | h5ad (spatial) → stats tables + plots | Py · squidpy, anndata, scanpy | Python | heatmap, lollipop, scatter, spatial-scatter |
| — | 574 | [574_stair_spatial_integration](02_spatial_transcriptomics/02_domains_svg_stats/574_stair_spatial_integration) | — | — | — | — | — |
| — | 575 | [575_scale_spatial_method](02_spatial_transcriptomics/02_domains_svg_stats/575_scale_spatial_method) | — | — | — | — | — |

### 02.03 · 解卷积与单细胞映射 — Deconvolution & mapping

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 074 | [074_tangram_sc_to_spatial.py](02_spatial_transcriptomics/03_deconvolution_mapping/074_tangram_sc_to_spatial.py) | Tangram single-cell → spatial mapping | sc + spatial h5ad → mapped h5ad | Py · tangram, scanpy | Python | — |
| 🔴 | 080 | [080_cell2location_squidpy_niche.py](02_spatial_transcriptomics/03_deconvolution_mapping/080_cell2location_squidpy_niche.py) | cell2location abundance + Squidpy neighborhood (GPU) | spatial h5ad + abundance → z-scores | Py · cell2location, squidpy | Python | — (downstream niche-map) |
| ✅ | 545 | [545_spotlight_deconvolution](02_spatial_transcriptomics/03_deconvolution_mapping/545_spotlight_deconvolution) | SPOTlight (NMF+NNLS) spot cell-type deconvolution vs known-mix baseline | scRNA ref + spatial → spot proportions | R · SPOTlight, SingleCellExperiment, scatterpie | R | scatterpie, heatmap, scatter, violin |

### 02.04 · 切片配准与三维重建 — Slice alignment & 3D

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 544 | [544_paste2_slice_alignment](02_spatial_transcriptomics/04_alignment_3d/544_paste2_slice_alignment) | PASTE optimal-transport spatial slice alignment / 3D stacking | two spatial h5ad → coupling + transform | Py · paste-bio, POT, anndata | Python | spatial-scatter, heatmap, violin, dumbbell |

### 02.05 · 细胞通讯 — Cell-cell communication

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 051 | [051_scrna_cellchat.R](02_spatial_transcriptomics/05_cell_communication/051_scrna_cellchat.R) | CellChat cell-communication network | annotated expr + labels → comm. figures | R · Seurat, CellChat | R | circos, bubble |
| 🔴 | 073 | [073_commot_spatial_communication.py](02_spatial_transcriptomics/05_cell_communication/073_commot_spatial_communication.py) | COMMOT spatial ligand–receptor communication | spatial h5ad + LR DB → comm. scores | Py · commot, scanpy | Python | — (downstream niche-map) |
| 🔴 | 077 | [077_nichenet_ligand_target.R](02_spatial_transcriptomics/05_cell_communication/077_nichenet_ligand_target.R) | NicheNet ligand activity + ligand–target links | receiver DE + prior → ligand activity | R · nichenetr | R | — (downstream heatmap) |
| ✅ | 509 | [509_communication_functional_loop](02_spatial_transcriptomics/05_cell_communication/509_communication_functional_loop) | Communication functional loop: ligand→UCell→enrich→Venn | receptor expr + prior + group → consensus | R · UCell, ggplot2 | R | lollipop, violin, venn |
| ✅ | 531 | [531_liana_consensus_cci](02_spatial_transcriptomics/05_cell_communication/531_liana_consensus_cci) | LIANA+ rank-aggregate consensus cell-cell communication (6 methods) | scRNA h5ad (celltype) → consensus L-R ranks | Py · liana, scanpy, anndata, plotnine | Python | dotplot, network, heatmap, lollipop, tile |
| — | 576 | [576_cellnest_spatial_ccc](02_spatial_transcriptomics/05_cell_communication/576_cellnest_spatial_ccc) | — | — | — | — | — |
| — | 577 | [577_spider_spatial_ccc](02_spatial_transcriptomics/05_cell_communication/577_spider_spatial_ccc) | — | — | — | — | — |

### 02.06 · 空间多组学 — Spatial multi-omics

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 521 | [521_spatialglue_multiomics](02_spatial_transcriptomics/06_spatial_multiomics/521_spatialglue_multiomics) | SpatialGlue spatial multi-omics domains (GNN; baseline local) | RNA + ADT grid → ARI + domains | Py · sklearn (baseline); SpatialGlue, torch-geometric | Python | spatial-scatter, lollipop |

## 03 · 虚拟扰动技术 — Virtual perturbation  (15)

### 03.01 · 虚拟敲除与扰动模拟 — In-silico knockout

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 026 | [026_scrna_seurat_sctenifoldknk_ko.R](03_virtual_perturbation/01_insilico_knockout/026_scrna_seurat_sctenifoldknk_ko.R) | Full QC/cluster/annotate + virtual KO pipeline | counts → object + full figure set | R · Seurat, monocle3, scTenifoldKnk, CellChat | R | UMAP, heatmap, PCA, violin, dot, feature-map |
| 🟡 | 067 | [067_scperturb_etest.R](03_virtual_perturbation/01_insilico_knockout/067_scperturb_etest.R) | scPerturb perturbation distance + E-test (light) | Seurat rds → edistance/etest | R · scperturbR | R | — |
| 🔴 | 068 | [068_gears_combo_perturbation.py](03_virtual_perturbation/01_insilico_knockout/068_gears_combo_perturbation.py) | GEARS single/combo perturbation prediction (GPU) | h5ad + perturb list → predictions | Py · GEARS, torch | Python | — |
| 🔴 | 069 | [069_celloracle_grn_perturbation.py](03_virtual_perturbation/01_insilico_knockout/069_celloracle_grn_perturbation.py) | CellOracle GRN virtual knockout (heavy) | Oracle pkl + genes → perturbed state | Py · celloracle | Python | — |
| 🔴 | 085 | [085_squidiff_diffusion_perturbation.py](03_virtual_perturbation/01_insilico_knockout/085_squidiff_diffusion_perturbation.py) | Squidiff/PerturbDiff diffusion perturbation (GPU) | h5ad + config → predictions | Py · Squidiff, torch | Python | — |
| 🔴 | 494 | [494_genki_vgae_ko.py](03_virtual_perturbation/01_insilico_knockout/494_genki_vgae_ko.py) | GenKI graph-VGAE virtual KO (KL ranking) | adata.h5ad + targets → KL ranking | Py · GenKI, torch-geometric | Python | — |
| ✅ | 495 | [495_bulkvgk_sctenifoldknk_di.R](03_virtual_perturbation/01_insilico_knockout/495_bulkvgk_sctenifoldknk_di.R) | Bulk co-expression virtual KO + differential influence | two-group expr → per-gene DI ranking | R · igraph | R | scatter (DE vs DI) |
| 🟡 | 507 | [507_geneformer_insilico](03_virtual_perturbation/01_insilico_knockout/507_geneformer_insilico) | Geneformer zero-shot embedding + in-silico deletion (baseline local) | counts/tokenized → KO ranking | Py · scanpy, sklearn (baseline); geneformer, torch | Python | UMAP, lollipop |
| — | 561 | [561_regvelo_grn_velocity](03_virtual_perturbation/01_insilico_knockout/561_regvelo_grn_velocity) | — | — | — | — | — |

### 03.02 · 基因调控网络推断 — GRN inference

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 047 | [047_rcistarget_tf_motif_network.R](03_virtual_perturbation/02_grn_inference/047_rcistarget_tf_motif_network.R) | RcisTarget motif/TF enrichment + regulatory network | gene list + motif DB → network/Sankey | R · RcisTarget, igraph, visNetwork | R | network, sankey |
| 🔴 | 081 | [081_pyscenic_regulon_tf_activity.py](03_virtual_perturbation/02_grn_inference/081_pyscenic_regulon_tf_activity.py) | pySCENIC GRN + ctx + AUCell wrapper | expr/loom → regulons + aucell | Py · pyscenic (GRNBoost) | Python | — (downstream UMAP/heatmap) |
| ✅ | 511 | [511_tf_convergence_depmap_jaspar](03_virtual_perturbation/02_grn_inference/511_tf_convergence_depmap_jaspar) | Three-evidence convergence to core TFs | tf_evidence.csv → convergence score | R · ggplot2, ggrepel | R | scatter, heatmap, lollipop |

### 03.04 · 药物扰动与响应 — Drug perturbation

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 070 | [070_chemcpa_drug_perturbation.py](03_virtual_perturbation/04_drug_perturbation/070_chemcpa_drug_perturbation.py) | chemCPA drug-perturbation expression prediction (GPU) | repo + config → train logs | Py · chemCPA, torch | Python | — |
| 🔴 | 071 | [071_scdrug_response_prediction.py](03_virtual_perturbation/04_drug_perturbation/071_scdrug_response_prediction.py) | scDrug single-cell drug response (heavy) | 10x/h5ad → cluster drug response | Py · scDrug, GDSC/PRISM | Python | — |
| 🟡 | 518 | [518_beyondcell_drug_response](03_virtual_perturbation/04_drug_perturbation/518_beyondcell_drug_response) | beyondcell core re-impl: BCS + therapeutic clusters | scRNA + drug signatures → BCS/ranking | R · UCell, ggplot2 | R | heatmap, lollipop, UMAP |

## 04 · 因果推断与遗传流行病 — Causal inference & genetics  (25)

### 04.01 · 工具变量准备 — Instrument preparation

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 028 | [028_mr_gwas_exposure_vcf_sig_snp.R](04_causal_inference_genetics/01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R) | Pick genome-wide significant SNPs from a VCF exposure | exposure vcf → significant SNPs | R · VariantAnnotation, gwasglue | R | — |
| 📄 | 029 | [029_mr_gwas_exposure_ld_clump.R](04_causal_inference_genetics/01_instrument_prep/029_mr_gwas_exposure_ld_clump.R) | LD-clump to independent instruments | candidate SNPs → independent IVs | R · gwasglue, TwoSampleMR | R | — |
| 📄 | 030 | [030_mr_gwas_exposure_add_eaf.R](04_causal_inference_genetics/01_instrument_prep/030_mr_gwas_exposure_add_eaf.R) | Add effect-allele frequency to exposure SNPs | SNPs → SNPs + EAF | R · ieugwasr | R | — |
| 📄 | 031 | [031_mr_gwas_exposure_weak_iv_filter.R](04_causal_inference_genetics/01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R) | F-statistic weak-instrument filter | SNPs → strong IVs (F≥10) | R · ieugwasr | R | — |

### 04.02 · 两样本孟德尔随机化 — Two-sample MR

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 032 | [032_mr_twosamplemr](04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr) | Self-contained MR causal inference (primary) | harmonized.csv → estimates + plots | R · ggplot2 (self-contained MR) | R | scatter, forest, funnel, leave-one-out |
| 📄 | 033 | [033_mr_gwas_finngen_outcome_backup.R](04_causal_inference_genetics/02_two_sample_mr/033_mr_gwas_finngen_outcome_backup.R) | GWAS exposure + FinnGen outcome MR template | exposure + FinnGen → MR table + plots | R · TwoSampleMR, qqman, RadialMR | R | scatter, forest, funnel, Manhattan, QQ, radial |
| 📄 | 043 | [043_MultimodalAD_MendelianRandomization_EN.R](04_causal_inference_genetics/02_two_sample_mr/043_MultimodalAD_MendelianRandomization_EN.R) | AD multi-omics MR main analysis | AD GWAS vcf → MR + sensitivity | R · TwoSampleMR, gwasglue | R | scatter, forest, funnel, Manhattan, QQ |
| 📄 | 055 | [055_immunecell_disease_mr_directionality.R](04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | Immune-cell ↔ disease bidirectional MR + Steiger | exposure + outcome → MR + Steiger | R · TwoSampleMR, RadialMR | R | — |
| ✅ | 519 | [519_local_mr_pipeline](04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline) | Fully local two-sample MR (no OpenGWAS API) | local exposure + outcome → estimates | R · TwoSampleMR, MRPRESSO, plinkbinr | R | scatter, forest, funnel, leave-one-out |

### 04.03 · cis-MR 与药靶 — cis-MR & drug targets

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 075 | [075_twosamplemr_coloc_drug_target.R](04_causal_inference_genetics/03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R) | MR + colocalization drug-target evidence chain | exposure/outcome/locus → MR + coloc | R · TwoSampleMR, coloc | R | — |
| ✅ | 535 | [535_mrbee_cis_mr](04_causal_inference_genetics/03_cis_mr_drug_target/535_mrbee_cis_mr) | MRBEE bias-corrected estimating-equation MR vs naive IVW | exposure/outcome GWAS summary → estimates | R · MRBEE, ggplot2 | R | lollipop, forest, scatter |
| 🟡 | 536 | [536_mrlink2_region_cis_mr](04_causal_inference_genetics/03_cis_mr_drug_target/536_mrlink2_region_cis_mr) | MR-link-2 single-region cis-MR (causal + pleiotropy) vs naive IVW | cis summary + LD → alpha/sigma_y + Type-I | Py · numpy, scipy, statsmodels; mrlink2 | Python | violin, forest, heatmap, scatter |

### 04.04 · 中介与多变量 MR — Mediation & MVMR

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 079 | [079_pqtl_mvmr_protein_mediation.R](04_causal_inference_genetics/04_mediation_mvmr/079_pqtl_mvmr_protein_mediation.R) | pQTL multivariable MR protein mediation | harmonised mvmr → MVMR mediation | R · (self-implemented MVMR) | R | — |
| ✅ | 499 | [499_lavaan_sem_mediation_path.R](04_causal_inference_genetics/04_mediation_mvmr/499_lavaan_sem_mediation_path.R) | SEM / path mediation with standardized-β diagram | composite scores → fit + path diagram | R · lavaan, semPlot | R | path-diagram |
| ✅ | 508 | [508_twostep_mediation_mr](04_causal_inference_genetics/04_mediation_mvmr/508_twostep_mediation_mr) | Two-step network mediation MR (Sobel/Delta/MC) | x + m instruments → mediation table | R · ggplot2 | R | path-diagram, forest |
| ✅ | 534 | [534_mvmr_cml_constrained](04_causal_inference_genetics/04_mediation_mvmr/534_mvmr_cml_constrained) | Constrained-ML multivariable MR (MVMR-cML-DP) vs IVW baseline | multi-exposure GWAS summary → direct effects | R · MendelianRandomization, ggplot2 | R | forest, dumbbell, heatmap |

### 04.05 · 共定位 — Colocalization

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 537 | [537_sharepro_coloc](04_causal_inference_genetics/05_colocalization/537_sharepro_coloc) | SharePro effect-group colocalization vs classic single-causal coloc | two-region summary + LD → group shares | Py · numpy, scipy, pandas; SharePro (vendored) | Python | scatter (locuscompare), lollipop, heatmap, dumbbell |

### 04.06 · TWAS 与单细胞 eQTL — TWAS & sc-eQTL

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 036 | [036_onek1k_twas_homo_hetero_fit.R](04_causal_inference_genetics/06_twas_sceqtl/036_onek1k_twas_homo_hetero_fit.R) | sc-eQTL homo/hetero elastic-net fit (step 1) | genotype + expr + cov → coefficients | R · glmnet, data.table | R | — |
| 🔴 | 037 | [037_onek1k_twas_component_fit.R](04_causal_inference_genetics/06_twas_sceqtl/037_onek1k_twas_component_fit.R) | Component-based prediction fit (step 2) | components → residuals + coef | R · glmnet, data.table | R | — |
| 🔴 | 038 | [038_onek1k_twas_weight_preprocess.R](04_causal_inference_genetics/06_twas_sceqtl/038_onek1k_twas_weight_preprocess.R) | Merge two-step coefficients for weights | step coefs → merged coefs | R · data.table | R | — |
| 🔴 | 039 | [039_onek1k_twas_fusion_weights.R](04_causal_inference_genetics/06_twas_sceqtl/039_onek1k_twas_fusion_weights.R) | Build per-cell-type FUSION weights | merged coefs → .RDat weights | R · data.table | R | — |
| 🔴 | 040 | [040_FUSION_TWAS_targetC.R](04_causal_inference_genetics/06_twas_sceqtl/040_FUSION_TWAS_targetC.R) | FUSION TWAS association (targetC weights) | GWAS sumstats + weights → TWAS | R · plink2R, glmnet | R | — |
| 🔴 | 041 | [041_FUSION_TWAS_S_targetC.R](04_causal_inference_genetics/06_twas_sceqtl/041_FUSION_TWAS_S_targetC.R) | FUSION TWAS association (S_targetC weights) | GWAS sumstats + weights → TWAS | R · plink2R, glmnet | R | — |
| 🔴 | 042 | [042_FUSION_TWAS_S_allC.R](04_causal_inference_genetics/06_twas_sceqtl/042_FUSION_TWAS_S_allC.R) | FUSION TWAS association (S_allC shared weights) | GWAS sumstats + weights → TWAS | R · plink2R, glmnet | R | — |

### 04.07 · 稳健 MR 估计量 — Robust MR estimators

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 533 | [533_mrcare_winnerscurse_mr](04_causal_inference_genetics/07_robust_mr_methods/533_mrcare_winnerscurse_mr) | Winner's-curse-corrected MR (CARE/RIVW) vs naive baseline | two-sample MR summary → estimates + plots | R · TwoSampleMR (baseline); MRcare | R | lollipop, forest, scatter, dumbbell |

## 05 · 机器学习 — Machine learning  (18)

### 05.01 · 特征筛选 — Feature selection

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 012 | [012_lasso_feature_selection](05_machine_learning/01_feature_selection/012_lasso_feature_selection) | LASSO logistic feature selection | expr + candidates → selected genes + plot | R · glmnet | R | scatter (CV / coef path) |
| ✅ | 013 | [013_svm_rfe_feature_selection](05_machine_learning/01_feature_selection/013_svm_rfe_feature_selection) | SVM-RFE recursive elimination ranking | expr + candidates → ranking + subset | R · e1071 | R | scatter, lollipop |
| ✅ | 014 | [014_randomforest_feature_selection](05_machine_learning/01_feature_selection/014_randomforest_feature_selection) | Random-forest Gini importance ranking | expr + candidates → importance + plot | R · randomForest | R | lollipop, scatter (OOB) |
| ✅ | 015 | [015_ml_feature_intersection_venn_upset](05_machine_learning/01_feature_selection/015_ml_feature_intersection_venn_upset) | Intersect feature genes across methods | gene-list dir → intersection + plot | R · UpSetR | R | venn, upset |
| ✅ | 034 | [034_multi_ml_feature_selection](05_machine_learning/01_feature_selection/034_multi_ml_feature_selection) | caret 10-method models → AUC + consensus features | expr → AUC table + consensus | R · caret, pROC, UpSetR | R | ROC, lollipop, upset |
| ✅ | 035 | [035_ml_combination_intersection](05_machine_learning/01_feature_selection/035_ml_combination_intersection) | Rank method combinations by intersection size | method feature lists → combo table | R · UpSetR | R | lollipop, upset |
| 📄 | 045 | [045_multimodalad_ml_models.R](05_machine_learning/01_feature_selection/045_multimodalad_ml_models.R) | Multi-ML integrated modelling (RSF/LASSO/GBM/BART…) | train/test txt → models + plots | R · randomForestSRC, glmnet, gbm, BART, xgboost | R | ROC, heatmap, lollipop |
| 📄 | 059 | [059_dual_disease_15ml_175combos.R](05_machine_learning/01_feature_selection/059_dual_disease_15ml_175combos.R) | Two-disease 15-ML × 175-combo screen + modelling | train/test → models + AUC heatmap | R · randomForestSRC, glmnet, ComplexHeatmap, sva | R | heatmap (AUC), ROC |
| ✅ | 502 | [502_biomarker_triple_vote](05_machine_learning/01_feature_selection/502_biomarker_triple_vote) | Topology × correlation × Boruta triple-vote shortlist | expr + group + candidates → vote/consensus | R · igraph, Boruta, Hmisc | R | heatmap (vote), lollipop |
| ✅ | 554 | [554_rra_consensus_features](05_machine_learning/01_feature_selection/554_rra_consensus_features) | Robust Rank Aggregation consensus across feature-selection methods | method×gene rank table → consensus + stability | R · RobustRankAggreg, ComplexHeatmap | R | lollipop, heatmap, upset, raincloud |

### 05.02 · 分类模型 — Classification models

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 550 | [550_tabpfn_tabular_classifier](05_machine_learning/02_classification_models/550_tabpfn_tabular_classifier) | TabPFN foundation model vs LASSO/GBDT honest incremental eval | expr + label → AUROC table + verdict | Py · tabpfn, scikit-learn, scipy | Python | dot (CV AUROC), ROC, PR, calibration, heatmap (confusion), lollipop |

### 05.03 · 生存机器学习 — Survival ML

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 496 | [496_mime_101combo_prognostic.R](05_machine_learning/03_survival_ml/496_mime_101combo_prognostic.R) | Mime 10-algorithm × 101-combo prognostic signature (usage snippet) | survival lists + genes → C-index table | R · Mime1 | R | lollipop, heatmap (C-index) |
| ✅ | 551 | [551_aorsf_oblique_survival](05_machine_learning/03_survival_ml/551_aorsf_oblique_survival) | Oblique random survival forest (aorsf) vs CoxPH / standard RSF baseline | survival table → C-index + risk strata | R · aorsf, survival, randomForestSRC, timeROC | R | time-dependent ROC, lollipop, KM |
| ✅ | 552 | [552_survex_survshap_explain](05_machine_learning/03_survival_ml/552_survex_survshap_explain) | survex time-dependent SurvSHAP(t) / SurvLIME explanation vs global baseline | survival table → time-varying importance | R · survex, survival, ranger | R | time-curve, dumbbell, lollipop, heatmap |
| ✅ | 553 | [553_riskregression_dca_calibration](05_machine_learning/03_survival_ml/553_riskregression_dca_calibration) | Honest survival-model eval: time-AUC + calibration + DCA + Brier/IBS | survival table → 4-axis eval + plots | R · riskRegression, dcurves, survival, prodlim | R | calibration, DCA, time-dependent ROC, line (Brier), lollipop |

### 05.04 · 可解释性 — Interpretability

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 052 | [052_shap_interpretation](05_machine_learning/04_interpretability/052_shap_interpretation) | Train best model + SHAP interpretation | geneexp.csv → SHAP table + plots | R · caret, kernelshap, shapviz | R | beeswarm, dependence, waterfall, force, ROC |

### 05.05 · 不确定性量化 — Uncertainty quantification

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 555 | [555_conformal_prediction_uq](05_machine_learning/05_uncertainty/555_conformal_prediction_uq) | Conformal prediction sets/intervals with finite-sample coverage vs naive baseline | table (target + features) → coverage diagnostics | Py · mapie, scikit-learn | Python | calibration, violin, line (efficiency), dumbbell |

### 05.06 · 泛化与外部验证 — Generalization & validation

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 503 | [503_generalization_robustness](05_machine_learning/06_generalization_validation/503_generalization_robustness) | Meta-analysis + LODO cross-cohort generalization | cohorts.rds → LODO/weight tables | R · metafor, glmnet, pROC | R | forest, lollipop, box |

## 06 · Bulk 组学 — Bulk omics  (19)

### 06.01 · 差异表达 — Differential expression

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 008 | [008_geo_expression_matrix](06_bulk_omics/01_differential_expression/008_geo_expression_matrix) | Annotate GEO probes → gene-level matrix | series_matrix + GPL → geneMatrix.csv | R base | R | — |
| ✅ | 009 | [009_geo_sample_grouping](06_bulk_omics/01_differential_expression/009_geo_sample_grouping) | Normalize matrix + attach group labels | geneMatrix + groups → labelled matrix | R · limma | R | — |
| ✅ | 010 | [010_geo_deg_volcano_heatmap_pca](06_bulk_omics/01_differential_expression/010_geo_deg_volcano_heatmap_pca) | limma two-group DE with three figures | expr matrix → DEG table + plots | R · limma, ComplexHeatmap | R | volcano, heatmap, PCA |
| ✅ | 056 | [056_geo_multicohort_batch_correction](06_bulk_omics/01_differential_expression/056_geo_multicohort_batch_correction) | Merge multi-cohort data + remove batch effect | cohort dir → corrected matrix + QC | R · limma/sva | R | PCA, box |

### 06.02 · 富集分析 — Enrichment

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 007 | [007_go_kegg_enrichment](06_bulk_omics/02_enrichment/007_go_kegg_enrichment) | GO/KEGG over-representation for a gene list | gene_list.csv → enrichment table + plots | R · clusterProfiler, org.Hs.eg.db, ggraph | R | dot, lollipop, network |
| ✅ | 546 | [546_enrichplot_emap_cnet_tree](06_bulk_omics/02_enrichment/546_enrichplot_emap_cnet_tree) | enrichGO + advanced plots replacing the plain bar (cnet/emap/tree) | gene_list.csv → enrichment table + plots | R · enrichplot, ggtangle, clusterProfiler, org.Hs.eg.db | R | dot, network (cnet/emap), tree, bar (baseline) |
| ✅ | 549 | [549_goplot_chord_enrichment](06_bulk_omics/02_enrichment/549_goplot_chord_enrichment) | GOplot gene×pathway chord / circle / membership heatmap | enrichment + logFC → relation matrix + plots | R · GOplot, ggplot2 | R | chord, circle, heatmap, lollipop, bar (baseline) |

### 06.03 · 共表达网络(WGCNA 家族) — Co-expression networks

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 054 | [054_wgcna_coexpression](06_bulk_omics/03_coexpression_networks/054_wgcna_coexpression) | Bulk WGCNA co-expression + module–trait | expr + traits → modules + figures | R · WGCNA, ComplexHeatmap | R | scale-free, dendrogram, module-trait heatmap |
| 🟡 | 504 | [504_hdwgcna_single_cell](06_bulk_omics/03_coexpression_networks/504_hdwgcna_single_cell) | hdWGCNA single-cell co-expression (metacell) | sc_counts.rds → modules + hubs | R · Seurat, hdWGCNA, igraph | R | soft-power, dendrogram, module feature-plot |
| ✅ | 538 | [538_netrep_module_preservation](06_bulk_omics/03_coexpression_networks/538_netrep_module_preservation) | NetRep permutation test of WGCNA module preservation across cohorts | discovery+test expr + modules → Zsummary/p | R · NetRep, ggplot2 | R | scatter, lollipop, density |
| ✅ | 539 | [539_smccnet_multiomics_network](06_bulk_omics/03_coexpression_networks/539_smccnet_multiomics_network) | SmCCNet trait-driven sparse multi-omics network vs unsupervised baseline | mRNA + miRNA + trait → subnetwork + hubs | R · SmCCNet, igraph, ggraph | R | network, heatmap, lollipop, dumbbell |
| 🟡 | 540 | [540_cwgcna_causal_module](06_bulk_omics/03_coexpression_networks/540_cwgcna_causal_module) | CWGCNA causal-direction (mediation) inference on WGCNA modules vs correlation baseline | expr + traits (driver) → causal directions | R · WGCNA (baseline); CWGCNA, ggraph | R | lollipop, dumbbell, network |

### 06.04 · 多组学整合与分型 — Multi-omics integration

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 083 | [083_mofa_diablo_multiomics.R](06_bulk_omics/04_multiomics_integration/083_mofa_diablo_multiomics.R) | MOFA2 / mixOmics DIABLO multi-omics latent integration | multi-view matrices → factors/heatmaps | R · MOFA2 (reticulate), mixOmics | R | factor plot, heatmap |
| ✅ | 084 | [084_nmf_consensus_clustering](06_bulk_omics/04_multiomics_integration/084_nmf_consensus_clustering) | NMF rank selection + consensus clustering subtyping | feature matrix → subtypes + figures | R · NMF, ConsensusClusterPlus, ComplexHeatmap | R | consensus, rank survey, subtype heatmap |

### 06.05 · 突变/甲基化/蛋白/代谢 — Mutation, methylation, proteome, metabolome

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 522 | [522_mutation_maftools_pipeline.R](06_bulk_omics/05_mutation_methylation_proteome/522_mutation_maftools_pipeline.R) | Somatic mutation summary template | MAF → oncoplot/summary | R · maftools | R | oncoplot |
| 📄 | 523 | [523_methylation_minfi_champ_pipeline.R](06_bulk_omics/05_mutation_methylation_proteome/523_methylation_minfi_champ_pipeline.R) | Methylation differential analysis template | beta + meta → M-value dist + heatmap | R · limma, minfi, ChAMP | R | heatmap |
| 📄 | 524 | [524_proteomics_limma_msstats_pipeline.R](06_bulk_omics/05_mutation_methylation_proteome/524_proteomics_limma_msstats_pipeline.R) | Proteomics differential analysis template | protein + meta → volcano + heatmap | R · limma, MSstats | R | volcano, heatmap |
| 📄 | 525 | [525_metabolomics_metaboanalystR_pipeline.R](06_bulk_omics/05_mutation_methylation_proteome/525_metabolomics_metaboanalystR_pipeline.R) | Metabolomics differential analysis template | metabolite + meta → volcano + heatmap | R · MetaboAnalystR | R | volcano, heatmap |
| 📄 | 526 | [526_cnv_gistic_or_cnvkit_pipeline.md](06_bulk_omics/05_mutation_methylation_proteome/526_cnv_gistic_or_cnvkit_pipeline.md) | CNV analysis entry note (GISTIC2/CNVkit/inferCNV) | — | — | md | — |

## 07 · 临床与转化 — Clinical & translational  (20)

### 07.01 · 诊断模型 — Diagnostic models

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 016 | [016_diagnostic_model_roc_calibration_dca](07_clinical_translational/01_diagnostic_models/016_diagnostic_model_roc_calibration_dca) | Logistic diagnostic model, full clinical evaluation | expr + genes → evaluation figures | R · rms, rmda, pROC | R | nomogram, calibration, DCA, ROC, forest, box |
| ✅ | 063 | [063_geo_diagnostic_validation](07_clinical_translational/01_diagnostic_models/063_geo_diagnostic_validation) | External-cohort validation of a diagnostic model | train + valid matrices → AUC + plot | R · rms, pROC | R | ROC, calibration |

### 07.02 · 预后与生存 — Prognosis & survival

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 048 | [048_tcga_single_gene_survival](07_clinical_translational/02_prognosis_survival/048_tcga_single_gene_survival) | Single-gene OS/DSS/DFI/PFI survival | gene_survival.csv → 4-endpoint KM | R · survival, survminer | R | KM (4 endpoints) |
| ✅ | 057 | [057_tcga_prognostic_risk_model](07_clinical_translational/02_prognosis_survival/057_tcga_prognostic_risk_model) | Prognostic risk model, five-figure panel | risk.csv → 5 figures + table | R · survival, timeROC, ComplexHeatmap | R | risk-plot, status, heatmap, KM, timeROC |
| ✅ | 060 | [060_tcga_immune_butterfly](07_clinical_translational/02_prognosis_survival/060_tcga_immune_butterfly) | Single-gene ↔ immune two-sided butterfly | expr + immune → correlation + plot | R · ggplot2 | R | butterfly (diverging) |

### 07.03 · 免疫浸润与解卷积 — Immune infiltration & deconvolution

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 017 | [017_immune_infiltration_source.R](07_clinical_translational/03_immune_infiltration/017_immune_infiltration_source.R) | CIBERSORT deconvolution engine (SVR source) | expr + LM22 → cell fractions | R · e1071 | R | — |
| 📄 | 018 | [018_immune_infiltration_scoring.R](07_clinical_translational/03_immune_infiltration/018_immune_infiltration_scoring.R) | Immune infiltration scoring + cell/function matrices | expr → score matrix + corr | R · e1071, preprocessCore | R | heatmap (corr) |
| ✅ | 021 | [021_immune_infiltration_viz](07_clinical_translational/03_immune_infiltration/021_immune_infiltration_viz) | Fraction-matrix difference/composition/correlation viz | CIBERSORT csv → 3 figures | R · ggpubr, ComplexHeatmap | R | box, stacked-bar, heatmap |
| 🟡 | 492 | [492_iobr_multimethod_deconvolution.R](07_clinical_translational/03_immune_infiltration/492_iobr_multimethod_deconvolution.R) | IOBR 7-method deconvolution + cross-method consistency | expr + group → merged matrix + plots | R · IOBR, tidyverse | R | stacked-bar, box, heatmap |
| ✅ | 520 | [520_bayesprism_deconvolution](07_clinical_translational/03_immune_infiltration/520_bayesprism_deconvolution) | BayesPrism Bayesian deconvolution with ground-truth check | scRNA ref + bulk → fractions + accuracy | R · BayesPrism | R | scatter, heatmap |

### 07.04 · 药物警戒 — Pharmacovigilance

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 078 | [078_faers_pharmacovigilance](07_clinical_translational/04_pharmacovigilance/078_faers_pharmacovigilance) | FAERS disproportionality (ROR/PRR/BCPNN/EBGM) | reports/2×2 counts → signals.csv | R · ggplot2 | R | forest, heatmap |

### 07.05 · 疾病负担与人群队列 — Disease burden & population cohorts

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 527 | [01_GBD/527_gbd_burden_trend](07_clinical_translational/05_epidemiology_burden/01_GBD/527_gbd_burden_trend) | GBD ASR trend + EAPC + Das Gupta decomposition + SDI | burden/pop/sdi csv → tables + figures | R · dplyr, ggplot2 | R | line-trend, forest, lollipop, diverging-lollipop, scatter |
| ✅ | 528 | [02_NHANES/528_nhanes_survey_weighted](07_clinical_translational/05_epidemiology_burden/02_NHANES/528_nhanes_survey_weighted) | NHANES survey-weighted means / regression / prevalence | nhanes.csv → svyglm + prevalence | R · survey, dplyr | R | dumbbell, forest, lollipop |
| ✅ | 529 | [03_CHARLS/529_charls_longitudinal_cohort](07_clinical_translational/05_epidemiology_burden/03_CHARLS/529_charls_longitudinal_cohort) | CHARLS trend + equipercentile equating + LMM + Cox/KM | panel.csv → trend/crosswalk/LMM/Cox | R · lme4, survival | R | line-trend, concordance, violin, forest, KM |
| ✅ | 530 | [04_comorbidity_network/530_comorbidity_network](07_clinical_translational/05_epidemiology_burden/04_comorbidity_network/530_comorbidity_network) | Disease-pair association → igraph → Louvain modules | patients.csv → network + metrics | R · igraph, ggraph | R | network, heatmap, lollipop |
| 🗃️ | — | [99_external_sources](07_clinical_translational/05_epidemiology_burden/99_external_sources) | GBD/NHANES/CHARLS 上游第三方源码树(git 忽略,仅本地参考) | — | — | R | — |
| 🗃️ | — | [comorbidity_paper_template_refs.ris](07_clinical_translational/05_epidemiology_burden/comorbidity_paper_template_refs.ris) | 共病选题的参考文献(仅本地) | — | — | — | — |
| 🗃️ | — | [literature_summary_comorbidity.md](07_clinical_translational/05_epidemiology_burden/literature_summary_comorbidity.md) | 共病文献综述草稿(仅本地) | — | — | — | — |
| 🗃️ | — | [sources_index.csv](07_clinical_translational/05_epidemiology_burden/sources_index.csv) | 疾病负担数据源索引(仅本地) | — | — | — | — |
| 🗃️ | — | [topic_candidates.md](07_clinical_translational/05_epidemiology_burden/topic_candidates.md) | 疾病负担选题候选(仅本地) | — | — | — | — |

## 08 · 结构与药物设计 — Structure & drug design  (5)

### 08.01 · 分子对接 — Molecular docking

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 022 | [022_docking_binding_energy_viz](08_structure_drug_design/01_docking/022_docking_binding_energy_viz) | Binding-energy heatmap + strongest-binding ranking | binding_energy.csv → figures | R · ComplexHeatmap | R | heatmap, lollipop |
| ✅ | 547 | [547_prolif_interaction_fingerprint](08_structure_drug_design/01_docking/547_prolif_interaction_fingerprint) | ProLIF protein-ligand interaction fingerprint + residue occupancy | pose/traj → fingerprint + occupancy | Py · prolif, MDAnalysis, rdkit | Python | barcode, heatmap, lollipop |
| ✅ | 556 | [556_posebusters_validity_panel](08_structure_drug_design/01_docking/556_posebusters_validity_panel) | PoseBusters physical-validity check panel for docking/AI poses | poses.sdf → check table + pass rates | Py · posebusters, rdkit, pandas | Python | heatmap (tick), lollipop, dumbbell |

### 08.02 · 分子动力学 — Molecular dynamics

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 086 | [086_vina_gromacs_mmpbsa_mdanalysis_pipeline.py](08_structure_drug_design/02_md_simulation/086_vina_gromacs_mmpbsa_mdanalysis_pipeline.py) | Vina docking + GROMACS MD + MM-PBSA pipeline | receptor/ligand → trajectory + ΔG | Py · Vina, GROMACS, gmx_MMPBSA, MDAnalysis | Python | scatter (RMSD/RMSF/Rg/SASA/energy) |
| ✅ | 548 | [548_bio3d_md_dccm_pca](08_structure_drug_design/02_md_simulation/548_bio3d_md_dccm_pca) | bio3d ensemble/MD: PCA + DCCM + RMSF with collectivity null | ensemble/PDB/traj → dynamics tables | R · bio3d, ggplot2 | R | heatmap (DCCM), scatter (PCA), vector-field (porcupine), lollipop, dumbbell |

## 09 · 网络药理学 — Network pharmacology  (8)

### 09.01 · 靶点数据库提取 — Target databases

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 001 | [001_ctd_compound_targets](09_network_pharmacology/01_target_databases/001_ctd_compound_targets) | Extract & de-dup compound targets from a CTD export | CTD csv → targets.csv | R base | R | — |
| ✅ | 002 | [002_swisstarget_compound_targets](09_network_pharmacology/01_target_databases/002_swisstarget_compound_targets) | Extract compound targets from a SwissTargetPrediction export | Swiss csv → targets.csv | R base | R | — |
| ✅ | 004 | [004_genecards_disease_targets](09_network_pharmacology/01_target_databases/004_genecards_disease_targets) | Extract disease targets from a GeneCards export | GeneCards csv → targets.csv | R base | R | — |

### 09.02 · 靶点交集与集合图 — Target intersection

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 003 | [003_ctd_swiss_target_union_venn](09_network_pharmacology/02_target_intersection/003_ctd_swiss_target_union_venn) | Union/intersection of CTD vs Swiss compound targets | target csvs → set table + plot | R · UpSetR | R | venn, upset, lollipop |
| ✅ | 005 | [005_omim_genecards_target_venn](09_network_pharmacology/02_target_intersection/005_omim_genecards_target_venn) | Union/intersection of OMIM vs GeneCards disease targets | target csvs → set table + plot | R · UpSetR | R | venn, upset, lollipop |
| ✅ | 006 | [006_disease_compound_target_venn](09_network_pharmacology/02_target_intersection/006_disease_compound_target_venn) | Disease ∩ compound targets → core targets | target csvs → intersection + plot | R · UpSetR | R | venn, lollipop |
| ✅ | 011 | [011_deg_drug_target_intersection](09_network_pharmacology/02_target_intersection/011_deg_drug_target_intersection) | Multi-set DEG ∩ drug ∩ disease target intersection | gene/target csvs → intersection + plot | R · UpSetR | R | venn, upset, lollipop |

### 09.03 · 成药性评分 — Druggability

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 493 | [493_opentargets_dgidb_chembl_druggability.py](09_network_pharmacology/03_druggability/493_opentargets_dgidb_chembl_druggability.py) | Composite druggability score for a target set (live APIs) | HGNC genes → score DataFrame | Py · requests, pandas, mygene | Python | — |

## 10 · 可视化 — Visualization  (13)

### 10.01 · 高级图型 — Advanced plot types

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 053 | [053_circlize_chromosome_circos](10_visualization/01_advanced_plots/053_circlize_chromosome_circos) | Plot genes onto chromosomes (genomic distribution) | gene_positions.csv → circos | R · circlize | R | circos |
| ✅ | 498 | [498_ggalluvial_sankey](10_visualization/01_advanced_plots/498_ggalluvial_sankey) | Multi-layer alluvial/Sankey flow | long table → alluvial figure | R · ggalluvial | R | sankey, alluvial |
| ✅ | 512 | [512_raincloud_plot](10_visualization/01_advanced_plots/512_raincloud_plot) | Raincloud (half-violin + box + jitter) vs bar charts | data.csv → stats + raincloud | R · ggdist | R | raincloud |
| ✅ | 513 | [513_ridgeline_plot](10_visualization/01_advanced_plots/513_ridgeline_plot) | Ridgeline distribution over an ordered factor | data.csv → summary + ridgeline | R · ggridges | R | ridgeline |
| ✅ | 514 | [514_dumbbell_slope_plot](10_visualization/01_advanced_plots/514_dumbbell_slope_plot) | Dumbbell + slope for paired change | data.csv → paired change + plots | R · ggrepel | R | dumbbell, slopegraph |
| ✅ | 515 | [515_chord_diagram](10_visualization/01_advanced_plots/515_chord_diagram) | Chord diagram for directed relations/flows | matrix.csv → flows + chord | R · circlize | R | chord |
| ✅ | 516 | [516_composite_multipanel](10_visualization/01_advanced_plots/516_composite_multipanel) | "Figure 1" multi-panel composite template | self-contained → composite figure | R · patchwork, ggrepel | R | composite (volcano+heatmap+forest+UMAP) |
| ✅ | 532 | [532_scpubr_publication_figures](10_visualization/01_advanced_plots/532_scpubr_publication_figures) | SCpubr publication-grade, colorblind-safe single-cell figure set | Seurat object → standardized figure set | R · SCpubr, Seurat, ggplot2 | R | UMAP, dotplot, feature-map, violin, alluvial |

### 10.02 · 模板与外部资源 — Templates & external resources

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | — | [advanced_figure_tools.csv](10_visualization/02_templates_resources/advanced_figure_tools.csv) | 高级图型工具清单 | — | — | — | — |
| 🗃️ | — | [ai_scientific_figures](10_visualization/02_templates_resources/ai_scientific_figures) | AutoFigure-Edit(ICLR'26)本地参考:方法描述 → 可编辑 SVG 示意图 | — | — | Python | schematic |
| 📄 | — | [download_advanced_figure_tools.ps1](10_visualization/02_templates_resources/download_advanced_figure_tools.ps1) | 按清单批量拉取高级图型工具 | 清单 csv → 本地仓库 | PowerShell | PS | — |
| 📄 | — | [literature_download_links_for_fdm.txt](10_visualization/02_templates_resources/literature_download_links_for_fdm.txt) | 高级图型的文献下载链接 | — | — | — | — |
| 📄 | — | [templates](10_visualization/02_templates_resources/templates) | 出图模板 | — | R · ggplot2 | R | — |

---

## 图类型 → 模块 反查表

想要某种图,直接查它由哪些模块产出。由各模块的「图型」字段自动生成。

- **Volcano · 火山图** → 010, 516, 524, 525, 558, 559
- **Heatmap · 热图** → 010, 018, 021, 022, 026, 045, 046, 049, 054, 057, 059, 076, 077, 078, 081, 083, 084, 492, 496, 502, 506, 510, 511, 516, 518, 520, 523, 524, 525, 530, 531, 534, 536, 537, 539, 543, 544, 545, 547, 548, 549, 550, 552, 554, 556, 559, 560
- **ROC / PR** → 016, 034, 045, 052, 057, 059, 063, 550, 551, 553
- **Calibration / DCA / nomogram · 校准与决策曲线** → 016, 063, 550, 553, 555
- **Forest · 森林图** → 016, 032, 033, 043, 078, 503, 508, 516, 519, 527, 528, 529, 533, 534, 535, 536
- **KM / survival · 生存曲线** → 048, 057, 497, 529, 551
- **UMAP / tSNE · 降维嵌入** → 026, 027, 044, 046, 049, 050, 058, 062, 081, 506, 507, 516, 518, 532, 541, 560
- **Violin · 小提琴** → 026, 027, 044, 046, 049, 050, 509, 529, 532, 536, 542, 544, 545, 555, 558
- **Raincloud · 云雨图** → 512, 554, 557, 559
- **Ridgeline · 山脊图** → 513
- **Box · 箱线图** → 016, 021, 056, 492, 503, 557
- **Dot / bubble · 点图气泡图** → 007, 026, 046, 049, 050, 051, 062, 082, 510, 531, 532, 546, 550, 557
- **Lollipop · 棒棒糖图(条形图替代)** → 003, 005, 006, 007, 011, 013, 014, 022, 034, 035, 045, 058, 496, 502, 503, 507, 509, 511, 518, 521, 527, 528, 530, 531, 533, 535, 537, 538, 539, 540, 541, 542, 543, 547, 548, 549, 550, 551, 552, 553, 554, 556, 557, 559, 560
- **Dumbbell / slopegraph · 哑铃图与斜率图** → 514, 528, 533, 534, 537, 539, 540, 544, 548, 552, 555, 556, 559
- **Venn / UpSet · 集合图** → 003, 005, 006, 011, 015, 034, 035, 509, 554
- **PCA** → 010, 026, 027, 056, 548
- **Scatter · 散点** → 012, 013, 014, 032, 033, 043, 086, 495, 506, 511, 519, 520, 521, 527, 533, 535, 536, 537, 538, 541, 542, 543, 544, 545, 548, 560
- **Network · 网络图** → 007, 047, 514, 530, 531, 539, 540, 546, 558
- **Chord / circos / alluvial · 弦图环图桑基** → 047, 051, 053, 498, 515, 532, 549
- **Trajectory / vector field · 轨迹与向量场** → 044, 049, 050, 062, 082, 517, 548
- **Spatial map · 空间分布图** → 027, 050, 073, 080, 505, 521, 541, 542, 543, 544, 545
- **Feature map · 基因表达投影** → 026, 027, 044, 046, 050, 532
- **Composite multi-panel · 多面板拼图** → 516

---

## 元信息待补

以下条目在目录中存在,但没有可靠的用途/输入输出记录 —— 用到时补,不臆造:

- `01_single_cell/01_pipeline_qc/562_mixhvg_hvg_selection`
- `01_single_cell/02_integration_batch/563_concord_contrastive_integration`
- `01_single_cell/02_integration_batch/564_scextract_prior_integration`
- `01_single_cell/02_integration_batch/565_scmultibench_integration_benchmark`
- `01_single_cell/03_annotation_typing/566_phispace_soft_annotation`
- `01_single_cell/05_differential_expression/567_glimes_mixed_effect_de`
- `01_single_cell/10_foundation_models/569_nicheformer_sc_spatial_fm`
- `01_single_cell/10_foundation_models/570_epiagent_scatac_fm`
- `01_single_cell/10_foundation_models/571_captain_rna_protein_fm`
- `01_single_cell/10_foundation_models/572_cellvq_discrete_cell_fm`
- `02_spatial_transcriptomics/01_pipeline_segmentation/573_proseg_cell_segmentation`
- `02_spatial_transcriptomics/02_domains_svg_stats/574_stair_spatial_integration`
- `02_spatial_transcriptomics/02_domains_svg_stats/575_scale_spatial_method`
- `02_spatial_transcriptomics/05_cell_communication/576_cellnest_spatial_ccc`
- `02_spatial_transcriptomics/05_cell_communication/577_spider_spatial_ccc`
- `03_virtual_perturbation/01_insilico_knockout/561_regvelo_grn_velocity`
