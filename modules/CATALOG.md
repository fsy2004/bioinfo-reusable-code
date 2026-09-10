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


## 01 · 单细胞分析 — Single-cell analysis  (33)

### 01.01 · 上游与质控 — Pipeline & QC

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 023 | [023_rds_structure_check.R](01_single_cell/01_pipeline_qc/023_rds_structure_check.R) | Inspect an RDS object's structure | RDS → structure print + csv | R base | R | — |
| 📄 | 024 | [024_scrna_rds_prep.R](01_single_cell/01_pipeline_qc/024_scrna_rds_prep.R) | Build a Seurat/RDS object from raw data | raw → Seurat object | R · Seurat, SingleR, celldex | R | — |
| 📄 | 025 | [025_scrna_data_prep.R](01_single_cell/01_pipeline_qc/025_scrna_data_prep.R) | Read 10x data into a Seurat object | 10x raw → Seurat object | R · Seurat, Matrix | R | — |
| 📄 | 061 | [061_scfocal_gui_input_prep.R](01_single_cell/01_pipeline_qc/061_scfocal_gui_input_prep.R) | Prepare scFOCAL GUI input + launch Shiny | RData + map csv → RDS + GUI | R · Seurat, scFOCAL, shiny | R | — (interactive) |
| 🟡 | 562 | [562_mixhvg_hvg_selection](01_single_cell/01_pipeline_qc/562_mixhvg_hvg_selection) | 混合多种 HVG 打分方法(按秩取max)选高变基因,自带 Seurat vst 基线 + ground-truth recall / ARI / silhouette 评测 | 输入 example_data/counts.csv(基因×细胞原始 count,首列基因名)+ 可选 cell_metadata.csv(细胞类型)+ 可选 ground_truth_hvg.csv(基因,Tier);输出 results/562_method_metrics.csv、562_selected_features.csv、562_summary.txt、562_sessionInfo.txt + 5 张 PDF/PNG(assets/ 存 PNG) | R 4.4.3 · Seurat 5.5.0, Matrix, ggplot2, ggrepel, scran, SingleCellExperiment, SummarizedExperiment, scuttle, mclust(可选,缺则 ARI=NA)· mixhvg(可选,装了才走 Step 6 官方一致性核对;本机未装) | R | 562_fig1_recall_lollipop(lollipop)、562_fig2_tier_slopegraph(slopegraph)、562_fig3_mean_variance_scatter(散点)、562_fig4_method_jaccard_heatmap(heatmap)、562_fig5_recall_vs_silhouette(散点);无条形图 |

### 01.02 · 整合与批次校正 — Integration & batch correction

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 506 | [506_scvi_scanvi_integration](01_single_cell/02_integration_batch/506_scvi_scanvi_integration) | scVI/scANVI integration + label transfer (vs PCA baseline) | h5ad (batch/label) → integration + labels | Py · scvi-tools, scanpy, sklearn | Python | UMAP, scatter, heatmap (confusion) |
| 🟡 | 563 | [563_concord_contrastive_integration](01_single_cell/02_integration_batch/563_concord_contrastive_integration) | 多批次单细胞整合模块:3 个本机可跑基线(PCA/批次中心化/ComBat+PCA)+ 守卫式 CONCORD 接口,用「批次混合熵 / 生物保真(kNN纯度+ARI/NMI) / 全局几何(trustworthiness+成对距离 Spearman)」三类指标同台对照评估。 | 输入:example_data/counts.csv(细胞×基因 raw count)+ cell_meta.csv(batch/cell_type)+ 可选 true_geometry.csv(合成数据的无批次参照系);或 --h5ad。输出:results/563_integration_metrics.csv、563_percell_batch_entropy.csv、563_summary.json、各方法 embedding_*.npy;assets/ 4 组 PNG+PDF。 | Python 3.12 · scanpy 1.12.1, anndata 0.12.14, scikit-learn 1.8.0, umap-learn 0.5.12, numpy, pandas, scipy, matplotlib(基线全部本机已装);可选 concord-sc(上游 v1.0.13,MIT;torch 需按 CUDA 自行先装,不在其 install_requires 内)——本机未装,已实测 import concord 报 ModuleNotFoundError。 | Python | fig1 UMAP 散点矩阵(方法×批次/细胞类型上色)、fig2 指标热图(格内原值+列内归一化配色)、fig3 批次混合 vs 生物保真权衡散点、fig4 每细胞批次混合熵 violin+抖动散点(raincloud 风格)。全程无条形图。 |
| 🟡 | 564 | [564_scextract_prior_integration](01_single_cell/02_integration_batch/564_scextract_prior_integration) | 多批次 scRNA 整合评测:以「批次混合熵 × 细胞类型 kNN 保真度 × 稀有类型保真度」双轴对比未校正 PCA / ComBat(/Harmony),并守卫式封装 scExtract 的 scanorama_prior / cellhint_prior 先验整合 | 输入 example_data/synthetic_3batch.h5ad(AnnData,1890 细胞 × 300 基因,需 obs['batch'] + obs['cell_type'],缺失自动按 seed=0 重建);输出 results/564_baseline_metrics.csv、results/564_summary.json、4 组 PDF+PNG 图,展示图复制到 assets/ | 基线(本机已具备):scanpy 1.12.1, anndata 0.12.14, scikit-learn 1.8.0, umap-learn 0.5.12, matplotlib, pandas, numpy;可选 harmonypy(装了自动多一个 Harmony 对照);完整方法需另装 scextract v0.2.0 + scanorama_prior / cellhint_prior fork(本机未装,守卫式跳过) | Python | fig1_umap_batch_vs_celltype(UMAP 散点网格)、fig2_mixing_vs_purity_tradeoff(权衡散点,点大小=稀有类型保真度)、fig3_shift_from_baseline(dumbbell)、fig4_crossbatch_celltype_similarity(热图)——全程无条形图 |
| 🟡 | 565 | [565_scmultibench_integration_benchmark](01_single_cell/02_integration_batch/565_scmultibench_integration_benchmark) | 把 scMultiBench 的 scIB 评测层封装成模块:给任意整合 embedding 打生物保留/批次校正/综合分,强制与朴素 PCA 基线对比,出热图+权衡散点+lollipop 排名。 | 输入 example_data/{rna_counts,adt_counts,metadata}.csv(行=细胞;metadata 须含 celltype/batch),或 --emb 名称=path.csv 传入外部 embedding;输出 results/metric.csv、results/565_summary.json、3 张 PDF+PNG(同步复制到 assets/)。 | 必需(本机全有):numpy pandas scipy scikit-learn networkx matplotlib + modules/_framework/pubstyle.py。可选:scib(--use-scib 上游路径,本机未装)、harmonypy(装了自动加 Harmony 对比)。kBET 另需 R 的 kBET 包 + rpy2。 | Python 3.12 | fig1_metric_heatmap(指标热图 viridis 0–1)、fig2_bio_vs_batch_scatter(生物保留 vs 批次校正权衡散点)、fig3_overall_lollipop(综合分 lollipop + 朴素基线参考虚线)。无条形图。 |

### 01.03 · 注释与细胞分型 — Annotation & cell typing

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 044 | [044_multimodalad_scrna.R](01_single_cell/03_annotation_typing/044_multimodalad_scrna.R) | AD brain scRNA pipeline + Monocle pseudotime | GSE157827 → object + trajectory | R · Seurat, SingleR, monocle | R | violin, UMAP, trajectory, feature-map |
| ✅ | 046 | [046_scrna_publication_figures](01_single_cell/03_annotation_typing/046_scrna_publication_figures) | Standard Seurat flow → publication figures | counts.csv → object + figures | R · Seurat, ggplot2 | R | UMAP, dotplot, heatmap, feature-map, violin |
| 📄 | 049 | [049_scrna_manual_annot_cellchat_trajectory.R](01_single_cell/03_annotation_typing/049_scrna_manual_annot_cellchat_trajectory.R) | Manual annotation + CellChat + trajectory | Seurat + markers → annotation + figures | R · Seurat, CellChat, monocle3 | R | violin, UMAP, heatmap, dotplot, trajectory |
| 🟡 | 566 | [566_phispace_soft_annotation](01_single_cell/03_annotation_typing/566_phispace_soft_annotation) | PhiSpace 连续表型软注释：把 query 细胞投影到参考类型张成的表型空间，给出每细胞×每类型的连续得分而非硬标签；自带质心相关 + PCA 回归两条本机可跑基线与已知混合比例真值评估，PhiSpace 本体为守卫式封装。 | 输入 example_data/reference_counts.csv + reference_metadata.csv(cell,cell_type) + query_counts.csv + query_metadata.csv(可选 true_group/true_*)，矩阵为行=基因列=细胞；输出 results/soft_scores_centroid_corr.csv、soft_scores_pca_regression.csv、hard_labels_argmax.csv、score_vs_truth_metrics.csv(+--run_phispace 时 phispace_scores.csv)，assets/ 4 图 PNG+PDF | R 4.4.3 · ggplot2(经 _framework/theme_pub.R)；可选 PhiSpace + SingleCellExperiment + S4Vectors(上游 DESCRIPTION 要求 R>=4.5.0，本机装不上，走守卫路径) | R | fig1 软得分热图(viridis tile) · fig2 真值vs预测散点(按方法分面+identity线) · fig3 小提琴+抖动点(按真实群体) · fig4 棒棒糖图(方法×类型 Pearson r)；全程无条形图 |

### 01.04 · 组成与丰度差异 — Composition / differential abundance

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 557 | [557_sccomp_composition_da](01_single_cell/04_composition_da/557_sccomp_composition_da) | sccomp Bayesian beta-binomial cell-composition DA vs 3 baselines | composition counts → DA table + plots | R · sccomp, voomCLR, limma, ggbeeswarm | R | boxplot, lollipop, raincloud, dot-matrix |
| 🟡 | 558 | [558_milo_neighborhood_da](01_single_cell/04_composition_da/558_milo_neighborhood_da) | Milo KNN-neighborhood differential abundance vs discrete-cluster baseline | SCE (reducedDim + condition) → DA table | R · miloR (or BiocNeighbors baseline), igraph, ggbeeswarm | R | beeswarm, network, volcano, violin |

### 01.05 · 差异表达(含 pseudobulk) — Differential expression

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| ✅ | 559 | [559_muscat_pseudobulk_ds](01_single_cell/05_differential_expression/559_muscat_pseudobulk_ds) | muscat multi-sample pseudobulk differential-state vs cell-level baseline | SCE (.rds) → DS table + plots | R · muscat, SingleCellExperiment, edgeR, limma | R | MDS, volcano, heatmap, lollipop, dumbbell, raincloud |
| 🟡 | 567 | [567_glimes_mixed_effect_de](01_single_cell/05_differential_expression/567_glimes_mixed_effect_de) | 多供体单细胞原始 UMI 计数上的 Poisson-GLMM(供体随机截距)差异表达,与朴素细胞级 t 检验、pseudobulk 两条基线同台对比,量化供体伪重复造成的一类错误膨胀。 | 输入 example_data/counts.csv(基因×细胞原始 UMI)+ cell_meta.csv(cell/donor/condition/exp_batch)+ 可选 truth.csv(合成金标准);输出 results/de_results_all_methods.csv、metrics_summary.csv、sessionInfo.txt,以及 assets/ 下 6 张 PNG+PDF;--use-glimes 且官方包已装时额外产出 glimes_official_poisson.csv | R 4.4.3 · MASS(glmmPQL)· nlme · ggplot2 · stats/utils;均随 R 分发或库内通用依赖,基线路径零安装。可选官方包 GLIMES(GitHub only,依赖 Bioconductor SummarizedExperiment/edgeR/MAST),本机未装,走守卫式跳过 | R | fig1 空基因 p 值 ECDF 阶梯折线(一类错误校准)· fig2 ROC 折线(含 AUC)· fig3 naive vs GLMM −log10 p 散点(按基因真值着色)· fig4 哑铃图(top20 供体驱动空基因显著性塌陷)· fig5 单基因供体级小提琴+抖动点+菱形供体均值 · fig6 GLMM σ² vs 朴素显著性散点(x 对数轴)。全程无条形图 |

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
| 🟡 | 568 | [568_scprint_foundation_grn](01_single_cell/10_foundation_models/568_scprint_foundation_grn) | scPRINT 单细胞基础模型三任务(GRN 推断 / 去噪 / 细胞嵌入与标签预测)的本机可跑朴素基线 + 守卫式官方 API 封装,评估口径逐条对齐上游源码 | 输入 example_data/{counts.csv (900 细胞 × 120 基因原始计数), cell_meta.csv (cell_id/cell_type/organism_ontology_term_id), true_grn_edges.csv (77 条真值边)};输出 results/{568_grn_benchmark.csv, 568_grn_top_edges.csv, 568_denoise_benchmark.csv, 568_celltype_cv_macro_f1.csv, 568_summary.json} + assets/ 4 图 (PNG+PDF) | 基线: numpy pandas scipy scikit-learn matplotlib(本机已装,零改动即跑);官方路径(不代装): scprint + scdataloader + lamindb 0.76.3 + bionty 0.49.0 + torch 2.2.0 + lightning + GRnnData + BenGRN + scib-metrics,python 3.10,可选 triton 2.0.0.dev20221202 走 GPU | Python 3.12 | fig1 PCA 嵌入散点 + 文库大小 violin+抖动 + 每折 macro-F1 lollipop;fig2 ROC/PR 曲线 + AUPRC 相对随机的 lollipop;fig3 TF→target 得分 heatmap(真值边红圈叠加)+ 真/假边得分 violin;fig4 Poisson NLL 随 k 曲线 + 相对 raw 的 MSE 改善 lollipop。全程无条形图 |
| ✅ | 569 | [569_nicheformer_sc_spatial_fm](01_single_cell/10_foundation_models/569_nicheformer_sc_spatial_fm) | Nicheformer(单细胞+空间联合预训练基础模型)守卫式封装 + 本机可跑的线性对照基线:比较 intrinsic 表达 PCA 与 niche-aware(⊕空间 kNN 邻域均值 PCA)在 niche / cell-type 标签上的同折 CV macro-F1,并做解离参考→空间 query 的跨模态标签迁移地板值。 | 输入:4 个 CSV —— 空间切片表达矩阵(cells×genes,原始计数)+ 元数据(cell_id, x, y, niche, cell_type),解离 scRNA 参考表达矩阵 + 元数据(cell_id, cell_type);示例 700×60 与 500×60,合成数据。输出:results/569_cv_macro_f1_per_fold.csv(task×representation×fold×macro_f1 长表)、results/569_summary.json(参数/均值±SD/delta/跨模态迁移指标/nicheformer 路径状态/依赖版本快照),assets/ 3 张 PNG+PDF。 | 基线(本机已装,零安装即跑):numpy, pandas, scikit-learn, matplotlib + 库内 modules/_framework/pubstyle.py。官方路径(本模块不代为安装,守卫式跳过):nicheformer(源码 pip install -e .,requires-python>=3.9)+ torch>=2.5.1 + pytorch-lightning>=2.0.0 + merlin-dataloader + dask-cuda + GPU + Mendeley Data 官方 checkpoint;上游钉死 numpy==1.26.4/pandas==1.5.3,与本机 numpy 2.x/pandas 2.x 冲突。 | Python | fig1_tissue_and_embeddings(2×2:组织切片按 niche / 按 cell type 的空间散点 + intrinsic 与 niche-aware 两套嵌入的 PC1/PC2 散点);fig2_representation_slopegraph(每折 macro-F1 的 slopegraph,intrinsic→niche-aware,niche 与 cell type 两任务分色);fig3_niche_confusion_heatmap(niche 标签行归一化混淆矩阵热图,两种表征并排 + 共享 colorbar)。全部 PDF 矢量 + 300dpi PNG 双出,图中文字英文,无条形图。 |
| 🟡 | 570 | [570_epiagent_scatac_fm](01_single_cell/10_foundation_models/570_epiagent_scatac_fm) | scATAC 细胞×cCRE 矩阵 → TF-IDF+SVD(LSI) 基线做嵌入/聚类/细胞类型预测/填补/批次混合评估,并守卫式封装 EpiAgent 基础模型路径(仅环境探测,不臆造调用) | 输入: example_data/cell_by_ccre_counts.csv(行=细胞,列=cCRE,0/1 计数)+ cell_metadata.csv(cell_type/batch);输出: results/570_summary.json、570_baseline_embedding.csv、4 张 PNG(同步 assets/) | Python 3.12 · numpy / pandas / scikit-learn / matplotlib(本机已具备,基线零安装);可选 EpiAgent 路径需 epiagent(PyPI v0.0.3)+ torch + flash-attn>=2.5.7 + NVIDIA GPU | python | 570_embedding_scatter.png(双 panel 散点)· 570_baseline_metrics_lollipop.png(lollipop + 朴素参照空心点)· 570_celltype_confusion.png(行归一化热图)· 570_coverage_violin.png(violin + 抖动散点);无条形图 |
| 🟡 | 571 | [571_captain_rna_protein_fm](01_single_cell/10_foundation_models/571_captain_rna_protein_fm) | 配对 CITE-seq「RNA→表面蛋白」填补基准台：matched-gene 与 PCA+Ridge 两条防泄漏基线 + CAPTAIN 本体守卫式探测（含 Drive 占位符识别） | 输入 example_data/citeseq_rna_counts.csv（细胞×基因 UMI）+ citeseq_adt_counts.csv（细胞×蛋白 ADT）+ protein_gene_map.csv（protein,cognate_gene）；输出 results/571_baseline_metrics.csv、571_pred_pca_ridge.csv、571_observed_clr.csv、571_summary.json 与 assets/ 三张 PNG+PDF | numpy, pandas, scikit-learn, matplotlib（全部本机已装，零安装即跑）；CAPTAIN 本体路径仅探测不 import，需 GPU + 上游 conda 环境（torch==2.1.2 等硬钉版本，与本机不兼容，须独立环境） | Python 3.12 | 571_protein_r_dumbbell.png（dumbbell，每蛋白两基线留出 Pearson r）、571_obs_vs_pred_scatter.png（散点小多图 + y=x 参考线）、571_specificity_heatmap.png（预测×观测相关矩阵，对角占优检验）；无条形图 |
| 🟡 | 572 | [572_cellvq_discrete_cell_fm](01_single_cell/10_foundation_models/572_cellvq_discrete_cell_fm) | 离散「细胞词表」(VQ codebook)单细胞表征模块：PCA+k-means 码本基线量化离散化的信息损失/码本坍缩/重构 R²，并对官方 CellVQ 提供逐符号核实过的守卫式封装。 | 输入 example_data/synthetic_counts.csv（cells×genes 计数，首列 cell_id、第二列 cell_type，550×160，synthetic）；输出 results/572_codebook_sweep.csv、572_cell_codes.csv、572_summary.json + assets/ 四图（PNG+PDF） | numpy, pandas, scikit-learn, matplotlib（全部本机已有）；torch 仅用于守卫式体检的可选探测，缺失时优雅降级 | Python | 572_codebook_sweep（折线+点三联，含连续 PCA 基线虚线）、572_code_usage_lollipop（lollipop）、572_code_celltype_heatmap（heatmap, viridis）、572_embedding_vs_code（散点双联）；无条形图 |

## 02 · 空间转录组 — Spatial transcriptomics  (24)

### 02.01 · 上游与细胞分割 — Pipeline & segmentation

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 📄 | 027 | [027_spatial_seurat_auto.R](02_spatial_transcriptomics/01_pipeline_segmentation/027_spatial_seurat_auto.R) | Spatial transcriptomics read/cluster/visualize | Visium h5 → spatial figures | R · Seurat, SingleR, glmGamPoi | R | PCA, violin, UMAP, feature-map, niche-map |
| 📄 | 050 | [050_spatial_cluster_annot_trajectory.R](02_spatial_transcriptomics/01_pipeline_segmentation/050_spatial_cluster_annot_trajectory.R) | Spatial cluster annotation + monocle3 pseudotime | Visium → spatial annotation + trajectory | R · Seurat, monocle3, patchwork | R | niche-map, violin, UMAP, feature-map, pseudotime |
| 🟡 | 573 | [573_proseg_cell_segmentation](02_spatial_transcriptomics/01_pipeline_segmentation/573_proseg_cell_segmentation) | 成像空间转录组转录本点云的细胞分割:自带「最近核外扩」可跑基线(半径扫描 + recall/precision/ambient-leak/ARI 评分),Proseg 本体为 Rust CLI 的守卫式命令行封装。 | 输入 example_data/transcripts.csv(合成,14,943 转录本 / 120 细胞 / 24 基因;必需列 x_location, y_location, feature_name, cell_id;基线评分另需 true_cell_id / true_cell_type)→ 输出 results/{baseline_radius_sweep.csv, baseline_cell_by_gene_counts.csv, baseline_cell_metadata.csv, 573_summary.json} + assets/ 5 图(PNG 300dpi + 矢量 PDF) | numpy, pandas, scipy, scikit-learn, matplotlib(均本机已装);可选外部 Rust 二进制 proseg 3.2.0(cargo install proseg / conda install -c bioconda rust-proseg),不在 PATH 时干净跳过 | Python 3.12 | fig1_segmentation_map(双 panel 空间散点,truth vs 基线)、fig2_radius_sweep(折线+点)、fig3_counts_per_cell(violin+抖动点)、fig4_count_recovery(散点+y=x 参考线)、fig5_celltype_heatmap(热图)—— 无条形图 |

### 02.02 · 空间域、空间可变基因与空间统计 — Domains, SVG & spatial statistics

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 505 | [505_spatial_advanced](02_spatial_transcriptomics/02_domains_svg_stats/505_spatial_advanced) | Spatial advanced: RCTD deconv + NMF niche + interface degree | sc ref + spatial rds → fractions/niche | R · spacexr, RcppML, mistyR | R | niche-map (×3) |
| ✅ | 541 | [541_banksy_spatial_domains](02_spatial_transcriptomics/02_domains_svg_stats/541_banksy_spatial_domains) | BANKSY neighbor-augmented spatial-domain segmentation vs non-spatial baseline | spatial csv → domains + ARI | R · Banksy, SpatialExperiment, aricode | R | spatial-scatter, lollipop, UMAP |
| ✅ | 542 | [542_nnsvg_spatial_svg](02_spatial_transcriptomics/02_domains_svg_stats/542_nnsvg_spatial_svg) | nnSVG spatially-variable genes (NNGP) vs non-spatial HVG baseline | counts + coords → SVG ranking | R · nnSVG, SpatialExperiment, scran | R | spatial-scatter, lollipop, scatter, violin |
| ✅ | 543 | [543_squidpy_spatial_statistics](02_spatial_transcriptomics/02_domains_svg_stats/543_squidpy_spatial_statistics) | squidpy spatial stats (Moran / nhood enrichment / co-occurrence / Ripley) | h5ad (spatial) → stats tables + plots | Py · squidpy, anndata, scanpy | Python | heatmap, lollipop, scatter, spatial-scatter |
| 🟡 | 574 | [574_stair_spatial_integration](02_spatial_transcriptomics/02_domains_svg_stats/574_stair_spatial_integration) | 多切片空间转录组整合模块：本机可跑的三级整合阶梯基线（PCA / ComBat+PCA / 空间平滑+ComBat+PCA，scIB 式双轴评分），外加签名逐行核对自上游源码的 STAIR HGAT 守卫式封装 | 输入 example_data/slices_expression.csv（spot×基因 raw count）+ slices_meta.csv（spot/slice/x/y/[true_domain]）；输出 results/574_baseline_metrics.csv、results/574_summary.json、assets/fig1-4（空间域散点矩阵 / 整合权衡散点 / 指标棒棒糖 / UMAP 批次着色，PNG+PDF） | Python 3.12：numpy pandas scipy scikit-learn matplotlib anndata scanpy umap-learn（均本机已装，基线零安装即跑）；可选 STAIR-tools 1.3.1 + torch 2.6.0 + torch_geometric + CUDA GPU（本机未装，守卫跳过）；clustering='mclust' 另需 rpy2 + R mclust | Python | fig1_spatial_domains（空间散点矩阵，行=方法含 Ground truth，列=切片）；fig2_integration_tradeoff（权衡散点，x=批次混合熵 y=ARI/NMI 均值）；fig3_metric_lollipop（棒棒糖图，四指标）；fig4_umap_by_slice（UMAP 按切片着色）。全部无条形图，色盲安全配色，统一 pubstyle |
| 🟡 | 575 | [575_scale_spatial_method](02_spatial_transcriptomics/02_domains_svg_stats/575_scale_spatial_method) | 空间组学多尺度空间域识别:空间平滑×Leiden分辨率网格上的跨种子稳定性搜索(朴素基线)+ SCALE 上游守卫式封装,含无空间信息对照。 | 输入 csv(行=spot,列 x,y,[domain_coarse,domain_fine],其余为基因);输出 results/ 下 stability_grid.csv、n_clusters_grid.csv、spatial_vs_nonspatial_ari.csv、baseline_domain_labels.csv、575_summary.json + 3 张 PNG/PDF | numpy, pandas, scikit-learn, scanpy, anndata, leidenalg, igraph, matplotlib(均本机已装);可选上游 scale + torch-geometric(未装,走守卫路径) | Python | 575_stability_grid (heatmap)、575_domains_multiscale (多面板空间散点)、575_spatial_gain_dumbbell (dumbbell);无条形图 |

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
| 🟡 | 576 | [576_cellnest_spatial_ccc](02_spatial_transcriptomics/05_cell_communication/576_cellnest_spatial_ccc) | 空间转录组细胞通讯:空间受限 LR 共表达乘积基线(本机 CPU 即跑)+ CellNEST(GATv2 图注意力)守卫式封装,输出与 CellNEST 官方 9 列 schema 对齐 | 输入 example_data/spatial_counts.csv(spot×gene 原始 counts)+ spatial_coordinates.csv(barcode,x,y[,region])+ lr_pairs.csv(Ligand,Receptor,Annotation,Reference);输出 results/ 下 baseline_top_ccc.csv、baseline_top_ccc_cellnest_schema.csv(CellNEST 9 列 headerless)、baseline_lr_summary.csv、baseline_permutation_control.csv、session_info.txt + 4 张 PNG/PDF;可选 --cellnest-csv 读入真实 CellNEST *_top20percent.csv 并出 baseline_vs_cellnest_concordance.csv | numpy, pandas, scipy, scikit-learn, matplotlib(本机全有,基线零安装);modules/_framework/pubstyle.py;CellNEST 本体需 Linux + NVIDIA GPU + singularity,本模块不 import 它 | Python | fig1_spatial_ccc_map.png(空间散点 + LineCollection 通讯边,viridis)、fig2_lr_lollipop.png(lollipop)、fig3_lr_component_heatmap.png(heatmap)、fig4_permutation_control.png(dumbbell + 误差棒);无条形图 |
| 🟡 | 577 | [577_spider_spatial_ccc](02_spatial_transcriptomics/05_cell_communication/577_spider_spatial_ccc) | 在相邻 spot 之间构建 interface(容量约束 Delaunay 图),对配体-受体对做 interface 层面 Moran's I 置换检验找空间可变互作(SVI),并与不建 interface 的 spot-level 共表达基线并排对照;官方 spider-st 包为守卫式可选路径。 | 输入 example_data/{expression.csv, spot_meta.csv(spot/x/y/cell_type), lr_pairs.csv(ligand/receptor), ground_truth.json(可选)};输出 results/{interfaces.csv, svi_results.csv(每 LR 对 spot/interface 两路的 Moran's I + p + BH-FDR), 577_summary.json, versions.txt} + assets/ 5 组 PNG+PDF。官方路径需 --adata 的 .h5ad。 | Python 3.12 · numpy / pandas / scipy(Delaunay) / scikit-learn(NearestNeighbors) / matplotlib;本机全部已装,零改动可跑。可选官方路径:spider-st(PyPI,python 3.8 环境 + somoclu/fa2/scgco + R 端 nnSVG/SPARK),本机未装,走守卫分支。 | Python | fig1 interface 空间网络散点(LineCollection 边着色);fig2 dumbbell(spot-level vs interface-level Moran's I,实心=FDR<0.05);fig3 lollipop(SVI 显著性排序);fig4 heatmap(LR 对 × 细胞类型对,行 z-score);fig5 散点对照图(top SVI vs null 对的空间分布)。无条形图。 |

### 02.06 · 空间多组学 — Spatial multi-omics

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 521 | [521_spatialglue_multiomics](02_spatial_transcriptomics/06_spatial_multiomics/521_spatialglue_multiomics) | SpatialGlue spatial multi-omics domains (GNN; baseline local) | RNA + ADT grid → ARI + domains | Py · sklearn (baseline); SpatialGlue, torch-geometric | Python | spatial-scatter, lollipop |
| 🟡 | 578 | [578_spatialex_omics_translation](02_spatial_transcriptomics/06_spatial_multiomics/578_spatialex_omics_translation) | 用 H&E 形态学做锚,把一张切片测到的组学 panel 跨切片翻译到另一张切片(SpatialEx/SpatialEx+ 的 panel 对角整合),自带 Ridge+空间平滑的可跑基线与均值地板对照 | 输入 example_data/slice{1,2}_{coords,he,panelA,panelB}.csv(细胞×坐标 / H&E embedding / 两个 log1p 表达 panel,合成);输出 results/578_metrics.csv、578_gene_pcc_{panelA_slice1_to_slice2,panelB_slice2_to_slice1}.csv、578_summary.json(含参数/种子/session_info),assets/ 4×PNG+PDF | Python 3.12 · numpy pandas scikit-learn matplotlib(本机已具备,零安装即跑);可选上游 SpatialEx(PyPI, MIT)+ torch + GPU + UNI 权重 | Python | 4 张:578_gene_pcc_dumbbell(dumbbell)、578_pred_vs_measured_scatter(散点)、578_spatial_translation(空间散点三联)、578_coexpression_heatmap(heatmap);无条形图 |
| 🟡 | 579 | [579_simo_spatial_multiomics](02_spatial_transcriptomics/06_spatial_multiomics/579_simo_spatial_multiomics) | 579 · SIMO —— 把无空间坐标的单细胞多组学(scRNA + 非转录组模态)通过最优传输概率性映射到空间转录组切片。模块并排跑三条路线:A 贪心相关性朴素基线(地板对照)、B 自写 POT fused-Gromov-Wasserstein 传输参照、C SIMO 正牌路线(守卫式,未装 simo-omics 即优雅退出打印安装命令,不静默降级)。合成数据自带 layer ground truth,四项外部指标(层准确率/中位位移/spot 占用率/单 spot 最大堆叠)全部落盘。 | 输入 example_data/ 五个 CSV:st_expression.csv(spot×gene counts)、st_coords.csv(spot,x,y,layer)、sc_rna_expression.csv + sc_rna_meta.csv(cell,layer)、sc_mod2_gene_activity.csv + sc_mod2_meta.csv([0,1] 活性分)。输出 results/:579_summary.json(全指标+依赖版本快照)、modality1_cell_to_spot.csv、modality2_cell_to_spot.csv、4 组 PNG+PDF;PNG 同步复制到 assets/。 | Python 3.12.5 · numpy 2.4.6 · pandas 2.3.3 · scipy 1.16.3 · scikit-learn 1.8.0 · POT 0.9.4 · matplotlib 3.11.0 · 库内 _framework/pubstyle;C 路线另需 simo-omics(PyPI 名 simo-omics,import 名 simo,上游要求 python=3.8;本机未安装,走守卫路径) | Python | fig1_spatial_mapping(散点三联:ground-truth / 贪心 / OT)、fig2_accuracy_slopegraph(每层准确率 slopegraph,双模态)、fig3_transport_heatmap(细胞层×spot 层传输质量热图,行归一)、fig4_spot_occupancy_dumbbell(占用率与最大堆叠 dumbbell)。全部无条形图,符合库内绘图铁律。 |

### 02.07 · 空间基础模型 — Spatial foundation models

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 580 | [580_novae_spatial_fm](02_spatial_transcriptomics/07_foundation_models/580_novae_spatial_fm) | 580 · Novae — 多切片空间转录组的空间域/niche 划分与跨切片可迁移性评估。默认跑三级朴素基线阶梯(expression-only PCA+KMeans / 空间 kNN niche 平滑 / niche 平滑+逐切片 z-score),并用复刻自上游 novae/monitor/eval.py 的 FIDE / JSD / heuristic 三项指标打分,另加 ARI 与诊断列 ARI_celltype。Novae 本体走守卫式封装(--run-novae),包未装、基因名非真实 symbol 或拿不到 HuggingFace 权重时如实返回 skipped/failed,不伪造结果。 | 输入:细胞级 CSV(cell_id, slide, x, y[, domain_true, cell_type_true], 其余列一律当基因原始 count),示例 example_data/spatial_demo_cells.csv(合成,972 cells × 40 genes × 3 slides)。输出:results/580_metrics.csv、580_domain_labels.csv、580_summary.json(跑 --run-novae 时另有 580_novae_status.json);图见 assets/。 | Python。基线路径仅需 numpy / pandas / scikit-learn / matplotlib(+ 本库 _framework/pubstyle)。Novae 路径另需 novae(pyproject 要求 Python>=3.11,依赖 scanpy>=1.9.8, anndata>=0.11.0, lightning>=2.2.1, torch>=2.2.1, torch-geometric>=2.5.2, huggingface-hub>=0.32.0, safetensors>=0.4.3, igraph>=0.11.8)+ anndata;本机未安装,Novae 路径未做端到端验证。 | Python | 3 张(均非条形图):580_domain_maps(空间散点矩阵,行=方法含 ground truth、列=切片)、580_metric_comparison(棒棒糖四联 FIDE/heuristic/ARI/JSD)、580_slide_composition(域×切片占比热图)。各含 .png + .pdf。 |

## 03 · 虚拟扰动技术 — Virtual perturbation  (26)

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
| ✅ | 561 | [561_regvelo_grn_velocity](03_virtual_perturbation/01_insilico_knockout/561_regvelo_grn_velocity) | GRN 约束的 RNA 速率 + 逐 TF in-silico 调控子敲除,经 CellRank 命运概率重分配打分筛查转录因子 | 输入 example_data/synthetic_velocity.h5ad(layers spliced/unspliced, obs cell_state/true_time)+ synthetic_prior_grn.csv(行=靶基因,列=调控因子,0/1);输出 results/tf_ko_fate_shift.csv、results/fate_probabilities_wt.csv、results/561_summary.json,及 assets/ 下 4 组 PNG+PDF | Python≥3.10 · numpy pandas scipy scikit-learn matplotlib anndata scanpy scvelo cellrank(2.3.2) statsmodels;官方路径另需 regvelo(拉 scvi-tools>=1.0.0,<1.2.1、torchode>=0.1.6、cellrank>=2.0.0)+ CUDA GPU,未装时守卫式跳过 | Python | fig1_grn_velocity_stream(扩散图上的速度流场)、fig2_tf_ko_fate_dumbbell(哑铃图,红点 FDR<0.05)、fig3_ko_effect_heatmap(TF×终末状态 AUROC 热图,★标 FDR<0.05)、fig4_pseudotime_vs_truth(双 panel 散点,分谱系 Spearman rho);全程无条形图 |
| 🟡 | 581 | [581_veloagent_velocity](03_virtual_perturbation/01_insilico_knockout/581_veloagent_velocity) | 空间信息驱动的 RNA velocity 与 in-silico 敲除:scVelo + 空间 kNN 平滑基线 + veloAgent 守卫式封装,出速度场 quiver / raincloud / lollipop | 输入 example_data/{spliced,unspliced}.csv(行=细胞 列=基因)+ spatial_meta.csv(cell,x,y,true_time,cluster)+ 可选 gene_truth.csv(gene,true_direction);输出 results/{581_summary.json, 581_knockout_scores.csv, 581_percell_consistency.csv} + assets/ 三张 PNG/PDF | 基线(本机已有):numpy pandas scipy scikit-learn matplotlib anndata scanpy scvelo;上游 veloagent 需独立 conda env(python>=3.11.8, mesa==2.1.5, veloae@git+VeloAE, PyTorch 单独装)+ STRING DB 三件套,本机未装 | Python | 581_spatial_velocity_field.png/.pdf(空间散点+quiver 速度场双 panel,点色=与真值方向余弦)、581_consistency_raincloud.png/.pdf(violin+box+jitter raincloud,两基线逐细胞方向准确度)、581_knockout_lollipop.png/.pdf(lollipop Top18,点色=mean|velocity|);无条形图 |

### 03.02 · 基因调控网络推断 — GRN inference

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 047 | [047_rcistarget_tf_motif_network.R](03_virtual_perturbation/02_grn_inference/047_rcistarget_tf_motif_network.R) | RcisTarget motif/TF enrichment + regulatory network | gene list + motif DB → network/Sankey | R · RcisTarget, igraph, visNetwork | R | network, sankey |
| 🔴 | 081 | [081_pyscenic_regulon_tf_activity.py](03_virtual_perturbation/02_grn_inference/081_pyscenic_regulon_tf_activity.py) | pySCENIC GRN + ctx + AUCell wrapper | expr/loom → regulons + aucell | Py · pyscenic (GRNBoost) | Python | — (downstream UMAP/heatmap) |
| ✅ | 511 | [511_tf_convergence_depmap_jaspar](03_virtual_perturbation/02_grn_inference/511_tf_convergence_depmap_jaspar) | Three-evidence convergence to core TFs | tf_evidence.csv → convergence score | R · ggplot2, ggrepel | R | scatter, heatmap, lollipop |
| 🟡 | 582 | [582_dspin_regulatory_network](03_virtual_perturbation/02_grn_inference/582_dspin_regulatory_network) | 从多重扰动 scRNA-seq 反推程序级 Ising 自旋网络(共享耦合矩阵 J + 每扰动场向量 h),自带混池相关 / 朴素平均场两条带真值评分的基线,D-SPIN 正式伪似然求解器走守卫式调用。 | 输入 example_data/expression.csv(行=细胞,列=基因,log 归一化)+ cell_meta.csv(必需列 sample_id/batch/if_control);可选 ground_truth_coupling.csv、ground_truth_response.csv。输出 results/ 下 baseline_correlation_network.csv、baseline_meanfield_network.csv、relative_response_vectors.csv、program_states.csv、program_top_genes.csv、582_summary.json;assets/ 下 6 张 PNG+PDF。 | 基线(本机全有):numpy pandas scipy scikit-learn matplotlib networkx。正式路径需 pip install dspin(上游 setup.py install_requires: anndata matplotlib scanpy tqdm igraph leidenalg — 这六个本机均已安装,仅 dspin 本体缺)。 | Python | 6 张(无条形图):fig1 并排耦合热图(真值 vs B0 vs B1)、fig2 边恢复散点(直接边 vs 间接对分色)、fig3 基线 dumbbell、fig4 扰动响应热图、fig5 响应 lollipop、fig6 力导向自旋网络图。均走 modules/_framework/pubstyle.py,PNG+PDF 双份。 |
| 🟡 | 583 | [583_kegni_knowledge_grn](03_virtual_perturbation/02_grn_inference/583_kegni_knowledge_grn) | 知识图增强 GRN 推断:5 种本机可跑基线(Pearson/Spearman/PCA 嵌入点积/纯知识先验/知识-表达秩融合)按 BEELINE 口径出 EPR/AUPRC/AUROC 榜单,并对上游 KEGNI 深度模型做守卫式 CLI 封装(不臆造 Python API)。 | 输入 example_data/expression.csv(基因×细胞,index_col=0)+ knowledge_graph.tsv(无表头 head/relation/tail)+ ground_truth_network.csv(表头 Gene1,Gene2);输出 results/benchmark_metrics.csv、edges_<method>.csv(Gene1,Gene2,EdgeWeight,可直接喂上游 eval.py)、predicted_edges_best.csv、583_summary.json(含 session_info + 守卫报告)、4 图 PNG+PDF(assets/ 存 PNG)。可选 --kegni-pred / --kegni-embedding 把上游结果并入同一张榜单。 | Python 3.12 · numpy 2.4.6 / pandas 2.3.3 / scipy 1.16.3 / scikit-learn 1.8.0 / networkx 3.6.1 / matplotlib 3.11.0(全部本机已有,基线路径零安装);上游 KEGNI 本体另需 torch + dgl + transformers(本机未装,故走守卫分支;上游无 setup.py/requirements.txt,版本上游从未给出) | Python | fig1_benchmark_dotplot.png(方法×指标 lollipop dot plot)、fig2_pr_curves.png(PR 曲线+随机基准)、fig3_evidence_complementarity.png(表达证据 vs 知识证据百分位秩散点)、fig4_knowledge_gain_slopegraph.png(纯表达→知识融合 slopegraph)。无条形图。 |
| 🟡 | 584 | [584_cellpolaris_grn_transfer](03_virtual_perturbation/02_grn_inference/584_cellpolaris_grn_transfer) | 在已有 GRN 上建高斯概率图模型做 TF 虚拟敲除(ΔX),并用 ΔX 与真实相邻状态表达差的余弦相似度沿分化轨迹排主控 TF;上游迁移学习生成 GRN 段为守卫式封装。 | 输入 example_data/: grn.txt(TF/TG/Score 制表符) + expr_metacell_<state>.csv(基因×metacell) + expr_pseudobulk_<state>.txt(基因<TAB>值,无表头) + trajectory.txt(状态制表符分隔);输出 results/: deltax_<state>.csv、pgm_params_<state>.csv、master_tf_scores.csv、584_summary.json;assets/: 3 张 PNG+PDF | Python 3.12 · numpy · pandas · torch(CPU) · matplotlib(本机全部已有,无需安装);GRN 迁移段另需 torch_geometric 2.3.1 + CUDA + sci-db PECA2 数据集 + 自训权重(未安装,守卫式检查) | Python | 584_fig1_pgm_fit.png(NLL 折线 + 耦合 k vs 边协方差散点)、584_fig2_deltax_heatmap.png(ΔX 发散色热图)、584_fig3_master_tf.png(主控 TF lollipop,点大小编码 ||观测ΔExpr||);无条形图 |
| ✅ | 585 | [585_ignite_grn_inference](03_virtual_perturbation/02_grn_inference/585_ignite_grn_inference) | 用非对称动力学 Ising 模型的反问题（IGNITE，PLoS Comput Biol 2026）从未扰动的伪时序单细胞数据反推有向有符号 GRN 并模拟基因敲除，自带 3 个本机可跑的朴素 GRN 基线（Pearson / GraphicalLassoCV 偏相关 / 滞后岭回归）做 AUROC-AUPRC 边恢复对照。 | 输入：`example_data/spins_pseudotime_ordered.csv`（行=基因，列=按伪时序排列的细胞，值 ±1 自旋；非 ±1 输入按每基因中位数自动二值化）+ 可选 `example_data/true_network.csv`（真值耦合矩阵 J[target, regulator]，仅评测用）。均为合成数据，文件头标注 synthetic。输出：`results/network_b1|b2|b3.csv`（各法耦合矩阵）、`results/edge_recovery_scores.csv`、`results/585_summary.json`、`results/ignite_J.csv` 与 `ignite_h.csv`（仅 --ignite-repo 成功时）、`assets/` 3 张图。参数：--spins --truth --ignite-repo --delta-t --lam --epochs --outdir --figdir --seed。 | Python 3.12；基线仅需 numpy / pandas / scipy / scikit-learn / matplotlib（本机全有）；IGNITE 路径需上游仓库 clone + numba（本机 0.65.1 已装，无需再装）。上游无 requirements 文件、无 LICENSE。 | Python | assets/585_fig1_networks.png（真值 vs 各法耦合矩阵 heatmap）、assets/585_fig2_edge_recovery.png（边恢复 ROC 曲线 + AUROC/AUPRC dumbbell）、assets/585_fig3_weight_agreement.png（真值-推断权重散点 + 真边/非边权重 violin）；各含配套 PDF。全部非条形图。 |
| 🟡 | 586 | [586_psgrn_grn_inference](03_virtual_perturbation/02_grn_inference/586_psgrn_grn_inference) | 586 · PSGRN — 从带干预标签的单细胞扰动矩阵(CRISPRi/Perturb-seq 风格,含 non-targeting 对照)推断有向基因调控网络。忠实复现上游算法:相关性造合成金标准 → 4 个扰动特征 → LightGBM 自训练重排全部有序基因对;内置两条朴素基线(共表达 |Pearson|、单变量干预效应)做强制对照,出 PR 曲线 / precision@K dumbbell / 打分热图。上游官方 CausalBench 评测入口做守卫式封装,缺包时打印真实命令、不伪造返回值。 | 输入:example_data/perturb_expression.csv(行=细胞,首列 intervention,'non-targeting'=对照,其余列=基因表达)+ 可选 ground_truth_edges.csv(From,To)。输出:results/edges_psgrn.csv、edges_perturbation.csv、edges_co-expression.csv、586_summary.json(参数/后端/AUPRC/P@K/依赖版本快照);assets/ 三张图(PDF+300dpi PNG)。 | Python · numpy 2.4.6 · pandas 2.3.3 · scikit-learn 1.8.0 · lightgbm 4.6.0 · matplotlib 3.11.0(全部本机已装,零安装可跑)。可选上游官方评测路径需另装 causalscbench(本机未装,已守卫式跳过)。 | Python | 3 张,全部无条形图:fig1_pr_curves(PR 曲线,三方法 + 随机基准线)、fig2_precision_at_k(dumbbell + dot,PSGRN vs 最佳基线在 top-10/25/50/100 的位移)、fig3_score_matrix(heatmap + 真值边红圈散点覆盖)。 |
| 🟡 | 587 | [587_regformer_grn_mamba_fm](03_virtual_perturbation/02_grn_inference/587_regformer_grn_mamba_fm) | RegFormer GRN 重建评测台：把上游「基因嵌入→余弦相似度 TF 有向图→谱聚类模块」下游链路本地复刻，用共表达/PCA 朴素嵌入作必跑基线，可插入任意外部 gene_embedding.npy 在同一口径下对比。 | 输入 example_data/expression_counts.csv（400 cells × 132 genes 合成计数）+ tf_list.txt（无表头 TSV 第1列，12 TF）+ 可选 ground_truth_edges.csv / ground_truth_modules.csv + 可选 --embedding/--gene-names 外部嵌入；输出 results/ 下 edges_<emb>.csv、modules_<emb>.csv、per_tf_recall.csv、metrics_table.csv、587_summary.json（含 session_info 版本快照），assets/ 下 4 张 PNG+PDF | numpy 2.4.6 / pandas 2.3.3 / scikit-learn 1.8.0 / networkx 3.6.1 / matplotlib 3.11.0（本机已装，基线路径零安装）；RegFormer 本体需 clone 仓库 + CUDA（torch==2.0.0+cuda11.6, dgl==1.1.2+cu118, mamba_ssm==1.1.1, causal_conv1d==1.1.1），不在 PyPI | Python 3.12 | fig1_pr_curve（PR 曲线+随机基线线）、fig2_metric_slopegraph（多指标 slopegraph）、fig3_similarity_heatmap（基因-基因余弦相似度热图，按真实调控程序排序）、fig4_tf_recall_lollipop（每 TF 靶基因召回 lollipop）；全部非条形图，grep 无 .bar/.barh/geom_bar |

### 03.03 · 因果表示与反事实 — Causal & counterfactual

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 588 | [588_sccausalvi_causal_perturbation](03_virtual_perturbation/03_causal_perturbation/588_sccausalvi_causal_perturbation) | 案例-对照单细胞扰动响应的因果解耦模块：默认跑可复现线性基线（条件中心化 PCA 背景表示 + 全局/细胞类型特异 Δ 反事实 + kNN 响应细胞打分），scCausalVI 深度模型为守卫式可选路径。 | 输入 example_data/synthetic_counts.npz（X 计数矩阵 960×300 + condition/cell_type/batch/gene_names/true_responder）；输出 results/disentanglement_metrics.csv、per_cell_perturbation_score.csv、celltype_treatment_deltas.csv、588_summary.json，assets/ 下 5 组 PNG+PDF；启用 --run-sccausalvi 且已装包时另出 sccausalvi_latent_bg.npy / sccausalvi_latent_te.npy / sccausalvi_responsive_cells.csv | 基线（本机已装，零安装即跑）：numpy 2.4.6 / pandas 2.3.3 / scikit-learn 1.8.0 / matplotlib 3.11.0，Python 3.12.5；可选完整方法：scCausalVI 0.0.11（上游 setup.py 声明 python>=3.9、numpy>=1.23.5,<2.0、scvi-tools>=0.16.1、torch>=2.0.0、anndata>=0.10.3、scanpy>=1.9.6、pytorch-lightning>=1.5.10、gdown>=5.2.0；MIT License），本机未安装、与本机 numpy 2.x 冲突需独立环境 | Python | fig1_background_embedding（2×2 散点：raw PCA vs 条件中心化 PCA，按 condition / cell type 上色）、fig2_disentanglement_slopegraph（slopegraph）、fig3_response_score_violin（小提琴+抖动点）、fig4_counterfactual_scatter（留出集 预测Δ vs 观测Δ 散点）、fig5_top_response_genes_lollipop（lollipop）；各出 PNG+PDF，无条形图 |

### 03.04 · 药物扰动与响应 — Drug perturbation

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🔴 | 070 | [070_chemcpa_drug_perturbation.py](03_virtual_perturbation/04_drug_perturbation/070_chemcpa_drug_perturbation.py) | chemCPA drug-perturbation expression prediction (GPU) | repo + config → train logs | Py · chemCPA, torch | Python | — |
| 🔴 | 071 | [071_scdrug_response_prediction.py](03_virtual_perturbation/04_drug_perturbation/071_scdrug_response_prediction.py) | scDrug single-cell drug response (heavy) | 10x/h5ad → cluster drug response | Py · scDrug, GDSC/PRISM | Python | — |
| 🟡 | 518 | [518_beyondcell_drug_response](03_virtual_perturbation/04_drug_perturbation/518_beyondcell_drug_response) | beyondcell core re-impl: BCS + therapeutic clusters | scRNA + drug signatures → BCS/ranking | R · UCell, ggplot2 | R | heatmap, lollipop, UMAP |
| 🟡 | 589 | [589_scdruglink_drug_response](03_virtual_perturbation/04_drug_perturbation/589_scdruglink_drug_response) | 按 scDrugLink 上游源码复刻的单细胞药物重定位打分:Drug2Cell 靶点臂(促进/抑制)× 扰动签名臂(敏感/耐药)在细胞类型层面 exp(weight) 相乘串联,输出全图谱与细胞类型两级治疗评分排序 + AUROC/AUPR 对照。 | 输入 example_data/:expr_lognorm.csv(基因×细胞,log-normalised)、cell_meta.csv(cell/cell_type/disease)、drug_targets.csv(drug_name/gene_names,";"分隔)、drug_perturb_signature.csv(基因×药物 z 值)、drug_labels.csv(drug_name/known_label,可选)。输出 results/ 8 个 csv(drug_target_d2c、prom_inh_weight、prom_inh_padj、sens_res_fdr、connectivity_score、drug_scores、cell_type_drug_scores、eval_metrics)+ assets/ 4 图(PNG+PDF)。CLI:--expr/--meta/--targets/--sig/--labels/--disease/--control/--outdir/--n_bins/--ctrl_size/--n_perm。 | 必需 ggplot2 + 共享 _framework/theme_pub.R;可选 Seurat(本机 5.5.0,有则走 FindMarkers 做 DEG,无则用内置 wilcox+Bonferroni 等价实现)、ggrepel(散点标签)。不需要上游 scDrugLink / Asgard / cmapR / effsize(算法按源码在 base R 内复刻)。完整 B 臂需服务器装 Asgard + cmapR + CMAP L1000 GSE70138/GSE92742 gctx。 | R | prom_inh_heatmap(A 臂权重 RdBu 发散热图,细胞类型×药物)、drug_score_lollipop(Top15 治疗评分 lollipop,log1p 轴)、linking_effect_scatter(unlinked vs linked 散点 + 对角参考线)、d2c_violin_top_drug(头名药 Drug2Cell 分的 violin+boxplot,细胞类型×病/对照)。全部非条形图。 |

### 03.05 · 扰动预测基准 — Perturbation benchmarks

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 590 | [590_scperturbench_generalization](03_virtual_perturbation/05_benchmark/590_scperturbench_generalization) | 把自己的单细胞扰动预测按 scPerturBench (Nat Methods 2026) 的同一套指标打分，并强制与上游口径的朴素基线 (controlBaseline / trainMean) 对照，判断深度模型是否真的赢过"什么都不预测"。 | 输入 example_data/observed.csv (cell_id, group∈{control,stimulated}, 基因列) + predicted.csv (cell_id, method, 基因列，可含多方法)；输出 results/performance.tsv (长表) + performance_wide.csv + rank_matrix.csv + summary.json (含 mean_rank / verdict_vs_naive_floor / caveats / 依赖版本快照) + 3 张 PNG+PDF 图 | numpy, pandas, scipy, matplotlib（全部本机已装，零额外安装）；可选 pertpy + anndata 走 --use-pertpy 守卫路径拿上游同名 Sinkhorn wasserstein | Python 3.12 | 3 张，均非条形图：fig1_rank_heatmap（方法×指标排名热图）、fig2_vs_baseline_dotplot（相对 controlBaseline 比值的点图，log 轴 + 基线虚线）、fig3_delta_scatter（每基因 delta 预测 vs 真实散点，每方法一 panel，标注 PCC-delta） |
| 🟡 | 591 | [591_scarchon_perturbation_benchmark](03_virtual_perturbation/05_benchmark/591_scarchon_perturbation_benchmark) | 留一批次(leave-one-batch-out)的单细胞扰动响应预测基准骨架:按 scArchon 口径评估预测的扰动后表达,强制与 control/mean 朴素基线同台对比 | In: .h5ad (AnnData) with obs columns `condition` (control/stimulated) and `batch`; defaults to bundled example_data/synthetic_perturb.h5ad (1800 cells x 300 genes, 4 donors). Out: results/591_scores_per_batch.csv, 591_scores_mean.csv, 591_per_gene_delta.csv, 591_summary.json, 4 PDF+PNG figures; display copies in assets/. | anndata, numpy, scipy, pandas, scikit-learn, matplotlib (all locally installed); no scArchon import — upstream is Snakemake-only. Optional server path: snakemake + singularity/apptainer. | Python | 4 figures, no bar charts: dumbbell (delta-R2 per held-out batch vs control floor), heatmap (method x metric, raw value in cell + direction-normalised colour), scatter (per-gene true vs predicted delta, 3 panels), raincloud (violin + jitter + box of per-gene absolute error). |

## 04 · 因果推断与遗传流行病 — Causal inference & genetics  (29)

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
| 🟡 | 594 | [594_colocboost_colocalization](04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | 同一基因座多性状(GWAS+eQTL/sQTL/pQTL)共定位:真包 coloc::coloc.abf 两两基线 + 守卫式 colocboost 多性状联合封装,出 dot/heatmap/dumbbell/lollipop 四图 | 输入 example_data/sumstat_<trait>.csv ×4 (variant,pos,maf,beta,se,z,n) + region_ld.csv(带 dimnames 的 LD 方阵) + true_causal.csv;输出 results/baseline_pairwise_coloc.csv、results/versions.txt(装了真包另出 colocboost_cos_summary.csv / colocboost_vcp.csv)、assets/ 4×PNG + 4×PDF | R ≥4.0;本机已装 coloc 5.2.3 + ggplot2 4.0.3(基线路径,始终可跑);可选 colocboost 1.0.9(CRAN,Imports: Rfast, matrixStats)—— 本机未装,守卫跳过 | R | region_multitrait_dots(散点/facet)、pairwise_pph4_heatmap(热图 viridis)、pph_dumbbell(哑铃 H3 vs H4)、variant_level_lollipop(棒棒糖 top15);无条形图 |

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
| 🟡 | 592 | [592_twist_transcriptome_wide_test](04_causal_inference_genetics/06_twas_sceqtl/592_twist_transcriptome_wide_test) | 拟时序(细胞状态)分辨的 TWAS:用 B-spline eQTL 权重矩阵沿拟时序逐点做 FUSION 式 burden 检验,与静态 TWAS 同框对照;正式 TWiST 三联检验为守卫式封装。 | 输入 example_data/(合成): gwas_sumstats_synth.txt(SNP/A1/A2/Z)、twist_weights_Wmat_synth.csv(ID/SNP/basis1..basis7)、wgtlist_synth.csv(ID/CHR/P0/P1/tss)、ld_reference_synth.csv、ngwas_synth.txt;输出 results/ 三个 CSV(每基因 z_static/p_static/t_peak/z_peak/p_scan_bonf/FDR、基因×拟时序 z 矩阵、按情景检出率)+ assets/ 四张 PNG/PDF | R 4.4+;基线仅 splines(base) + ggplot2 + modules/_framework/theme_pub.R;正式路径需 TWiST + plink2R + fda + grpreg(不自动安装) | R | fig1_static_vs_scan_dumbbell(dumbbell)、fig2_pseudotime_z_heatmap(发散色 heatmap)、fig3_effect_trajectories(轨迹曲线)、fig4_detection_rate_slopegraph(slopegraph);无条形图 |
| 🟡 | 593 | [593_case_celltype_eqtl_finemap](04_causal_inference_genetics/06_twas_sceqtl/593_case_celltype_eqtl_finemap) | 多细胞类型 eQTL 联合精细定位:区分跨细胞类型共享效应与细胞类型特异效应,内置「完全特异」「完全共享」两条纯 base R 极端基线 + 守卫式 CASE 上游调用 | 输入 example_data/genotypes.csv(样本×SNP 剂量 0/1/2)+ expression.csv(样本×细胞类型伪bulk 表达)+ 可选 true_eqtl.csv;输出 results/593_pip_table.csv、593_credible_sets.csv、593_truth_pip.csv、593_zscores.csv、593_sessionInfo.txt + assets/ 4 图(PNG+PDF) | R 4.4 base + ggplot2(经 _framework/theme_pub.R);上游 CASE 0.3.1 可选(其自身 Imports: magrittr/MASS/mvtnorm/stats,Depends R>=4.0.0,GPL(>=3)) | R | fig1_zscore_scatter(分面散点+Bonferroni 线)、fig2_pip_heatmap(热图)、fig3_sharing_slopegraph(slopegraph)、fig4_credible_set_size(棒棒糖点图);无条形图 |

### 04.07 · 稳健 MR 估计量 — Robust MR estimators

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 533 | [533_mrcare_winnerscurse_mr](04_causal_inference_genetics/07_robust_mr_methods/533_mrcare_winnerscurse_mr) | Winner's-curse-corrected MR (CARE/RIVW) vs naive baseline | two-sample MR summary → estimates + plots | R · TwoSampleMR (baseline); MRcare | R | lollipop, forest, scatter, dumbbell |
| 🟡 | 595 | [595_mreills_robust_mr](04_causal_inference_genetics/07_robust_mr_methods/595_mreills_robust_mr) | 不变性(EILLS)稳健 MR：整合多个异质 GWAS summary 数据集，对含水平多效性的无效工具做筛选并给出单/多暴露因果估计，与 MVMR-IVW / MR-Egger 同数据对照 | 输入 example_data/mreills_multi_gwas_summary.csv（synthetic 长表：dataset, SNP, beta_X*/se_X*, beta_Y/se_Y, true_invalid_IV；3 数据集 × 150 SNP × 2 暴露）；输出 results/MR_estimates_eills_vs_baseline.csv、IV_selection_QSj.csv、lambda_path.csv、sessionInfo.txt + assets/ 4 图（PDF+PNG） | R 4.4.3；ggplot2 + 模块框架 _framework/theme_pub.R（bio_args/pal_pub/theme_pub/save_fig/scale_fill_diverge 均已在框架 line 81/94/140/153/204 确认存在）；基线与本地转写实现仅用 base R，无需装包；可选官方包 MREILLS（devtools::install_github("hhoulei/MREILLS")，本机未装） | R | fig1_estimates_dotwhisker（dot-and-whisker，四方法估计+95%CI，按暴露分面）；fig2_iv_screening_violin（violin+jitter，QSj 按真实工具有效性分层+lambda 阈值线）；fig3_lambda_sensitivity（折线+点，点大小=入选工具数）；fig4_bias_heatmap（RdBu 发散热图，有符号偏倚）。无条形图 |

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

## 08 · 结构与药物设计 — Structure & drug design  (6)

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

### 08.03 · 虚拟筛选与打分 — Virtual screening & scoring

| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |
|----|---|------|------|------------|------|------|------|
| 🟡 | 596 | [596_scorch2_virtual_screening](08_structure_drug_design/03_virtual_screening/596_scorch2_virtual_screening) | SCORCH2 双视图共识 ML 重打分的守卫式封装 + 本机可跑的虚拟筛选富集评测骨架(EF1%/EF5%/BEDROC/AUROC,按靶点分层、GroupKFold 防泄漏) | 输入 pose 级特征表 CSV(target/compound_id/docking_score/label + 8 个 PS 相互作用特征 + 6 个 PB 理化位姿特征;缺省自动生成合成样例);输出 results/ 下 overall_metrics.csv、per_target_metrics.csv、consensus_weight_sweep.csv、compound_level_scores.csv、596_summary.json 及 4 张图(PNG+PDF),展示图同步 assets/ | 本机基线零安装:numpy / pandas / scikit-learn / matplotlib(均已装)。上游 SCORCH2 本体需 conda env(Python>=3.10, xgboost, rdkit, openbabel)+ Zenodo 权重 + ADFRsuite(仅 pdbqt 转换用) | Python | fig1 富集曲线(折线) · fig2 每靶点 EF1% slopegraph · fig3 active/decoy 分数分离 violin+抖动散点 · fig4 共识权重扫描点线 + 两视图散点。无条形图 |

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
- **Heatmap · 热图** → 010, 018, 021, 022, 026, 045, 046, 049, 054, 057, 059, 076, 077, 078, 081, 083, 084, 492, 496, 502, 506, 510, 511, 516, 518, 520, 523, 524, 525, 530, 531, 534, 536, 537, 539, 543, 544, 545, 547, 548, 549, 550, 552, 554, 556, 559, 560, 561, 562, 565, 566, 568, 569, 571, 572, 573, 575, 576, 577, 578, 579, 584, 585, 586, 587, 589, 590, 591, 592, 593, 594, 595
- **ROC / PR** → 016, 034, 045, 052, 057, 059, 063, 550, 551, 553, 561, 567, 568, 585
- **Calibration / DCA / nomogram · 校准与决策曲线** → 016, 063, 550, 553, 555
- **Forest · 森林图** → 016, 032, 033, 043, 078, 503, 508, 516, 519, 527, 528, 529, 533, 534, 535, 536
- **KM / survival · 生存曲线** → 048, 057, 497, 529, 551
- **UMAP / tSNE · 降维嵌入** → 026, 027, 044, 046, 049, 050, 058, 062, 081, 506, 507, 516, 518, 532, 541, 560, 563, 564, 569, 570, 572, 574, 588
- **Violin · 小提琴** → 026, 027, 044, 046, 049, 050, 509, 529, 532, 536, 542, 544, 545, 555, 558, 563, 568, 570, 573, 581, 585, 588, 589, 591, 595, 596
- **Raincloud · 云雨图** → 512, 554, 557, 559, 563, 581, 591
- **Ridgeline · 山脊图** → 513
- **Box · 箱线图** → 016, 021, 056, 492, 503, 557, 581, 589, 591
- **Dot / bubble · 点图气泡图** → 007, 026, 046, 049, 050, 051, 062, 082, 510, 531, 532, 546, 550, 557, 561, 583, 586, 590, 592, 594, 595
- **Lollipop · 棒棒糖图(条形图替代)** → 003, 005, 006, 007, 011, 013, 014, 022, 034, 035, 045, 058, 496, 502, 503, 507, 509, 511, 518, 521, 527, 528, 530, 531, 533, 535, 537, 538, 539, 540, 541, 542, 543, 547, 548, 549, 550, 551, 552, 553, 554, 556, 557, 559, 560, 562, 565, 568, 570, 572, 574, 576, 577, 581, 582, 583, 584, 587, 588, 589, 594
- **Dumbbell / slopegraph · 哑铃图与斜率图** → 514, 528, 533, 534, 537, 539, 540, 544, 548, 552, 555, 556, 559, 561, 562, 564, 569, 571, 575, 576, 577, 578, 579, 582, 583, 585, 586, 587, 588, 591, 592, 593, 594, 596
- **Venn / UpSet · 集合图** → 003, 005, 006, 011, 015, 034, 035, 509, 554
- **PCA** → 010, 026, 027, 056, 548, 568, 572, 588
- **Scatter · 散点** → 012, 013, 014, 032, 033, 043, 086, 495, 506, 511, 519, 520, 521, 527, 533, 535, 536, 537, 538, 541, 542, 543, 544, 545, 548, 560, 562, 565, 570, 571, 578, 588, 589, 590, 591, 593
- **Network · 网络图** → 007, 047, 514, 530, 531, 539, 540, 546, 558, 562, 569, 579, 583, 585, 587, 588, 592, 593, 596
- **Chord / circos / alluvial · 弦图环图桑基** → 047, 051, 053, 498, 515, 532, 549
- **Trajectory / vector field · 轨迹与向量场** → 044, 049, 050, 062, 082, 517, 548, 561, 581, 592
- **Spatial map · 空间分布图** → 027, 050, 073, 080, 505, 521, 541, 542, 543, 544, 545, 569, 574, 575, 576, 578, 579, 581
- **Feature map · 基因表达投影** → 026, 027, 044, 046, 050, 532
- **Composite multi-panel · 多面板拼图** → 516
