# Module catalog

One line per module: purpose, input → output, key dependencies, language, and the
figure types it produces. For the high-level category overview see the
[top-level README](../README.md); for the new-project scaffold and shared style see
[`_framework/`](_framework/).

> **Reuse first, never from scratch.** Pick a module here (or a real tool) and adapt
> it — do not hand-write analysis code from memory (that invites hallucinated APIs and
> wrong parameters). See [`_framework/CONVENTIONS.md` §0](_framework/CONVENTIONS.md).

## Status legend

| Mark | Meaning |
|------|---------|
| ✅ | Turnkey — runs locally on bundled/synthetic example data, no edits |
| 🟡 | Core re-implementation or honest baseline runs locally; the full method needs the package on the analysis server (see [`_framework/SERVER_DEPENDENCIES.md`](_framework/SERVER_DEPENDENCIES.md)) |
| 🔴 | Heavy / GPU / external toolchain (FUSION, GROMACS, deep-learning FMs) — reference wrapper, not locally rendered |
| 📄 | Template or upstream-derived script — bring your own data + install; no bundled example |
| 📦 | Vendored third-party package — kept manifest-only / local |
| 🗃️ | Local-only, git-ignored — not in the public repo (kept for the author's drafts) |

---

## 01 · Network pharmacology & target databases

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 001 | [001_ctd_compound_targets](01_network_pharmacology/001_ctd_compound_targets) | Extract & de-dup compound targets from a CTD export | CTD csv → targets.csv | R base | R | — |
| ✅ | 002 | [002_swisstarget_compound_targets](01_network_pharmacology/002_swisstarget_compound_targets) | Extract compound targets from a SwissTargetPrediction export | Swiss csv → targets.csv | R base | R | — |
| ✅ | 003 | [003_ctd_swiss_target_union_venn](01_network_pharmacology/003_ctd_swiss_target_union_venn) | Union/intersection of CTD vs Swiss compound targets | target csvs → set table + plot | R · UpSetR | R | venn, upset, lollipop |
| ✅ | 004 | [004_genecards_disease_targets](01_network_pharmacology/004_genecards_disease_targets) | Extract disease targets from a GeneCards export | GeneCards csv → targets.csv | R base | R | — |
| ✅ | 005 | [005_omim_genecards_target_venn](01_network_pharmacology/005_omim_genecards_target_venn) | Union/intersection of OMIM vs GeneCards disease targets | target csvs → set table + plot | R · UpSetR | R | venn, upset, lollipop |
| ✅ | 006 | [006_disease_compound_target_venn](01_network_pharmacology/006_disease_compound_target_venn) | Disease ∩ compound targets → core targets | target csvs → intersection + plot | R · UpSetR | R | venn, lollipop |
| ✅ | 011 | [011_deg_drug_target_intersection](01_network_pharmacology/011_deg_drug_target_intersection) | Multi-set DEG ∩ drug ∩ disease target intersection | gene/target csvs → intersection + plot | R · UpSetR | R | venn, upset, lollipop |
| 📄 | 493 | [493_opentargets_dgidb_chembl_druggability.py](01_network_pharmacology/493_opentargets_dgidb_chembl_druggability.py) | Composite druggability score for a target set (live APIs) | HGNC genes → score DataFrame | Py · requests, pandas, mygene | Python | — |

## 02 · GO / KEGG enrichment

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 007 | [007_go_kegg_enrichment](02_enrichment/007_go_kegg_enrichment) | GO/KEGG over-representation for a gene list | gene_list.csv → enrichment table + plots | R · clusterProfiler, org.Hs.eg.db, ggraph | R | dot, lollipop, network |

## 03 · Transcriptomics (GEO) & differential expression

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 008 | [008_geo_expression_matrix](03_transcriptomics_deg/008_geo_expression_matrix) | Annotate GEO probes → gene-level matrix | series_matrix + GPL → geneMatrix.csv | R base | R | — |
| ✅ | 009 | [009_geo_sample_grouping](03_transcriptomics_deg/009_geo_sample_grouping) | Normalize matrix + attach group labels | geneMatrix + groups → labelled matrix | R · limma | R | — |
| ✅ | 010 | [010_geo_deg_volcano_heatmap_pca](03_transcriptomics_deg/010_geo_deg_volcano_heatmap_pca) | limma two-group DE with three figures | expr matrix → DEG table + plots | R · limma, ComplexHeatmap | R | volcano, heatmap, PCA |
| ✅ | 056 | [056_geo_multicohort_batch_correction](03_transcriptomics_deg/056_geo_multicohort_batch_correction) | Merge multi-cohort data + remove batch effect | cohort dir → corrected matrix + QC | R · limma/sva | R | PCA, box |

## 04 · Machine-learning feature selection

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 012 | [012_lasso_feature_selection](04_ml_feature_selection/012_lasso_feature_selection) | LASSO logistic feature selection | expr + candidates → selected genes + plot | R · glmnet | R | scatter (CV / coef path) |
| ✅ | 013 | [013_svm_rfe_feature_selection](04_ml_feature_selection/013_svm_rfe_feature_selection) | SVM-RFE recursive elimination ranking | expr + candidates → ranking + subset | R · e1071 | R | scatter, lollipop |
| ✅ | 014 | [014_randomforest_feature_selection](04_ml_feature_selection/014_randomforest_feature_selection) | Random-forest Gini importance ranking | expr + candidates → importance + plot | R · randomForest | R | lollipop, scatter (OOB) |
| ✅ | 015 | [015_ml_feature_intersection_venn_upset](04_ml_feature_selection/015_ml_feature_intersection_venn_upset) | Intersect feature genes across methods | gene-list dir → intersection + plot | R · UpSetR | R | venn, upset |
| ✅ | 034 | [034_multi_ml_feature_selection](04_ml_feature_selection/034_multi_ml_feature_selection) | caret 10-method models → AUC + consensus features | expr → AUC table + consensus | R · caret, pROC, UpSetR | R | ROC, lollipop, upset |
| ✅ | 035 | [035_ml_combination_intersection](04_ml_feature_selection/035_ml_combination_intersection) | Rank method combinations by intersection size | method feature lists → combo table | R · UpSetR | R | lollipop, upset |
| 📄 | 045 | [045_multimodalad_ml_models.R](04_ml_feature_selection/045_multimodalad_ml_models.R) | Multi-ML integrated modelling (RSF/LASSO/GBM/BART…) | train/test txt → models + plots | R · randomForestSRC, glmnet, gbm, BART, xgboost | R | ROC, heatmap, lollipop |
| ✅ | 052 | [052_shap_interpretation](04_ml_feature_selection/052_shap_interpretation) | Train best model + SHAP interpretation | geneexp.csv → SHAP table + plots | R · caret, kernelshap, shapviz | R | beeswarm, dependence, waterfall, force, ROC |
| 📄 | 059 | [059_dual_disease_15ml_175combos.R](04_ml_feature_selection/059_dual_disease_15ml_175combos.R) | Two-disease 15-ML × 175-combo screen + modelling | train/test → models + AUC heatmap | R · randomForestSRC, glmnet, ComplexHeatmap, sva | R | heatmap (AUC), ROC |
| 📄 | 496 | [496_mime_101combo_prognostic.R](04_ml_feature_selection/496_mime_101combo_prognostic.R) | Mime 10-algorithm × 101-combo prognostic signature (usage snippet) | survival lists + genes → C-index table | R · Mime1 | R | lollipop, heatmap (C-index) |
| ✅ | 502 | [502_biomarker_triple_vote](04_ml_feature_selection/502_biomarker_triple_vote) | Topology × correlation × Boruta triple-vote shortlist | expr + group + candidates → vote/consensus | R · igraph, Boruta, Hmisc | R | heatmap (vote), lollipop |

## 05 · Diagnostic models & validation

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 016 | [016_diagnostic_model_roc_calibration_dca](05_diagnostic_models/016_diagnostic_model_roc_calibration_dca) | Logistic diagnostic model, full clinical evaluation | expr + genes → evaluation figures | R · rms, rmda, pROC | R | nomogram, calibration, DCA, ROC, forest, box |
| ✅ | 063 | [063_geo_diagnostic_validation](05_diagnostic_models/063_geo_diagnostic_validation) | External-cohort validation of a diagnostic model | train + valid matrices → AUC + plot | R · rms, pROC | R | ROC, calibration |
| ✅ | 503 | [503_generalization_robustness](05_diagnostic_models/503_generalization_robustness) | Meta-analysis + LODO cross-cohort generalization | cohorts.rds → LODO/weight tables | R · metafor, glmnet, pROC | R | forest, lollipop, box |

## 06 · Immune infiltration

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 📄 | 017 | [017_immune_infiltration_source.R](06_immune_infiltration/017_immune_infiltration_source.R) | CIBERSORT deconvolution engine (SVR source) | expr + LM22 → cell fractions | R · e1071 | R | — |
| 📄 | 018 | [018_immune_infiltration_scoring.R](06_immune_infiltration/018_immune_infiltration_scoring.R) | Immune infiltration scoring + cell/function matrices | expr → score matrix + corr | R · e1071, preprocessCore | R | heatmap (corr) |
| ✅ | 021 | [021_immune_infiltration_viz/](06_immune_infiltration/021_immune_infiltration_viz) | Fraction-matrix difference/composition/correlation viz | CIBERSORT csv → 3 figures | R · ggpubr, ComplexHeatmap | R | box, stacked-bar, heatmap |
| 🟡 | 492 | [492_iobr_multimethod_deconvolution.R](06_immune_infiltration/492_iobr_multimethod_deconvolution.R) | IOBR 7-method deconvolution + cross-method consistency | expr + group → merged matrix + plots | R · IOBR, tidyverse | R | stacked-bar, box, heatmap |
| ✅ | 520 | [520_bayesprism_deconvolution](06_immune_infiltration/520_bayesprism_deconvolution) | BayesPrism Bayesian deconvolution with ground-truth check | scRNA ref + bulk → fractions + accuracy | R · BayesPrism | R | scatter, heatmap |

## 07 · Molecular docking & dynamics

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 022 | [022_docking_binding_energy_viz](07_molecular_docking/022_docking_binding_energy_viz) | Binding-energy heatmap + strongest-binding ranking | binding_energy.csv → figures | R · ComplexHeatmap | R | heatmap, lollipop |
| 🔴 | 086 | [086_vina_gromacs_mmpbsa_mdanalysis_pipeline.py](07_molecular_docking/086_vina_gromacs_mmpbsa_mdanalysis_pipeline.py) | Vina docking + GROMACS MD + MM-PBSA pipeline | receptor/ligand → trajectory + ΔG | Py · Vina, GROMACS, gmx_MMPBSA, MDAnalysis | Python | scatter (RMSD/RMSF/Rg/SASA/energy) |

## 08 · Single-cell / spatial / trajectory

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 📄 | 023 | [023_rds_structure_check.R](08_singlecell_spatial_trajectory/023_rds_structure_check.R) | Inspect an RDS object's structure | RDS → structure print + csv | R base | R | — |
| 📄 | 024 | [024_scrna_rds_prep.R](08_singlecell_spatial_trajectory/024_scrna_rds_prep.R) | Build a Seurat/RDS object from raw data | raw → Seurat object | R · Seurat, SingleR, celldex | R | — |
| 📄 | 025 | [025_scrna_data_prep.R](08_singlecell_spatial_trajectory/025_scrna_data_prep.R) | Read 10x data into a Seurat object | 10x raw → Seurat object | R · Seurat, Matrix | R | — |
| 📄 | 026 | [026_scrna_seurat_sctenifoldknk_ko.R](08_singlecell_spatial_trajectory/026_scrna_seurat_sctenifoldknk_ko.R) | Full QC/cluster/annotate + virtual KO pipeline | counts → object + full figure set | R · Seurat, monocle3, scTenifoldKnk, CellChat | R | UMAP, heatmap, PCA, violin, dot, feature-map |
| 📄 | 027 | [027_spatial_seurat_auto.R](08_singlecell_spatial_trajectory/027_spatial_seurat_auto.R) | Spatial transcriptomics read/cluster/visualize | Visium h5 → spatial figures | R · Seurat, SingleR, glmGamPoi | R | PCA, violin, UMAP, feature-map, niche-map |
| 📄 | 044 | [044_multimodalad_scrna.R](08_singlecell_spatial_trajectory/044_multimodalad_scrna.R) | AD brain scRNA pipeline + Monocle pseudotime | GSE157827 → object + trajectory | R · Seurat, SingleR, monocle | R | violin, UMAP, trajectory, feature-map |
| ✅ | 046 | [046_scrna_publication_figures](08_singlecell_spatial_trajectory/046_scrna_publication_figures) | Standard Seurat flow → publication figures | counts.csv → object + figures | R · Seurat, ggplot2 | R | UMAP, dotplot, heatmap, feature-map, violin |
| 📄 | 049 | [049_scrna_manual_annot_cellchat_trajectory.R](08_singlecell_spatial_trajectory/049_scrna_manual_annot_cellchat_trajectory.R) | Manual annotation + CellChat + trajectory | Seurat + markers → annotation + figures | R · Seurat, CellChat, monocle3 | R | violin, UMAP, heatmap, dotplot, trajectory |
| 📄 | 050 | [050_spatial_cluster_annot_trajectory.R](08_singlecell_spatial_trajectory/050_spatial_cluster_annot_trajectory.R) | Spatial cluster annotation + monocle3 pseudotime | Visium → spatial annotation + trajectory | R · Seurat, monocle3, patchwork | R | niche-map, violin, UMAP, feature-map, pseudotime |
| 📄 | 051 | [051_scrna_cellchat.R](08_singlecell_spatial_trajectory/051_scrna_cellchat.R) | CellChat cell-communication network | annotated expr + labels → comm. figures | R · Seurat, CellChat | R | circos, bubble |
| 📄 | 058 | [058_scrna_scissor.R](08_singlecell_spatial_trajectory/058_scrna_scissor.R) | Scissor — link bulk phenotype to disease-relevant cells | bulk pheno + Seurat → Scissor cells | R · Seurat, Scissor | R | UMAP, lollipop |
| 📄 | 061 | [061_scfocal_gui_input_prep.R](08_singlecell_spatial_trajectory/061_scfocal_gui_input_prep.R) | Prepare scFOCAL GUI input + launch Shiny | RData + map csv → RDS + GUI | R · Seurat, scFOCAL, shiny | R | — (interactive) |
| 🟡 | 062 | [062_sctour_pseudotime_vectorfield.py](08_singlecell_spatial_trajectory/062_sctour_pseudotime_vectorfield.py) | scTour pseudotime + latent space + vector field | h5ad → pseudotime + vector field | Py · scanpy, sctour | Python | UMAP, pseudotime, vector-field |
| 📄 | 082 | [082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R](08_singlecell_spatial_trajectory/082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R) | Slingshot / tradeSeq / CytoTRACE2 trajectory consensus | Seurat RDS → pseudotime table | R · slingshot, tradeSeq, CytoTRACE2 | R | pseudotime |
| 📄 | 087 | [087_palantir_branch_probability.py](08_singlecell_spatial_trajectory/087_palantir_branch_probability.py) | Palantir pseudotime + branch probability + entropy | h5ad + root → pseudotime/branch csv | Py · palantir, scanpy | Python | — (tables) |
| 🟡 | 506 | [506_scvi_scanvi_integration](08_singlecell_spatial_trajectory/506_scvi_scanvi_integration) | scVI/scANVI integration + label transfer (vs PCA baseline) | h5ad (batch/label) → integration + labels | Py · scvi-tools, scanpy, sklearn | Python | UMAP, scatter, heatmap (confusion) |
| ✅ | 517 | [517_vector_trajectory_direction](08_singlecell_spatial_trajectory/517_vector_trajectory_direction) | VECTOR expression-potential differentiation direction | embedding + expr → potential + field | R · ggplot2 | R | vector-field |

> `491_sctour_extra_files/` is a support folder (tutorial run scripts + env notes) for module 062, not a standalone module.

## 09 · Mendelian randomization & GWAS

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 📄 | 028 | [028_mr_gwas_exposure_vcf_sig_snp.R](09_mendelian_randomization/028_mr_gwas_exposure_vcf_sig_snp.R) | Pick genome-wide significant SNPs from a VCF exposure | exposure vcf → significant SNPs | R · VariantAnnotation, gwasglue | R | — |
| 📄 | 029 | [029_mr_gwas_exposure_ld_clump.R](09_mendelian_randomization/029_mr_gwas_exposure_ld_clump.R) | LD-clump to independent instruments | candidate SNPs → independent IVs | R · gwasglue, TwoSampleMR | R | — |
| 📄 | 030 | [030_mr_gwas_exposure_add_eaf.R](09_mendelian_randomization/030_mr_gwas_exposure_add_eaf.R) | Add effect-allele frequency to exposure SNPs | SNPs → SNPs + EAF | R · ieugwasr | R | — |
| 📄 | 031 | [031_mr_gwas_exposure_weak_iv_filter.R](09_mendelian_randomization/031_mr_gwas_exposure_weak_iv_filter.R) | F-statistic weak-instrument filter | SNPs → strong IVs (F≥10) | R · ieugwasr | R | — |
| ✅ | 032 | [032_mr_twosamplemr](09_mendelian_randomization/032_mr_twosamplemr) | Self-contained MR causal inference (primary) | harmonized.csv → estimates + plots | R · ggplot2 (self-contained MR) | R | scatter, forest, funnel, leave-one-out |
| 📄 | 033 | [033_mr_gwas_finngen_outcome_backup.R](09_mendelian_randomization/033_mr_gwas_finngen_outcome_backup.R) | GWAS exposure + FinnGen outcome MR template | exposure + FinnGen → MR table + plots | R · TwoSampleMR, qqman, RadialMR | R | scatter, forest, funnel, Manhattan, QQ, radial |
| 📄 | 043 | [043_MultimodalAD_MendelianRandomization_EN.R](09_mendelian_randomization/043_MultimodalAD_MendelianRandomization_EN.R) | AD multi-omics MR main analysis | AD GWAS vcf → MR + sensitivity | R · TwoSampleMR, gwasglue | R | scatter, forest, funnel, Manhattan, QQ |
| 📄 | 055 | [055_immunecell_disease_mr_directionality.R](09_mendelian_randomization/055_immunecell_disease_mr_directionality.R) | Immune-cell ↔ disease bidirectional MR + Steiger | exposure + outcome → MR + Steiger | R · TwoSampleMR, RadialMR | R | — |
| 📄 | 075 | [075_twosamplemr_coloc_drug_target.R](09_mendelian_randomization/075_twosamplemr_coloc_drug_target.R) | MR + colocalization drug-target evidence chain | exposure/outcome/locus → MR + coloc | R · TwoSampleMR, coloc | R | — |
| 📄 | 079 | [079_pqtl_mvmr_protein_mediation.R](09_mendelian_randomization/079_pqtl_mvmr_protein_mediation.R) | pQTL multivariable MR protein mediation | harmonised mvmr → MVMR mediation | R · (self-implemented MVMR) | R | — |
| ✅ | 499 | [499_lavaan_sem_mediation_path.R](09_mendelian_randomization/499_lavaan_sem_mediation_path.R) | SEM / path mediation with standardized-β diagram | composite scores → fit + path diagram | R · lavaan, semPlot | R | path-diagram |
| ✅ | 508 | [508_twostep_mediation_mr](09_mendelian_randomization/508_twostep_mediation_mr) | Two-step network mediation MR (Sobel/Delta/MC) | x + m instruments → mediation table | R · ggplot2 | R | path-diagram, forest |
| ✅ | 519 | [519_local_mr_pipeline](09_mendelian_randomization/519_local_mr_pipeline) | Fully local two-sample MR (no OpenGWAS API) | local exposure + outcome → estimates | R · TwoSampleMR, MRPRESSO, plinkbinr | R | scatter, forest, funnel, leave-one-out |

## 10 · TWAS (single-cell eQTL weights)

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 🔴 | 036 | [036_onek1k_twas_homo_hetero_fit.R](10_twas_sceqtl/036_onek1k_twas_homo_hetero_fit.R) | sc-eQTL homo/hetero elastic-net fit (step 1) | genotype + expr + cov → coefficients | R · glmnet, data.table | R | — |
| 🔴 | 037 | [037_onek1k_twas_component_fit.R](10_twas_sceqtl/037_onek1k_twas_component_fit.R) | Component-based prediction fit (step 2) | components → residuals + coef | R · glmnet, data.table | R | — |
| 🔴 | 038 | [038_onek1k_twas_weight_preprocess.R](10_twas_sceqtl/038_onek1k_twas_weight_preprocess.R) | Merge two-step coefficients for weights | step coefs → merged coefs | R · data.table | R | — |
| 🔴 | 039 | [039_onek1k_twas_fusion_weights.R](10_twas_sceqtl/039_onek1k_twas_fusion_weights.R) | Build per-cell-type FUSION weights | merged coefs → .RDat weights | R · data.table | R | — |
| 🔴 | 040 | [040_FUSION_TWAS_targetC.R](10_twas_sceqtl/040_FUSION_TWAS_targetC.R) | FUSION TWAS association (targetC weights) | GWAS sumstats + weights → TWAS | R · plink2R, glmnet | R | — |
| 🔴 | 041 | [041_FUSION_TWAS_S_targetC.R](10_twas_sceqtl/041_FUSION_TWAS_S_targetC.R) | FUSION TWAS association (S_targetC weights) | GWAS sumstats + weights → TWAS | R · plink2R, glmnet | R | — |
| 🔴 | 042 | [042_FUSION_TWAS_S_allC.R](10_twas_sceqtl/042_FUSION_TWAS_S_allC.R) | FUSION TWAS association (S_allC shared weights) | GWAS sumstats + weights → TWAS | R · plink2R, glmnet | R | — |

## 11 · WGCNA co-expression

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 054 | [054_wgcna_coexpression](11_wgcna/054_wgcna_coexpression) | Bulk WGCNA co-expression + module–trait | expr + traits → modules + figures | R · WGCNA, ComplexHeatmap | R | scale-free, dendrogram, module-trait heatmap |
| 🟡 | 504 | [504_hdwgcna_single_cell](11_wgcna/504_hdwgcna_single_cell) | hdWGCNA single-cell co-expression (metacell) | sc_counts.rds → modules + hubs | R · Seurat, hdWGCNA, igraph | R | soft-power, dendrogram, module feature-plot |

## 12 · TCGA prognosis (reference only)

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 048 | [048_tcga_single_gene_survival](12_tcga_prognosis/048_tcga_single_gene_survival) | Single-gene OS/DSS/DFI/PFI survival | gene_survival.csv → 4-endpoint KM | R · survival, survminer | R | KM (4 endpoints) |
| ✅ | 057 | [057_tcga_prognostic_risk_model](12_tcga_prognosis/057_tcga_prognostic_risk_model) | Prognostic risk model, five-figure panel | risk.csv → 5 figures + table | R · survival, timeROC, ComplexHeatmap | R | risk-plot, status, heatmap, KM, timeROC |
| ✅ | 060 | [060_tcga_immune_butterfly](12_tcga_prognosis/060_tcga_immune_butterfly) | Single-gene ↔ immune two-sided butterfly | expr + immune → correlation + plot | R · ggplot2 | R | butterfly (diverging) |
| 📦 | 497 | [497_scsurvival_cohort](12_tcga_prognosis/497_scsurvival_cohort) | scSurvival — single-cell cohort survival (vendored pkg) | sc cohort + survival → risk model | Py · PyTorch, scanpy, lifelines | Python | cohort-survival |

## 13 · Transcription-factor regulation / circos

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 🔴 | 047 | [047_rcistarget_tf_motif_network.R](13_tf_regulation_circos/047_rcistarget_tf_motif_network.R) | RcisTarget motif/TF enrichment + regulatory network | gene list + motif DB → network/Sankey | R · RcisTarget, igraph, visNetwork | R | network, sankey |
| ✅ | 053 | [053_circlize_chromosome_circos](13_tf_regulation_circos/053_circlize_chromosome_circos) | Plot genes onto chromosomes (genomic distribution) | gene_positions.csv → circos | R · circlize | R | circos |
| 🔴 | 081 | [081_pyscenic_regulon_tf_activity.py](13_tf_regulation_circos/081_pyscenic_regulon_tf_activity.py) | pySCENIC GRN + ctx + AUCell wrapper | expr/loom → regulons + aucell | Py · pyscenic (GRNBoost) | Python | — (downstream UMAP/heatmap) |
| ✅ | 511 | [511_tf_convergence_depmap_jaspar](13_tf_regulation_circos/511_tf_convergence_depmap_jaspar) | Three-evidence convergence to core TFs | tf_evidence.csv → convergence score | R · ggplot2, ggrepel | R | scatter, heatmap, lollipop |

## 14 · Single-cell in-silico perturbation

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 🟡 | 067 | [067_scperturb_etest.R](14_singlecell_perturbation/067_scperturb_etest.R) | scPerturb perturbation distance + E-test (light) | Seurat rds → edistance/etest | R · scperturbR | R | — |
| 🔴 | 068 | [068_gears_combo_perturbation.py](14_singlecell_perturbation/068_gears_combo_perturbation.py) | GEARS single/combo perturbation prediction (GPU) | h5ad + perturb list → predictions | Py · GEARS, torch | Python | — |
| 🔴 | 069 | [069_celloracle_grn_perturbation.py](14_singlecell_perturbation/069_celloracle_grn_perturbation.py) | CellOracle GRN virtual knockout (heavy) | Oracle pkl + genes → perturbed state | Py · celloracle | Python | — |
| 🔴 | 085 | [085_squidiff_diffusion_perturbation.py](14_singlecell_perturbation/085_squidiff_diffusion_perturbation.py) | Squidiff/PerturbDiff diffusion perturbation (GPU) | h5ad + config → predictions | Py · Squidiff, torch | Python | — |
| 🔴 | 494 | [494_genki_vgae_ko.py](14_singlecell_perturbation/494_genki_vgae_ko.py) | GenKI graph-VGAE virtual KO (KL ranking) | adata.h5ad + targets → KL ranking | Py · GenKI, torch-geometric | Python | — |
| ✅ | 495 | [495_bulkvgk_sctenifoldknk_di.R](14_singlecell_perturbation/495_bulkvgk_sctenifoldknk_di.R) | Bulk co-expression virtual KO + differential influence | two-group expr → per-gene DI ranking | R · igraph | R | scatter (DE vs DI) |
| 🟡 | 507 | [507_geneformer_insilico](14_singlecell_perturbation/507_geneformer_insilico) | Geneformer zero-shot embedding + in-silico deletion (baseline local) | counts/tokenized → KO ranking | Py · scanpy, sklearn (baseline); geneformer, torch | Python | UMAP, lollipop |

## 15 · Drug perturbation / repurposing

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 🔴 | 070 | [070_chemcpa_drug_perturbation.py](15_drug_perturbation/070_chemcpa_drug_perturbation.py) | chemCPA drug-perturbation expression prediction (GPU) | repo + config → train logs | Py · chemCPA, torch | Python | — |
| 🔴 | 071 | [071_scdrug_response_prediction.py](15_drug_perturbation/071_scdrug_response_prediction.py) | scDrug single-cell drug response (heavy) | 10x/h5ad → cluster drug response | Py · scDrug, GDSC/PRISM | Python | — |
| ✅ | 078 | [078_faers_pharmacovigilance](15_drug_perturbation/078_faers_pharmacovigilance) | FAERS disproportionality (ROR/PRR/BCPNN/EBGM) | reports/2×2 counts → signals.csv | R · ggplot2 | R | forest, heatmap |
| 🟡 | 518 | [518_beyondcell_drug_response](15_drug_perturbation/518_beyondcell_drug_response) | beyondcell core re-impl: BCS + therapeutic clusters | scRNA + drug signatures → BCS/ranking | R · UCell, ggplot2 | R | heatmap, lollipop, UMAP |

## 16 · Spatial communication / cell fate

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 🔴 | 072 | [072_cellrank_fate_drivers.py](16_spatial_communication/072_cellrank_fate_drivers.py) | scVelo + CellRank fate probabilities & drivers | h5ad (spliced) → driver genes | Py · scanpy, scvelo, cellrank | Python | — |
| 🔴 | 073 | [073_commot_spatial_communication.py](16_spatial_communication/073_commot_spatial_communication.py) | COMMOT spatial ligand–receptor communication | spatial h5ad + LR DB → comm. scores | Py · commot, scanpy | Python | — (downstream niche-map) |
| 🔴 | 074 | [074_tangram_sc_to_spatial.py](16_spatial_communication/074_tangram_sc_to_spatial.py) | Tangram single-cell → spatial mapping | sc + spatial h5ad → mapped h5ad | Py · tangram, scanpy | Python | — |
| 🟡 | 076 | [076_decoupler_tf_pathway_activity.py](16_spatial_communication/076_decoupler_tf_pathway_activity.py) | decoupler TF / pathway activity inference | h5ad → activity scores | Py · decoupler, scanpy | Python | — (downstream heatmap) |
| 🔴 | 077 | [077_nichenet_ligand_target.R](16_spatial_communication/077_nichenet_ligand_target.R) | NicheNet ligand activity + ligand–target links | receiver DE + prior → ligand activity | R · nichenetr | R | — (downstream heatmap) |
| 🔴 | 080 | [080_cell2location_squidpy_niche.py](16_spatial_communication/080_cell2location_squidpy_niche.py) | cell2location abundance + Squidpy neighborhood (GPU) | spatial h5ad + abundance → z-scores | Py · cell2location, squidpy | Python | — (downstream niche-map) |
| 🟡 | 505 | [505_spatial_advanced](16_spatial_communication/505_spatial_advanced) | Spatial advanced: RCTD deconv + NMF niche + interface degree | sc ref + spatial rds → fractions/niche | R · spacexr, RcppML, mistyR | R | niche-map (×3) |
| ✅ | 509 | [509_communication_functional_loop](16_spatial_communication/509_communication_functional_loop) | Communication functional loop: ligand→UCell→enrich→Venn | receptor expr + prior + group → consensus | R · UCell, ggplot2 | R | lollipop, violin, venn |
| 🟡 | 521 | [521_spatialglue_multiomics](16_spatial_communication/521_spatialglue_multiomics) | SpatialGlue spatial multi-omics domains (GNN; baseline local) | RNA + ADT grid → ARI + domains | Py · sklearn (baseline); SpatialGlue, torch-geometric | Python | spatial-scatter, lollipop |

## 17 · Advanced result figures

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 498 | [498_ggalluvial_sankey](17_advanced_figures/498_ggalluvial_sankey) | Multi-layer alluvial/Sankey flow | long table → alluvial figure | R · ggalluvial | R | sankey, alluvial |
| ✅ | 512 | [512_raincloud_plot](17_advanced_figures/512_raincloud_plot) | Raincloud (half-violin + box + jitter) vs bar charts | data.csv → stats + raincloud | R · ggdist | R | raincloud |
| ✅ | 513 | [513_ridgeline_plot](17_advanced_figures/513_ridgeline_plot) | Ridgeline distribution over an ordered factor | data.csv → summary + ridgeline | R · ggridges | R | ridgeline |
| ✅ | 514 | [514_dumbbell_slope_plot](17_advanced_figures/514_dumbbell_slope_plot) | Dumbbell + slope for paired change | data.csv → paired change + plots | R · ggrepel | R | dumbbell, slopegraph |
| ✅ | 515 | [515_chord_diagram](17_advanced_figures/515_chord_diagram) | Chord diagram for directed relations/flows | matrix.csv → flows + chord | R · circlize | R | chord |
| ✅ | 516 | [516_composite_multipanel](17_advanced_figures/516_composite_multipanel) | "Figure 1" multi-panel composite template | self-contained → composite figure | R · patchwork, ggrepel | R | composite (volcano+heatmap+forest+UMAP) |

> Also: `advanced_figure_tools.csv` (tool index), `download_advanced_figure_tools.ps1` (fetch external tools), `templates/` (closed-loop figure concept notes).

## 18 · External method sources

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 📦🗃️ | — | ai_scientific_figures/ | Vendored AutoFigure-Edit (ICLR'26): method text → editable SVG schematics. Local reference only, git-ignored. | — | — | — | — |

## 19 · Multi-omics integration & subtyping

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 🔴 | 083 | [083_mofa_diablo_multiomics.R](19_multiomics_integration/083_mofa_diablo_multiomics.R) | MOFA2 / mixOmics DIABLO multi-omics latent integration | multi-view matrices → factors/heatmaps | R · MOFA2 (reticulate), mixOmics | R | factor plot, heatmap |
| ✅ | 084 | [084_nmf_consensus_clustering](19_multiomics_integration/084_nmf_consensus_clustering) | NMF rank selection + consensus clustering subtyping | feature matrix → subtypes + figures | R · NMF, ConsensusClusterPlus, ComplexHeatmap | R | consensus, rank survey, subtype heatmap |

## 20 · Mutation / CNV / methylation / proteome / metabolome (templates)

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| 📄 | 522 | [522_mutation_maftools_pipeline.R](20_mutation_methylation_proteome/522_mutation_maftools_pipeline.R) | Somatic mutation summary template | MAF → oncoplot/summary | R · maftools | R | oncoplot |
| 📄 | 523 | [523_methylation_minfi_champ_pipeline.R](20_mutation_methylation_proteome/523_methylation_minfi_champ_pipeline.R) | Methylation differential analysis template | beta + meta → M-value dist + heatmap | R · limma, minfi, ChAMP | R | heatmap |
| 📄 | 524 | [524_proteomics_limma_msstats_pipeline.R](20_mutation_methylation_proteome/524_proteomics_limma_msstats_pipeline.R) | Proteomics differential analysis template | protein + meta → volcano + heatmap | R · limma, MSstats | R | volcano, heatmap |
| 📄 | 525 | [525_metabolomics_metaboanalystR_pipeline.R](20_mutation_methylation_proteome/525_metabolomics_metaboanalystR_pipeline.R) | Metabolomics differential analysis template | metabolite + meta → volcano + heatmap | R · MetaboAnalystR | R | volcano, heatmap |
| 📄 | 526 | [526_cnv_gistic_or_cnvkit_pipeline.md](20_mutation_methylation_proteome/526_cnv_gistic_or_cnvkit_pipeline.md) | CNV analysis entry note (GISTIC2/CNVkit/inferCNV) | — | — | md | — |

## 21 · Disease burden (GBD / NHANES / CHARLS / comorbidity)

Each sub-folder now ships a turnkey module grounded in the real cloned upstream
tools under `99_external_sources/` (which, with the topic/literature drafts, stays
local-only — git-ignored). Synthetic example data regenerates on run (not committed).

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 527 | [01_GBD/527_gbd_burden_trend](21_disease_burden_gbd/01_GBD/527_gbd_burden_trend) | GBD ASR trend + EAPC + Das Gupta decomposition + SDI | burden/pop/sdi csv → tables + figures | R · dplyr, ggplot2 | R | line-trend, forest, lollipop, diverging-lollipop, scatter |
| ✅ | 528 | [02_NHANES/528_nhanes_survey_weighted](21_disease_burden_gbd/02_NHANES/528_nhanes_survey_weighted) | NHANES survey-weighted means / regression / prevalence | nhanes.csv → svyglm + prevalence | R · survey, dplyr | R | dumbbell, forest, lollipop |
| ✅ | 529 | [03_CHARLS/529_charls_longitudinal_cohort](21_disease_burden_gbd/03_CHARLS/529_charls_longitudinal_cohort) | CHARLS trend + equipercentile equating + LMM + Cox/KM | panel.csv → trend/crosswalk/LMM/Cox | R · lme4, survival | R | line-trend, concordance, violin, forest, KM |
| ✅ | 530 | [04_comorbidity_network/530_comorbidity_network](21_disease_burden_gbd/04_comorbidity_network/530_comorbidity_network) | Disease-pair association → igraph → Louvain modules | patients.csv → network + metrics | R · igraph, ggraph | R | network, heatmap, lollipop |

## 22 · Single-cell metabolism

| St | # | Module | Purpose | Input → Output | Deps | Lang | Figures |
|----|---|--------|---------|----------------|------|------|---------|
| ✅ | 510 | [510_scmetabolism_pathway_activity](22_metabolism/510_scmetabolism_pathway_activity) | Single-cell metabolic pathway activity (AUCell/UCell-style) | expr + meta (+gmt) → activity + figures | R · ggplot2 (self-contained scorer) | R | dotplot, row-z heatmap, distribution |

---

## Figure type → module reverse index

Look up a figure you want; the listed modules can produce it.

- **Volcano** → 010, 516, 524, 525
- **Heatmap / ComplexHeatmap** → 010, 018, 022, 026, 034, 045, 054, 057, 059, 063, 083, 084, 492, 496, 502, 506, 510, 516, 518, 520
- **ROC** → 016, 034, 045, 052, 059, 063
- **Calibration / DCA / nomogram** → 016, 063
- **Forest** → 016, 032, 033, 043, 078, 503, 508
- **KM / survival** → 048, 057, 497(sc), 529
- **Line trend / time-series ribbon** → 527, 529
- **Concordance / crosswalk curve** → 529
- **Time-dependent ROC / risk-plot** → 057
- **UMAP / tSNE** → 026, 027, 044, 046, 049, 050, 058, 062, 504, 506, 507, 518
- **Violin (split)** → 026, 027, 044, 046, 049, 050, 509
- **Box** → 016, 021, 056, 492, 503
- **Bubble / dot plot** → 007, 022, 046, 049, 051, 510
- **Lollipop (bar replacement)** → 003, 005, 006, 011, 013, 014, 034, 035, 058, 496, 503, 509, 511, 518
- **Venn / UpSet** → 003, 005, 006, 011, 015, 034, 035, 509
- **PCA** → 010, 026, 027, 044, 056
- **Scatter** → 012, 013, 014, 032, 033, 043, 086, 495, 506, 511, 520, 521
- **Circos / chord** → 051, 053, 515
- **Network (igraph / visNetwork / ggraph)** → 007, 047, 530
- **Sankey / alluvial** → 047, 498
- **Correlation matrix** → 018, 021, 060, 492, 502
- **MR funnel / leave-one-out / Manhattan / QQ / radial** → 032, 033, 043, 519
- **SEM / mediation path diagram** → 499(lavaan), 508
- **Docking energy / MD curves (RMSD/RMSF/Rg/SASA/ΔG)** → 022, 086
- **Raincloud** → 512
- **Ridgeline** → 513
- **Dumbbell / slopegraph** → 514, 528
- **Butterfly (diverging)** → 060
- **Spatial feature / niche map** → 027, 050, 073, 080, 505, 521
- **Trajectory / pseudotime** → 044, 049, 050, 062, 082, 087
- **Vector field** → 062, 517
- **SHAP (beeswarm / dependence / waterfall / force)** → 052
- **Oncoplot** → 522
- **Composite multi-panel "Figure 1"** → 516
- **Stacked composition bar** (legitimate, not a ranking bar) → 021, 026, 044, 049, 492

---

## Cleanup changelog (resolved 2026-06-26)

All previously-flagged duplicates and numbering inconsistencies have been actioned:

| Item | Action taken |
|------|--------------|
| `06/019_immune_infiltration_source_generic.R` | **Deleted** — was byte-identical to 017 (only the header `# 编号` line differed) |
| `06/020_immune_infiltration_scoring_generic.R` | **Deleted** — was byte-identical to 018 |
| `06/021_immune_infiltration_viz.R` (loose) | **Deleted** — stale pre-modularization copy of the `021_immune_infiltration_viz/` directory module |
| `08/082_palantir_branch_probability.py` | **Renumbered → 087** (resolves the 082 collision with the Slingshot/tradeSeq/CytoTRACE2 module) |
| `09/497_lavaan_sem_mediation_path.R` | **Renumbered → 499** (resolves the cross-category 497 collision with 12's scSurvival) |
| `18/14_ai_scientific_figures/` | **Renamed → `ai_scientific_figures/`** (dropped the stale `14_` category-number prefix) |
| `20/*` five templates | **Renumbered → 522–526** to match the repo-wide `NNN_` convention |
| `04/045_multimodalad_ml_models.R` | **Marked TEMPLATE** in its header — it `source()`s a project-specific `refer.ML.R` helper that is not bundled, so it is a reference, not a turnkey run. For a turnkey multi-method run use 034; for prognostic combos use 059/496 |

> 045, 059, 496 are heavy/upstream-derived scripts that sit in category 04 but are
> really modelling/prognostic work — kept here for provenance; see status marks.
> New module numbers continue at **527+**.
