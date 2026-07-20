# 04 · 因果推断与遗传流行病 — Causal inference & genetics

本域共 25 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 工具变量准备 — Instrument preparation

- [028_mr_gwas_exposure_vcf_sig_snp.R](01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R) — Pick genome-wide significant SNPs from a VCF exposure
- [029_mr_gwas_exposure_ld_clump.R](01_instrument_prep/029_mr_gwas_exposure_ld_clump.R) — LD-clump to independent instruments
- [030_mr_gwas_exposure_add_eaf.R](01_instrument_prep/030_mr_gwas_exposure_add_eaf.R) — Add effect-allele frequency to exposure SNPs
- [031_mr_gwas_exposure_weak_iv_filter.R](01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R) — F-statistic weak-instrument filter

## 两样本孟德尔随机化 — Two-sample MR

- [032_mr_twosamplemr](02_two_sample_mr/032_mr_twosamplemr) — Self-contained MR causal inference (primary)
- [033_mr_gwas_finngen_outcome_backup.R](02_two_sample_mr/033_mr_gwas_finngen_outcome_backup.R) — GWAS exposure + FinnGen outcome MR template
- [043_MultimodalAD_MendelianRandomization_EN.R](02_two_sample_mr/043_MultimodalAD_MendelianRandomization_EN.R) — AD multi-omics MR main analysis
- [055_immunecell_disease_mr_directionality.R](02_two_sample_mr/055_immunecell_disease_mr_directionality.R) — Immune-cell ↔ disease bidirectional MR + Steiger
- [519_local_mr_pipeline](02_two_sample_mr/519_local_mr_pipeline) — Fully local two-sample MR (no OpenGWAS API)

## cis-MR 与药靶 — cis-MR & drug targets

- [075_twosamplemr_coloc_drug_target.R](03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R) — MR + colocalization drug-target evidence chain
- [535_mrbee_cis_mr](03_cis_mr_drug_target/535_mrbee_cis_mr) — MRBEE bias-corrected estimating-equation MR vs naive IVW
- [536_mrlink2_region_cis_mr](03_cis_mr_drug_target/536_mrlink2_region_cis_mr) — MR-link-2 single-region cis-MR (causal + pleiotropy) vs naive IVW

## 中介与多变量 MR — Mediation & MVMR

- [079_pqtl_mvmr_protein_mediation.R](04_mediation_mvmr/079_pqtl_mvmr_protein_mediation.R) — pQTL multivariable MR protein mediation
- [499_lavaan_sem_mediation_path.R](04_mediation_mvmr/499_lavaan_sem_mediation_path.R) — SEM / path mediation with standardized-β diagram
- [508_twostep_mediation_mr](04_mediation_mvmr/508_twostep_mediation_mr) — Two-step network mediation MR (Sobel/Delta/MC)
- [534_mvmr_cml_constrained](04_mediation_mvmr/534_mvmr_cml_constrained) — Constrained-ML multivariable MR (MVMR-cML-DP) vs IVW baseline

## 共定位 — Colocalization

- [537_sharepro_coloc](05_colocalization/537_sharepro_coloc) — SharePro effect-group colocalization vs classic single-causal coloc

## TWAS 与单细胞 eQTL — TWAS & sc-eQTL

- [036_onek1k_twas_homo_hetero_fit.R](06_twas_sceqtl/036_onek1k_twas_homo_hetero_fit.R) — sc-eQTL homo/hetero elastic-net fit (step 1)
- [037_onek1k_twas_component_fit.R](06_twas_sceqtl/037_onek1k_twas_component_fit.R) — Component-based prediction fit (step 2)
- [038_onek1k_twas_weight_preprocess.R](06_twas_sceqtl/038_onek1k_twas_weight_preprocess.R) — Merge two-step coefficients for weights
- [039_onek1k_twas_fusion_weights.R](06_twas_sceqtl/039_onek1k_twas_fusion_weights.R) — Build per-cell-type FUSION weights
- [040_FUSION_TWAS_targetC.R](06_twas_sceqtl/040_FUSION_TWAS_targetC.R) — FUSION TWAS association (targetC weights)
- [041_FUSION_TWAS_S_targetC.R](06_twas_sceqtl/041_FUSION_TWAS_S_targetC.R) — FUSION TWAS association (S_targetC weights)
- [042_FUSION_TWAS_S_allC.R](06_twas_sceqtl/042_FUSION_TWAS_S_allC.R) — FUSION TWAS association (S_allC shared weights)

## 稳健 MR 估计量 — Robust MR estimators

- [533_mrcare_winnerscurse_mr](07_robust_mr_methods/533_mrcare_winnerscurse_mr) — Winner's-curse-corrected MR (CARE/RIVW) vs naive baseline
