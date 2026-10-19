# 方法 → 论文 → GitHub → 模块

范围：纳入论文中逐项提取的方法，加上明确标为扩展的设计环节。每项资源都有实现状态；相关模块只是复用线索，不能据此声称该方法已完整实现。

七个新入口 604–610 已完成合成输入行为核查，均需显式输入/配置；其余旧模块保留历史状态，仅参考资源没有安装、运行或被包装成空壳模块。原始文献方法不构成质量背书。

机器读取用 [methods.json](methods.json)，来源版本/许可用 [github_sources.json](github_sources.json)。论文 ID 到 DOI、读取级别和参数边界分别见 [观察性文献](observational_papers.json) 与 [MR/GWAS 文献](mr_gwas_papers.json)。

状态：`new_adapter_synthetic_checked`=本次新增且合成接口已核查；`reference_only_with_related_entries`=已核来源与相关旧入口，具体方法仍需适配；`reference_only`=仅来源登记；其余状态按行解释。

## GWAS质控

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `gwas_regional_ld_audit` | 设计扩展 | [Cloufield/gwaslab](https://github.com/Cloufield/gwaslab) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc) | `new_adapter_synthetic_checked`；区域覆盖/build/祖源/效应与signed LD顺序审计；不运行GWAS或clumping。 |

## MR 主估计与诊断

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `bidirectional_two_sample_mr` | OBS2026_02 | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `cochran_q` | OBS2026_03 | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `egger_intercept` | OBS2026_03 | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `heterogeneity` | mr2026_glucocorticoid_nsclc | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `ivw` | OBS2026_02, OBS2026_03, mr2026_glucocorticoid_nsclc, background2021_cortisol_cancer | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `leave_one_out` | OBS2026_03, mr2026_glucocorticoid_nsclc | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `mr_egger` | OBS2026_02, OBS2026_03, mr2026_glucocorticoid_nsclc | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `reverse_mr` | DIRECT2024_01 | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `steiger` | DIRECT2024_01, mr2026_glucocorticoid_nsclc | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `two_sample_mr` | DIRECT2024_01, mr2026_glucocorticoid_nsclc, background2021_cortisol_cancer | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `two_sample_mr_unspecified` | OBS2026_01 | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 来源短文没有报告完整估计器/参数；不能填补为已核实IVW。 |
| `wald_ratio` | mr2026_glucocorticoid_nsclc | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_median` | OBS2026_02, OBS2026_03 | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr), [519](../../../modules/04_causal_inference_genetics/02_two_sample_mr/519_local_mr_pipeline), [055](../../../modules/04_causal_inference_genetics/02_two_sample_mr/055_immunecell_disease_mr_directionality.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## MRZ

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `mrz` | method2026_mrz | 未确定独立 GitHub；见论文/设计说明 | — | `code_not_identified`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 仅摘要/早期accepted版本支持；没有核实到作者代码，不实现估计器。 |

## SMR/HEIDI

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `smr_heidi` | mr2026_glucocorticoid_nsclc, mr2026_gastric_proteome | [JianYang-Lab/SMR](https://github.com/JianYang-Lab/SMR) | [075](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 SMR源码GPLv2与二进制MIT许可分别记录；SMR/HEIDI不等同Bayesian coloc。 |

## 交互与联合暴露

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `joint_exposure_four_groups` | OBS2026_05 | [cran/survey](https://github.com/cran/survey), [therneau/survival](https://github.com/therneau/survival) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `multiplicative_interaction` | OBS2026_05 | [cran/survey](https://github.com/cran/survey), [therneau/survival](https://github.com/therneau/survival) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `subgroup_analysis` | OBS2026_07, OBS2026_08 | [cran/survey](https://github.com/cran/survey), [therneau/survival](https://github.com/therneau/survival) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `subgroup_interaction` | OBS2026_01, OBS2026_03, OBS2026_04, OBS2026_05, OBS2026_06 | [cran/survey](https://github.com/cran/survey), [therneau/survival](https://github.com/therneau/survival) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 倾向评分

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `propensity_score_matching` | DIRECT2024_01 | [ngreifer/MatchIt](https://github.com/ngreifer/MatchIt) | [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 607实现PS权重诊断；匹配需另适配MatchIt，二者不等同。 |

## 偏倚敏感性

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `unmeasured_confounding_sensitivity` | 设计扩展 | [mayamathur/evalue_package](https://github.com/mayamathur/evalue_package) | — | `reference_only`；E-value不能校正混杂；OR/HR到RR的解释需要适用条件。 |

## 全基因组定位

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `bayesian_fdr` | method2026_coloc_design | [jianzeng/GCTB](https://github.com/jianzeng/GCTB), [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `functional_annotation` | method2026_gwfm | [jianzeng/GCTB](https://github.com/jianzeng/GCTB), [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `genome_wide_finemapping` | method2026_gwfm | [jianzeng/GCTB](https://github.com/jianzeng/GCTB), [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `local_global_credible_sets` | method2026_gwfm | [jianzeng/GCTB](https://github.com/jianzeng/GCTB), [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `polygenic_score` | preprint2026_breast_cross_ancestry | [jianzeng/GCTB](https://github.com/jianzeng/GCTB), [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `sbayesrc` | method2026_gwfm | [jianzeng/GCTB](https://github.com/jianzeng/GCTB), [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 决策曲线

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `decision_curve_analysis` | OBS2026_08 | [ddsjoberg/dcurves](https://github.com/ddsjoberg/dcurves) | [016](../../../modules/07_clinical_translational/01_diagnostic_models/016_diagnostic_model_roc_calibration_dca), [553](../../../modules/05_machine_learning/03_survival_ml/553_riskregression_dca_calibration) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 分子与药靶 MR

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `cis_eqtl_mr` | mr2026_glucocorticoid_nsclc | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [WSpiller/MVMR](https://github.com/WSpiller/MVMR) | [075](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R), [535](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/535_mrbee_cis_mr), [536](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/536_mrlink2_region_cis_mr), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `cis_pqtl_mr` | mr2026_gastric_proteome, preprint2026_breast_cross_ancestry | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [WSpiller/MVMR](https://github.com/WSpiller/MVMR) | [075](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R), [535](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/535_mrbee_cis_mr), [536](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/536_mrlink2_region_cis_mr), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `proteome_wide_mr` | mr2026_bcell_proteome | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [WSpiller/MVMR](https://github.com/WSpiller/MVMR) | [075](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R), [535](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/535_mrbee_cis_mr), [536](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/536_mrlink2_region_cis_mr), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `target_prioritization` | mr2026_bcell_proteome, mr2026_gastric_proteome | [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [WSpiller/MVMR](https://github.com/WSpiller/MVMR) | [075](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/075_twosamplemr_coloc_drug_target.R), [535](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/535_mrbee_cis_mr), [536](../../../modules/04_causal_inference_genetics/03_cis_mr_drug_target/536_mrlink2_region_cis_mr), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 加性交互

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `additive_interaction_reri_ap_si` | OBS2026_05 | [tunsmart/interactionR](https://github.com/tunsmart/interactionR) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 因果研究设计

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `dag_adjustment_sets` | 设计扩展 | [jtextor/dagitty](https://github.com/jtextor/dagitty), [r-causal/ggdag](https://github.com/r-causal/ggdag) | — | `reference_only`；先定义适应证、共同原因、碰撞点与中介，调整集来自假设。 |

## 复杂抽样

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `complex_survey_design` | OBS2026_02, OBS2026_08 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_chi_square` | OBS2026_06, OBS2026_07 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_cox` | OBS2026_07, OBS2026_08 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_descriptive_statistics` | OBS2026_06, OBS2026_07 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_linear_regression` | OBS2026_02 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_logistic_regression` | OBS2026_06, OBS2026_08, DIRECT2024_01 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `weighted_t_test` | OBS2026_07 | [cran/survey](https://github.com/cran/survey) | [528](../../../modules/17_population_epidemiology/02_complex_surveys/528_nhanes_survey_weighted), [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `survey_prevalence_and_mortality` | 设计扩展 | [cran/survey](https://github.com/cran/survey), [therneau/survival](https://github.com/therneau/survival) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `new_adapter_synthetic_checked`；主分析保留完整抽样设计；癌症史、死亡与新发结局分开。 |

## 复核与证据合成

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `cross_ancestry_meta_analysis` | preprint2026_breast_cross_ancestry | [wviechtb/metafor](https://github.com/wviechtb/metafor), [guido-s/meta](https://github.com/guido-s/meta) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `replication_meta_analysis` | DIRECT2024_01, background2021_cortisol_cancer | [wviechtb/metafor](https://github.com/wviechtb/metafor), [guido-s/meta](https://github.com/guido-s/meta) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 多变量MR

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `mvmr_conditional_strength` | 设计扩展 | [WSpiller/MVMR](https://github.com/WSpiller/MVMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract), [534](../../../modules/04_causal_inference_genetics/04_mediation_mvmr/534_mvmr_cml_constrained) | `reference_only`；需要联合工具和暴露关联协方差；610只审计已算条件F及来源。 |

## 多性状定位

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `multi_trait_finemapping` | method2026_mvsusie | [stephenslab/mvsusieR](https://github.com/stephenslab/mvsusieR) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `mvsusie` | method2026_mvsusie | [stephenslab/mvsusieR](https://github.com/stephenslab/mvsusieR) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 多祖源共定位

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `ecaviar` | method2026_maca | [fhormoz/caviar](https://github.com/fhormoz/caviar) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 已登记独立方法的作者源码，未运行；跨祖源LD与模型条件须逐方法核对。 |
| `mscaviar` | method2026_maca | [nlapier2/MsCAVIAR](https://github.com/nlapier2/MsCAVIAR) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 已登记独立方法的作者源码，未运行；跨祖源LD与模型条件须逐方法核对。 |
| `multi_ancestry_colocalization` | method2026_maca | [cthyshn/MACA](https://github.com/cthyshn/MACA) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `susiex` | method2026_maca | [getian107/SuSiEx](https://github.com/getian107/SuSiEx) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 已登记独立方法的作者源码，未运行；跨祖源LD与模型条件须逐方法核对。 |

## 多重检验

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `bonferroni` | OBS2026_02, DIRECT2024_01 | [wch/r-source](https://github.com/wch/r-source) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 工具变量

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `allele_harmonization` | OBS2026_03 | [MRCIEU/ieugwasr](https://github.com/MRCIEU/ieugwasr), [chrchang/plink-ng](https://github.com/chrchang/plink-ng), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [028](../../../modules/04_causal_inference_genetics/01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R), [029](../../../modules/04_causal_inference_genetics/01_instrument_prep/029_mr_gwas_exposure_ld_clump.R), [030](../../../modules/04_causal_inference_genetics/01_instrument_prep/030_mr_gwas_exposure_add_eaf.R), [031](../../../modules/04_causal_inference_genetics/01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R), [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `f_statistic` | OBS2026_03 | [MRCIEU/ieugwasr](https://github.com/MRCIEU/ieugwasr), [chrchang/plink-ng](https://github.com/chrchang/plink-ng), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [028](../../../modules/04_causal_inference_genetics/01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R), [029](../../../modules/04_causal_inference_genetics/01_instrument_prep/029_mr_gwas_exposure_ld_clump.R), [030](../../../modules/04_causal_inference_genetics/01_instrument_prep/030_mr_gwas_exposure_add_eaf.R), [031](../../../modules/04_causal_inference_genetics/01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R), [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `ld_clumping` | OBS2026_03 | [MRCIEU/ieugwasr](https://github.com/MRCIEU/ieugwasr), [chrchang/plink-ng](https://github.com/chrchang/plink-ng), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [028](../../../modules/04_causal_inference_genetics/01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R), [029](../../../modules/04_causal_inference_genetics/01_instrument_prep/029_mr_gwas_exposure_ld_clump.R), [030](../../../modules/04_causal_inference_genetics/01_instrument_prep/030_mr_gwas_exposure_add_eaf.R), [031](../../../modules/04_causal_inference_genetics/01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R), [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `snp_instrument_selection` | OBS2026_03 | [MRCIEU/ieugwasr](https://github.com/MRCIEU/ieugwasr), [chrchang/plink-ng](https://github.com/chrchang/plink-ng), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [028](../../../modules/04_causal_inference_genetics/01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R), [029](../../../modules/04_causal_inference_genetics/01_instrument_prep/029_mr_gwas_exposure_ld_clump.R), [030](../../../modules/04_causal_inference_genetics/01_instrument_prep/030_mr_gwas_exposure_add_eaf.R), [031](../../../modules/04_causal_inference_genetics/01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R), [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `summary_harmonization` | method2026_coloc_design | [MRCIEU/ieugwasr](https://github.com/MRCIEU/ieugwasr), [chrchang/plink-ng](https://github.com/chrchang/plink-ng), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [028](../../../modules/04_causal_inference_genetics/01_instrument_prep/028_mr_gwas_exposure_vcf_sig_snp.R), [029](../../../modules/04_causal_inference_genetics/01_instrument_prep/029_mr_gwas_exposure_ld_clump.R), [030](../../../modules/04_causal_inference_genetics/01_instrument_prep/030_mr_gwas_exposure_add_eaf.R), [031](../../../modules/04_causal_inference_genetics/01_instrument_prep/031_mr_gwas_exposure_weak_iv_filter.R), [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 常规统计与模型设计

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `chi_square` | OBS2026_04 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `complete_case_analysis` | OBS2026_06 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `complete_case_sensitivity` | OBS2026_05 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `likelihood_ratio_test` | OBS2026_04 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `mann_whitney_u` | OBS2026_04 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `multivariable_logistic` | OBS2026_01 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `progressive_covariate_models` | OBS2026_04, OBS2026_06, OBS2026_07 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `shapiro_wilk` | OBS2026_04 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `variance_inflation_factor` | OBS2026_04, OBS2026_07 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `winsorization_sensitivity` | OBS2026_08 | [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer), [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 机器学习算法

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `decision_tree_random_forest_xgboost_svm_neural_network_logistic_qda_lda_naive_bayes_knn` | OBS2026_07 | [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [034](../../../modules/05_machine_learning/01_feature_selection/034_multi_ml_feature_selection), [045](../../../modules/05_machine_learning/01_feature_selection/045_multimodalad_ml_models.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 框架入口不能代替各engine配置和跨验证泄漏检查；未新增全部算法实现。 |
| `logistic_ridge_lasso_elasticnet_tree_xgboost_lightgbm` | OBS2026_08 | [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [034](../../../modules/05_machine_learning/01_feature_selection/034_multi_ml_feature_selection), [045](../../../modules/05_machine_learning/01_feature_selection/045_multimodalad_ml_models.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 框架入口不能代替各engine配置和跨验证泄漏检查；未新增全部算法实现。 |

## 机器学习验证

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `grid_search` | OBS2026_07 | [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [034](../../../modules/05_machine_learning/01_feature_selection/034_multi_ml_feature_selection), [503](../../../modules/05_machine_learning/06_generalization_validation/503_generalization_robustness) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `ten_fold_cv` | OBS2026_07, OBS2026_08 | [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [034](../../../modules/05_machine_learning/01_feature_selection/034_multi_ml_feature_selection), [503](../../../modules/05_machine_learning/06_generalization_validation/503_generalization_robustness) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `train_validation_split` | OBS2026_07, OBS2026_08 | [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [034](../../../modules/05_machine_learning/01_feature_selection/034_multi_ml_feature_selection), [503](../../../modules/05_machine_learning/06_generalization_validation/503_generalization_robustness) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 模型解释

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `shap` | OBS2026_08 | [ModelOriented/DALEX](https://github.com/ModelOriented/DALEX) | [052](../../../modules/05_machine_learning/04_interpretability/052_shap_interpretation) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 特征选择

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `abess_feature_selection` | OBS2026_07 | [abess-team/abess](https://github.com/abess-team/abess) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `boruta_feature_selection` | OBS2026_07 | [cran/Boruta](https://github.com/cran/Boruta) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 生存与纵向

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `early_event_exclusion` | OBS2026_05 | [therneau/survival](https://github.com/therneau/survival) | [529](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/529_charls_longitudinal_cohort), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `kaplan_meier` | OBS2026_03, OBS2026_04, OBS2026_07 | [therneau/survival](https://github.com/therneau/survival) | [529](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/529_charls_longitudinal_cohort), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `log_rank` | OBS2026_04, OBS2026_07 | [therneau/survival](https://github.com/therneau/survival) | [529](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/529_charls_longitudinal_cohort), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `multivariable_cox` | OBS2026_01, OBS2026_03, OBS2026_04, OBS2026_05 | [therneau/survival](https://github.com/therneau/survival) | [529](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/529_charls_longitudinal_cohort), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `schoenfeld_test` | OBS2026_05, OBS2026_08 | [therneau/survival](https://github.com/therneau/survival) | [529](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/529_charls_longitudinal_cohort), [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `interval_censoring_and_first_cancer` | 设计扩展 | [therneau/survival](https://github.com/therneau/survival) | [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `new_adapter_synthetic_checked`；原因别Cox、区间删失AFT或规则person-period离散模型；不是自动选择模型。 |

## 稀疏结局

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `rare_event_firth` | 设计扩展 | [georgheinze/logistf](https://github.com/georgheinze/logistf) | — | `reference_only`；仅邻近方法补充；不声称已纳入论文使用，复杂抽样适配另核。 |

## 稳健 MR

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `mr_presso` | OBS2026_02 | [rondolab/MR-PRESSO](https://github.com/rondolab/MR-PRESSO) | [032](../../../modules/04_causal_inference_genetics/02_two_sample_mr/032_mr_twosamplemr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `mr_raps` | OBS2026_02, method2026_mrsimss | [qingyuanzhao/mr.raps](https://github.com/qingyuanzhao/mr.raps) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 稳健 MR 扩展

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `pleiotropy_mixture` | method2026_mrz | [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [595](../../../modules/04_causal_inference_genetics/07_robust_mr_methods/595_mreills_robust_mr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `robust_mr` | method2026_mrz | [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization) | [595](../../../modules/04_causal_inference_genetics/07_robust_mr_methods/595_mreills_robust_mr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 竞争风险

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `competing_risk_cif_fine_gray` | 设计扩展 | [cran/cmprsk](https://github.com/cran/cmprsk), [tagteam/riskRegression](https://github.com/tagteam/riskRegression) | [606](../../../modules/17_population_epidemiology/03_longitudinal_cohorts/606_longitudinal_incident_cancer) | `reference_only`；606原因别HR不能替代CIF/Fine-Gray亚分布HR。 |

## 精细定位与共定位

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `coloc` | method2026_maca | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `colocalization` | mr2026_glucocorticoid_nsclc, mr2026_bcell_proteome, mr2026_gastric_proteome, method2026_coloc_design, preprint2026_breast_cross_ancestry | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `coregulation_audit` | method2026_coloc_design | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `credible_sets` | method2026_mvsusie | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `pip` | method2026_mvsusie, method2026_gwfm | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `prior_sensitivity` | method2026_coloc_design | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `variant_coverage_qc` | method2026_coloc_design | [chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR) | [608](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/608_gwas_summary_qc), [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc), [537](../../../modules/04_causal_inference_genetics/05_colocalization/537_sharepro_coloc), [594](../../../modules/04_causal_inference_genetics/05_colocalization/594_colocboost_colocalization) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `susie_rss_multisignal_coloc` | 设计扩展 | [stephenslab/susieR](https://github.com/stephenslab/susieR), [chr1swallace/coloc](https://github.com/chr1swallace/coloc) | [609](../../../modules/04_causal_inference_genetics/08_gwas_finemapping/609_susie_coloc) | `new_adapter_synthetic_checked`；真实包接口的共享/独立信号与方向守卫检查，prior敏感性与H0-H4。 |

## 细胞语境证据

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `cell_type_context` | method2026_coloc_design | 未确定独立 GitHub；见论文/设计说明 | [046](../../../modules/01_single_cell/03_annotation_typing/046_scrna_publication_figures) | `existing_entry_or_design_reference`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `single_cell_context` | mr2026_bcell_proteome, mr2026_gastric_proteome | 未确定独立 GitHub；见论文/设计说明 | [046](../../../modules/01_single_cell/03_annotation_typing/046_scrna_publication_figures) | `existing_entry_or_design_reference`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 缺失与混杂

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `missingness_and_ps_diagnostics` | 设计扩展 | [amices/mice](https://github.com/amices/mice), [ngreifer/WeightIt](https://github.com/ngreifer/WeightIt), [ngreifer/cobalt](https://github.com/ngreifer/cobalt) | [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `new_adapter_synthetic_checked`；607实际运行MICE与base-R PS诊断；WeightIt/cobalt是已核参考，未安装。 |

## 缺失数据

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `multiple_imputation_mice` | OBS2026_05, OBS2026_08 | [amices/mice](https://github.com/amices/mice) | [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `rubin_pooling` | OBS2026_08 | [amices/mice](https://github.com/amices/mice) | [607](../../../modules/17_population_epidemiology/06_confounding_missingness/607_missingness_confounding) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 药物暴露

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `atc_drug_classification` | DIRECT2024_01 | [python/cpython](https://github.com/python/cpython) | [604](../../../modules/17_population_epidemiology/05_medication_exposure/604_prescription_exposure_audit) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 604只审计用户提供药物词典，不内置完整ATC数据库或临床滥用阈值。 |
| `exposure_outcome_timing_reconstruction` | DIRECT2024_01 | [python/cpython](https://github.com/python/cpython) | [604](../../../modules/17_population_epidemiology/05_medication_exposure/604_prescription_exposure_audit) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `exposure_record_audit` | 设计扩展 | [python/cpython](https://github.com/python/cpython) | [604](../../../modules/17_population_epidemiology/05_medication_exposure/604_prescription_exposure_audit) | `new_adapter_synthetic_checked`；药名/途径/时长/重复与缺失审计，使用不等于滥用。 |

## 药物流行病设计

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `target_trial_msm` | 设计扩展 | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation) | — | `reference_only`；新使用者、活动对照、时间零点与时变混杂，需要满足数据条件后适配。 |

## 跨队列设计

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `ancestry_stratified_analysis` | mr2026_gastric_proteome | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `assumption_assessment` | perspective2026_mr_future | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `estimand_definition` | perspective2026_mr_future | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `external_replication` | mr2026_glucocorticoid_nsclc, preprint2026_breast_cross_ancestry | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `multi_cohort_harmonization` | OBS2026_05 | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `multi_cohort_replication` | OBS2026_01, OBS2026_07 | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `power_planning` | method2026_gwfm | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `subtype_analysis` | preprint2026_breast_cross_ancestry | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `triangulation` | perspective2026_mr_future | [Causal-LDA/TrialEmulation](https://github.com/Causal-LDA/TrialEmulation), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 通路富集

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `pathway_enrichment` | OBS2026_01 | [YuLab-SMU/clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) | [051](../../../modules/02_spatial_transcriptomics/05_cell_communication/051_scrna_cellchat.R) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 遗传相关

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `genetic_correlation_ldsc` | DIRECT2024_01 | [CBIIT/ldsc](https://github.com/CBIIT/ldsc), [bulik/ldsc](https://github.com/bulik/ldsc) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 遗传研究设计

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `mr_exposure_construct_audit` | 设计扩展 | [WSpiller/MVMR](https://github.com/WSpiller/MVMR) | [610](../../../modules/04_causal_inference_genetics/01_instrument_prep/610_mr_exposure_contract) | `new_adapter_synthetic_checked`；区分用药倾向、内源皮质醇和分子工具；边际F不能替代条件F。 |

## 重叠与选择偏倚修正

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `mr_simss` | method2026_mrsimss | [amandaforde/mr.simss](https://github.com/amandaforde/mr.simss) | [533](../../../modules/04_causal_inference_genetics/07_robust_mr_methods/533_mrcare_winnerscurse_mr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `sample_overlap_correction` | method2026_mrsimss | [amandaforde/mr.simss](https://github.com/amandaforde/mr.simss) | [533](../../../modules/04_causal_inference_genetics/07_robust_mr_methods/533_mrcare_winnerscurse_mr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `simulation_benchmark` | method2026_mrsimss, method2026_mrz | [amandaforde/mr.simss](https://github.com/amandaforde/mr.simss) | [533](../../../modules/04_causal_inference_genetics/07_robust_mr_methods/533_mrcare_winnerscurse_mr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `winner_curse_correction` | method2026_mrsimss | [amandaforde/mr.simss](https://github.com/amandaforde/mr.simss) | [533](../../../modules/04_causal_inference_genetics/07_robust_mr_methods/533_mrcare_winnerscurse_mr) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 非线性与阈值

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `restricted_cubic_spline` | OBS2026_04, OBS2026_07, OBS2026_08, DIRECT2024_01 | [harrelfe/rms](https://github.com/harrelfe/rms), [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 605实现natural cubic spline基函数；未声称逐字复现论文RCS结点参数。 |
| `two_piecewise_regression` | OBS2026_04 | [harrelfe/rms](https://github.com/harrelfe/rms), [wch/r-source](https://github.com/wch/r-source) | [605](../../../modules/17_population_epidemiology/02_complex_surveys/605_survey_cancer_models) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 预测增益

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `nri_idi` | OBS2026_06 | [cran/nricens](https://github.com/cran/nricens), [cran/PredictABEL](https://github.com/cran/PredictABEL) | — | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 预测评估

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `brier_calibration` | OBS2026_08 | [xrobin/pROC](https://github.com/xrobin/pROC), [tagteam/riskRegression](https://github.com/tagteam/riskRegression), [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [016](../../../modules/07_clinical_translational/01_diagnostic_models/016_diagnostic_model_roc_calibration_dca), [553](../../../modules/05_machine_learning/03_survival_ml/553_riskregression_dca_calibration) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `pr_auc` | OBS2026_08 | [xrobin/pROC](https://github.com/xrobin/pROC), [tagteam/riskRegression](https://github.com/tagteam/riskRegression), [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [016](../../../modules/07_clinical_translational/01_diagnostic_models/016_diagnostic_model_roc_calibration_dca), [553](../../../modules/05_machine_learning/03_survival_ml/553_riskregression_dca_calibration) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `roc_auc` | OBS2026_06, OBS2026_07, OBS2026_08 | [xrobin/pROC](https://github.com/xrobin/pROC), [tagteam/riskRegression](https://github.com/tagteam/riskRegression), [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | [016](../../../modules/07_clinical_translational/01_diagnostic_models/016_diagnostic_model_roc_calibration_dca), [553](../../../modules/05_machine_learning/03_survival_ml/553_riskregression_dca_calibration) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |

## 风险评分与中介

| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |
|---|---|---|---|---|
| `mediation_analysis` | OBS2026_01 | [cran/lavaan](https://github.com/cran/lavaan), [therneau/survival](https://github.com/therneau/survival) | [499](../../../modules/04_causal_inference_genetics/04_mediation_mvmr/499_lavaan_sem_mediation_path.R), [057](../../../modules/07_clinical_translational/02_prognosis_survival/057_tcga_prognostic_risk_model) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
| `proteomic_metabolomic_risk_scores` | OBS2026_01 | [cran/lavaan](https://github.com/cran/lavaan), [therneau/survival](https://github.com/therneau/survival) | [499](../../../modules/04_causal_inference_genetics/04_mediation_mvmr/499_lavaan_sem_mediation_path.R), [057](../../../modules/07_clinical_translational/02_prognosis_survival/057_tcga_prognostic_risk_model) | `reference_only_with_related_entries`；GitHub实现是可选复用来源；未声称为每篇论文实际使用的包。相关旧模块需要具体项目适配，未在本次重跑。 |
