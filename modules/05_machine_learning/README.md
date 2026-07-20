# 05 · 机器学习 — Machine learning

本域共 18 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 特征筛选 — Feature selection

- [012_lasso_feature_selection](01_feature_selection/012_lasso_feature_selection) — LASSO logistic feature selection
- [013_svm_rfe_feature_selection](01_feature_selection/013_svm_rfe_feature_selection) — SVM-RFE recursive elimination ranking
- [014_randomforest_feature_selection](01_feature_selection/014_randomforest_feature_selection) — Random-forest Gini importance ranking
- [015_ml_feature_intersection_venn_upset](01_feature_selection/015_ml_feature_intersection_venn_upset) — Intersect feature genes across methods
- [034_multi_ml_feature_selection](01_feature_selection/034_multi_ml_feature_selection) — caret 10-method models → AUC + consensus features
- [035_ml_combination_intersection](01_feature_selection/035_ml_combination_intersection) — Rank method combinations by intersection size
- [045_multimodalad_ml_models.R](01_feature_selection/045_multimodalad_ml_models.R) — Multi-ML integrated modelling (RSF/LASSO/GBM/BART…)
- [059_dual_disease_15ml_175combos.R](01_feature_selection/059_dual_disease_15ml_175combos.R) — Two-disease 15-ML × 175-combo screen + modelling
- [502_biomarker_triple_vote](01_feature_selection/502_biomarker_triple_vote) — Topology × correlation × Boruta triple-vote shortlist
- [554_rra_consensus_features](01_feature_selection/554_rra_consensus_features) — Robust Rank Aggregation consensus across feature-selection methods

## 分类模型 — Classification models

- [550_tabpfn_tabular_classifier](02_classification_models/550_tabpfn_tabular_classifier) — TabPFN foundation model vs LASSO/GBDT honest incremental eval

## 生存机器学习 — Survival ML

- [496_mime_101combo_prognostic.R](03_survival_ml/496_mime_101combo_prognostic.R) — Mime 10-algorithm × 101-combo prognostic signature (usage snippet)
- [551_aorsf_oblique_survival](03_survival_ml/551_aorsf_oblique_survival) — Oblique random survival forest (aorsf) vs CoxPH / standard RSF baseline
- [552_survex_survshap_explain](03_survival_ml/552_survex_survshap_explain) — survex time-dependent SurvSHAP(t) / SurvLIME explanation vs global baseline
- [553_riskregression_dca_calibration](03_survival_ml/553_riskregression_dca_calibration) — Honest survival-model eval: time-AUC + calibration + DCA + Brier/IBS

## 可解释性 — Interpretability

- [052_shap_interpretation](04_interpretability/052_shap_interpretation) — Train best model + SHAP interpretation

## 不确定性量化 — Uncertainty quantification

- [555_conformal_prediction_uq](05_uncertainty/555_conformal_prediction_uq) — Conformal prediction sets/intervals with finite-sample coverage vs naive baseline

## 泛化与外部验证 — Generalization & validation

- [503_generalization_robustness](06_generalization_validation/503_generalization_robustness) — Meta-analysis + LODO cross-cohort generalization
