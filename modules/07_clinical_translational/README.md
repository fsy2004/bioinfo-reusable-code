# 07 · 临床与转化 — Clinical & translational

本域共 11 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 诊断模型 — Diagnostic models

- [016_diagnostic_model_roc_calibration_dca](01_diagnostic_models/016_diagnostic_model_roc_calibration_dca) — Logistic diagnostic model, full clinical evaluation
- [063_geo_diagnostic_validation](01_diagnostic_models/063_geo_diagnostic_validation) — External-cohort validation of a diagnostic model

## 预后与生存 — Prognosis & survival

- [048_tcga_single_gene_survival](02_prognosis_survival/048_tcga_single_gene_survival) — Single-gene OS/DSS/DFI/PFI survival
- [057_tcga_prognostic_risk_model](02_prognosis_survival/057_tcga_prognostic_risk_model) — Prognostic risk model, five-figure panel
- [060_tcga_immune_butterfly](02_prognosis_survival/060_tcga_immune_butterfly) — Single-gene ↔ immune two-sided butterfly

## 免疫浸润与解卷积 — Immune infiltration & deconvolution

- [017_immune_infiltration_source.R](03_immune_infiltration/017_immune_infiltration_source.R) — CIBERSORT deconvolution engine (SVR source)
- [018_immune_infiltration_scoring.R](03_immune_infiltration/018_immune_infiltration_scoring.R) — Immune infiltration scoring + cell/function matrices
- [021_immune_infiltration_viz](03_immune_infiltration/021_immune_infiltration_viz) — Fraction-matrix difference/composition/correlation viz
- [492_iobr_multimethod_deconvolution.R](03_immune_infiltration/492_iobr_multimethod_deconvolution.R) — IOBR 7-method deconvolution + cross-method consistency
- [520_bayesprism_deconvolution](03_immune_infiltration/520_bayesprism_deconvolution) — BayesPrism Bayesian deconvolution with ground-truth check

## 药物警戒 — Pharmacovigilance

- [078_faers_pharmacovigilance](04_pharmacovigilance/078_faers_pharmacovigilance) — FAERS disproportionality (ROR/PRR/BCPNN/EBGM)
