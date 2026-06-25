# 12 · TCGA tumor prognostic survival (reference)

Prognostic analysis for tumor cohorts (single gene or risk signature). Also applicable to non-tumor studies with follow-up data.

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [048 Single-gene multi-endpoint survival](048_TCGA单基因生存曲线/) | OS/DSS/DFI/PFI KM | R | KM for 4 endpoints | Done |
| [057 Prognostic risk model](057_TCGA预后风险模型/) | Risk signature, five outputs | R | Risk distribution, status, heatmap, KM, time-dependent ROC | Done |
| [060 Immune dual butterfly plot](060_TCGA免疫双蝴蝶图/) | Gene-immune correlation butterfly plot | R | butterfly | Done |
| 497 scSurvival | Joint single-cell and cohort survival (external package) | Python | Cohort survival | External package |

```bash
Rscript 048_TCGA单基因生存曲线/048_single_gene_survival.R --gene TP53   # single gene
Rscript 057_TCGA预后风险模型/057_prognostic_risk_model.R                  # risk signature
```

048 and 057 follow the [unified framework conventions](../_framework/CONVENTIONS.md). 497_scSurvival is an external Python package, kept for reference.
