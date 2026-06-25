# 05 · Diagnostic models and validation

Build logistic diagnostic models from the feature genes selected in category 04, with internal and external validation.

| Module | Purpose | Language | Output figures |
|------|------|------|--------|
| [016 Diagnostic model](016_diagnostic_model_roc_calibration_dca/) | Modeling + internal evaluation | R | Nomogram, calibration, DCA, ROC, OR forest plot, boxplot |
| [063 External validation](063_geo_diagnostic_validation/) | Independent cohort validation | R | Training/validation ROC, calibration |

```bash
Rscript 016_diagnostic_model_roc_calibration_dca/016_diagnostic_model.R                      # Internal: modeling + evaluation
Rscript 063_geo_diagnostic_validation/063_diagnostic_validation.R                     # External: independent cohort validation
```

Follows the [unified framework conventions](../_framework/CONVENTIONS.md). Takes feature genes from category 04 as upstream input.
