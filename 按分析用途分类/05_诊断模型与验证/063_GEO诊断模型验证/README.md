# 063 · GEO Diagnostic Model External Validation

External (independent-cohort) validation of a diagnostic model, comparing training and validation ROC curves and plotting a calibration curve on the validation set to assess generalization.

| | |
|---|---|
| **Language / Main dependencies** | R · `rms` `pROC` `ggplot2` |
| **Purpose** | External (independent-cohort) validation of a diagnostic model |
| **Input** | `example_data/` (train + validation matrices + gene list) |
| **Output** | `results/` AUC and figures; display figures in `assets/` |

## Input

| File | Required | Description |
|------|:---:|------|
| `--train` training matrix csv | Yes | First column is gene; sample-name suffix encodes group (`*_con`/`*_dis`) |
| `--valid` validation matrix csv | Optional | Same format, independent cohort; if omitted, the training set is self-evaluated |
| `--genes` diagnostic gene csv | Yes | Model genes |

## Method

`rms::lrm` fits a logistic model on the training cohort, predicts on the validation cohort, and `pROC` computes AUC for both cohorts and overlays the ROC curves. A validation-set calibration curve is drawn using quantile binning.

## Purpose

Complementary to 016 (internal evaluation): uses an independent cohort to test whether the diagnostic model is overfit and whether it generalizes, which is key evidence for publication.

## Features

- Runs on two matrices plus a gene list; automatically aligns the shared genes.
- Training vs validation ROC overlay (direct AUC comparison) plus a validation-set calibration curve.

## Outputs

| File | Figure type | Description |
|------|------|------|
| `assets/ROC_train_vs_valid.png` | ROC | Training/validation AUC comparison |
| `assets/Calibration_valid.png` | Calibration curve | Validation-set agreement |
| `results/AUC.csv` | Table | AUC per cohort |

![ROC](assets/ROC_train_vs_valid.png)

## Usage

```bash
Rscript 063_diagnostic_validation.R                                              # 示例
Rscript 063_diagnostic_validation.R --train data/train.csv --valid data/valid.csv --genes data/genes.csv
```

## Dependencies

```r
install.packages(c("rms","pROC","ggplot2"))
```
