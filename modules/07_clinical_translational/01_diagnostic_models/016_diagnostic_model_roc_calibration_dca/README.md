# 016 · Diagnostic Model — ROC / Calibration / DCA / Nomogram

Builds a logistic diagnostic model from an expression matrix and a diagnostic gene list, and produces a set of clinical evaluation plots (nomogram, calibration, DCA, ROC, OR forest, boxplot).

| | |
|---|---|
| Language / main dependencies | R · `rms` `rmda` `pROC` `ggplot2` |
| Purpose | Build and evaluate a gene-based diagnostic model from multiple angles |
| Input | `example_data/Sample_Type_Matrix.csv` + `diagnostic_genes.csv` |
| Output | Tables and figures in `results/`; display figures in `assets/` |

## Input

| File | Required | Description |
|------|:---:|------|
| `--input` expression matrix csv | yes | First column is gene; sample names carry a group suffix (`*_con` / `*_dis`) |
| `--genes` diagnostic gene csv | yes | First column is the selected diagnostic genes (typically from category 04 feature selection) |

## Method

`rms::lrm` multi-gene logistic regression; `nomogram` for the nomogram; `calibrate` (bootstrap) for the calibration curve; `rmda::decision_curve` for the decision curve; `pROC` for combined and single-gene ROC; `glm` to extract OR (95% CI) for the forest plot.

Method citations: Harrell, *rms* package; Vickers & Elkin, *Med Decis Making* 2006 (DCA).

## Usage

Maps selected feature genes into a usable diagnostic scoring model and evaluates it with ROC (discrimination), calibration (agreement), and DCA (clinical utility).

## Notes

- Runs the example without edits; detects groups automatically and handles complete separation.
- Outputs: nomogram, calibration curve, DCA, ROC (combined and single-gene), OR forest plot, and gene boxplot; base plots rendered with ggplot/theme_pub.

## Outputs

| File | Plot type | Description |
|------|------|------|
| `assets/Nomogram.png` | Nomogram | Risk score |
| `assets/Calibration.png` | Calibration curve | Predicted vs. observed agreement |
| `assets/DCA.png` | Decision curve | Clinical net benefit |
| `assets/ROC.png` | ROC | Combined model and single gene |
| `assets/OR_forest.png` · `Gene_boxplot.png` | Forest / boxplot | OR and expression difference |

![ROC](assets/ROC.png)
![DCA](assets/DCA.png)

## Run

```bash
Rscript 016_diagnostic_model.R                              # 示例
Rscript 016_diagnostic_model.R --input data/expr.csv --genes data/genes.csv --case _dis
```

## Dependencies

```r
install.packages(c("rms","rmda","pROC","ggplot2"))
```
