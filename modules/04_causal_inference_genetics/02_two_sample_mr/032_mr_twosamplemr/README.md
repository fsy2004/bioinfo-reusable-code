# 032 · Mendelian Randomization (MR) Analysis

Causal inference from harmonized instrument data: IVW, MR-Egger, and weighted median estimates plus scatter, forest, funnel, and leave-one-out plots.

| | |
|---|---|
| **Language / main dependency** | R · `ggplot2` (MR core is self-contained, does **not** depend on TwoSampleMR) |
| **Purpose** | Causal inference: MR estimation of an exposure on an outcome with sensitivity analysis |
| **Input** | `example_data/harmonized_data.csv` |
| **Output** | Estimate tables and plots in `results/`; display figures in `assets/` |

## Input

Harmonized CSV containing: `SNP`, `beta.exposure`, `se.exposure`, `beta.outcome`, `se.outcome` (the output format of `TwoSampleMR::harmonise_data`; can be produced by the GWAS processing pipelines in modules 028-031).

## Method

- **IVW** (fixed effect): inverse-variance weighted regression through the origin, giving the primary causal estimate.
- **MR-Egger**: weighted regression with intercept; the intercept tests for directional pleiotropy.
- **Weighted median**: weighted median of the Wald ratios; robust to some invalid instruments.
- **Sensitivity**: leave-one-out analysis and funnel plot.

Method citation: Burgess *et al.*, *Eur J Epidemiol* 2017. The core MR implementation is self-contained, portable, and runs offline.

## Use

Infer the causal effect of an exposure on an outcome using genetic instruments (avoiding confounding and reverse causation), and test robustness with multiple methods plus sensitivity analysis.

## Outputs

| File | Plot | Description |
|------|------|------|
| `assets/MR_scatter.png` | Scatter | SNP effects with causal slopes from three methods |
| `assets/MR_forest.png` | Forest | Single-SNP Wald ratios with pooled estimate |
| `assets/MR_funnel.png` · `MR_leaveoneout.png` | Funnel / leave-one-out | Pleiotropy and robustness |

![scatter](assets/MR_scatter.png)

## Usage

```bash
Rscript 032_MR_analysis.R                                       # 示例
Rscript 032_MR_analysis.R --input data/harmonized.csv
```

## Dependencies

```r
install.packages("ggplot2")   # 核心 MR 为自包含实现,无需 TwoSampleMR
```
