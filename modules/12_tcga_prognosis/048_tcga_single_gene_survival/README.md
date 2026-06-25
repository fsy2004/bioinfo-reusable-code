# 048 · TCGA single-gene multi-endpoint survival curves

Generates high/low expression Kaplan-Meier curves for a single gene across the OS, DSS, DFI, and PFI endpoints from expression and survival data.

| | |
|---|---|
| **Language / main dependencies** | R · `survival` `survminer` |
| **Purpose** | Multi-endpoint KM evaluation of a single gene's prognostic value |
| **Input** | `example_data/gene_survival.csv` |
| **Output** | `results/` summary and `assets/` per-endpoint KM |

## Input

CSV containing the target gene expression column plus paired columns for each endpoint: `<EP>.time` (days) and `<EP>` (0/1), where EP is one of OS/DSS/DFI/PFI. Only the endpoints present are plotted.

## Method

Samples are split into high/low groups by the median gene expression. For each endpoint, `survfit` computes the KM estimate and `coxph` computes HR, 95% CI, and p value.

## Usage

```bash
Rscript 048_single_gene_survival.R                              # 示例
Rscript 048_single_gene_survival.R --input data/gene_survival.csv --gene TP53
```

## Outputs

Quickly assesses the prognostic relevance of a single gene across multiple survival endpoints, a standard TCGA pan-cancer prognostic analysis. Each endpoint produces an independent KM plot with HR/p and a risk table. Available endpoints are detected automatically.

| File | Plot type |
|------|------|
| `assets/KM_OS.png` / `KM_DSS.png` / `KM_DFI.png` / `KM_PFI.png` | Per-endpoint KM |

![OS](assets/KM_OS.png)

## Dependencies

```r
install.packages(c("survival","survminer"))
```
