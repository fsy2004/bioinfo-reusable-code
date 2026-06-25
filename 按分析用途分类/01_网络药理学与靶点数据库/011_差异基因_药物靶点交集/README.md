# 011 · DEG x Drug Target Intersection (Venn / UpSet)

Computes multi-set intersections of differentially expressed genes, drug targets, and disease targets, and renders Venn and UpSet plots.

## Input

A directory passed via `--input` containing multiple gene/target lists (CSV; the `Gene` column or first column is detected automatically). When 3 or more sets are present, an UpSet plot is also generated.

## Method

Compute the multi-set intersection (DEG ∩ drug targets ∩ disease targets), then draw the result with `venn_pub` (for 3 sets), UpSet, and a set-size bar chart.

## Usage

Identify core genes that are differentially expressed, are drug targets, and are disease-associated, for use as candidates in mechanistic studies and druggability assessment.

## Features

Provides two views (3-set Venn and UpSet). The Venn plot has no external dependencies.

## Outputs

| File | Plot type |
|------|-----------|
| `assets/Target_Venn.png` | 3-set Venn |
| `assets/Target_UpSet.png` | UpSet |
| `assets/Set_size_bar.png` | Set sizes |

![Venn](assets/Target_Venn.png)

## Run

```bash
Rscript 011_DEG_drug_target_venn.R
```

Dependencies: R, `theme_pub`, `UpSetR` (`install.packages("UpSetR")`).
