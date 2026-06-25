# 005 · OMIM and GeneCards disease target Venn

Compute the intersection and union of multiple disease target lists (OMIM / GeneCards) and draw a Venn diagram and set-size bar chart.

## Input

An `--input` directory containing at least two disease target lists in CSV format. The gene column is detected automatically (`Gene` or the first column). Example data is provided in `example_data/` (OMIM and GeneCards disease target CSVs).

## Method

Compute set union and intersection, then render a Venn diagram via `venn_pub`, a set-size bar chart, and an UpSet plot. Uses the same engine as [003](../003_ctd_swiss_target_union_venn/).

## Usage

```bash
Rscript 005_disease_target_venn.R
```

## Outputs

- `results/`: intersection and union tables
- `assets/Target_Venn.png`: Venn diagram
- `assets/Set_size_bar.png`: set-size bar chart

The merged target set provides the full / high-confidence list of disease-related genes, for intersection with compound targets in module 006.

![Venn](assets/Target_Venn.png)

## Dependencies

R, `theme_pub`, `UpSetR`.

```r
install.packages("UpSetR")
```
