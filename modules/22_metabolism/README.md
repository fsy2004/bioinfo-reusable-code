# 22 · Single-cell metabolism

Per-cell metabolic pathway activity scoring (scMetabolism / AUCell / UCell style)
and cell-type metabolic-preference comparison.

| Module | Purpose | Language | Output figures |
|--------|---------|----------|----------------|
| [510 scMetabolism pathway activity](510_scmetabolism_pathway_activity/) | Score KEGG/Reactome metabolic pathways per cell, compare across cell types | R | Activity dotplot, row-z heatmap, distribution |

## Usage

```bash
Rscript 510_scmetabolism_pathway_activity/510_scmetabolism_pathway_activity.R \
        --expr data/expr.csv --meta data/meta.csv
```

Self-contained scorer (no heavy dependency); pass an optional `.gmt` to use custom
gene sets. Follows the [shared framework conventions](../_framework/CONVENTIONS.md).
Pathway-activity scores pair naturally with 076 (decoupler activity) and the
single-cell modules in category 08.
