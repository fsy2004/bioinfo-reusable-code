# 11 · WGCNA Co-expression Network

WGCNA co-expression network analysis: soft-threshold selection, module detection, and module-trait correlation.

| Module | Purpose | Language | Output figures |
|--------|---------|----------|----------------|
| [054 WGCNA co-expression network](054_wgcna_coexpression/) | Soft threshold, modules, module-trait | R | Scale-free fit, module dendrogram, module-trait heatmap |

## Usage

```bash
Rscript 054_wgcna_coexpression/054_WGCNA_coexpression.R --input data/expr.csv --traits data/traits.csv
```

Follows the [shared framework conventions](../_framework/CONVENTIONS.md). Key module genes can feed into 007 enrichment or 047 TF regulation.
