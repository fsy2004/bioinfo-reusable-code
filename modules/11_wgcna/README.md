# 11 · WGCNA Co-expression Network

WGCNA co-expression network analysis: soft-threshold selection, module detection, and module-trait correlation.

| Module | Purpose | Language | Output figures |
|--------|---------|----------|----------------|
| [054 WGCNA co-expression network](054_wgcna_coexpression/) | Soft threshold, modules, module-trait | R | Scale-free fit, module dendrogram, module-trait heatmap |
| [540 CWGCNA causal module inference 🟡](540_cwgcna_causal_module/) | 在 WGCNA 框架内判"模块→性状"还是"性状→模块"(双向中介);内置诚实基线对照普通模块-性状相关 | R | 基线相关 lollipop、因果方向 dumbbell、方向指数 lollipop、因果拓扑网络 |

## Usage

```bash
Rscript 054_wgcna_coexpression/054_WGCNA_coexpression.R --input data/expr.csv --traits data/traits.csv
```

Follows the [shared framework conventions](../_framework/CONVENTIONS.md). Key module genes can feed into 007 enrichment or 047 TF regulation.
