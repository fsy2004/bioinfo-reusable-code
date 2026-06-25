# 02 · GO / KEGG Enrichment Analysis

Annotate a set of candidate genes (differential genes, module genes, machine-learning features, target intersections, etc.) to GO functions and KEGG pathways.

| Module | Purpose | Language | Output figures |
|--------|---------|----------|----------------|
| [007 GO/KEGG enrichment](007_GO_KEGG富集分析/) | GO (BP/CC/MF) and KEGG pathway enrichment | R | GO faceted dot plot, KEGG lollipop, gene-pathway network |

```bash
Rscript 007_GO_KEGG富集分析/007_GO_KEGG_enrichment.R                 # 跑示例
Rscript 007_GO_KEGG富集分析/007_GO_KEGG_enrichment.R --input data/genes.csv
```

Follows the [shared framework conventions](../_framework/CONVENTIONS.md). Upstream is typically the 03 series (differential genes) or the 04 series (feature genes).
