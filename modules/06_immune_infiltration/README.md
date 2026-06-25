# 06 · Immune Infiltration and Immune Visualization

Estimate immune cell infiltration from an expression matrix, then visualize group differences and correlations.

## Modules

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [021 Immune visualization](021_immune_infiltration_viz/) | Proportion matrix to boxplot, stacked composition, correlation heatmap | R | Group boxplot, stacked composition, correlation heatmap | Available |
| 017 / 019 Deconvolution source functions | CIBERSORT / general deconvolution functions | R | None (engine) | Engine |
| 018 / 020 Immune scoring | Deconvolution, scoring, correlation | R | Correlation matrix | Engine |
| 492 IOBR multi-algorithm deconvolution | IOBR multi-method combination | R | Heatmap, correlation | Heavy environment |

## Pipeline

```
Expression matrix ── 017/019 deconvolution functions + 018/020 scoring ── cell proportion matrix ── 021 visualization (three figures)
```

## Notes

- 021: Takes a cell proportion matrix (standard CIBERSORT output) and produces the three figures with a single command.
- 017-020: Immune deconvolution engine. 017/019 are CIBERSORT source functions requiring the LM22 signature matrix plus raw expression; 018/020 are the scoring pipeline. Original scripts are kept as upstream engines.
- 492 IOBR: Requires installing IOBR from GitHub (dozens of dependencies). Not rendered locally; kept for reference.
- All modules follow the [unified framework conventions](../_framework/CONVENTIONS.md).
