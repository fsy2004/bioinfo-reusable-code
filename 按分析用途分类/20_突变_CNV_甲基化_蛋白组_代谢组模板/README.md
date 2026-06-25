# 20 · Mutation / CNV / Methylation / Proteomics / Metabolomics templates

Command-line templates (argparse-style flags such as `--maf` / `--outdir`, no hardcoded paths) for standard single-omics analyses. Each template depends on its own Bioconductor packages and must be run in an environment where those packages are installed.

## Templates

| Template | Purpose | Main dependencies | Typical output figures |
|------|------|--------|-----------|
| `mutation_maftools_pipeline.R` | Somatic mutation summary | maftools | oncoplot, summary |
| `methylation_minfi_champ_pipeline.R` | Methylation differential analysis | minfi, ChAMP, limma | M-value distribution, heatmap |
| `proteomics_limma_msstats_pipeline.R` | Proteomics differential analysis | limma, MSstats | volcano, heatmap |
| `metabolomics_metaboanalystR_pipeline.R` | Metabolomics differential analysis | MetaboAnalystR | volcano, heatmap |
| `cnv_gistic_or_cnvkit_pipeline.md` | CNV analysis workflow notes | GISTIC2 / CNVKit / inferCNV | — |

## Usage

```bash
Rscript mutation_maftools_pipeline.R --maf cohort.maf --outdir results/mutation
```

## Suggested workflow

Multi-omics matrix, differential screening, [083 MOFA/DIABLO integration](../19_多组学整合_分型模板/), [084 NMF subtyping](../19_多组学整合_分型模板/084_NMF共识聚类分型/), then immune / spatial / prognostic interpretation.

## Dependencies

The Bioconductor packages used here (maftools, minfi, ChAMP, MSstats, MetaboAnalystR) are large with many dependencies, so no example figures are rendered locally. The template structure runs directly once the corresponding packages are installed. For volcano plots and heatmaps, the code in category 03 ([010](../03_GEO转录组整理与差异分析/)) and the ComplexHeatmap usage there can be reused.
