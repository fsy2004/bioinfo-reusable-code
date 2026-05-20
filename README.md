# Reusable Bioinformatics Code Library

This repository is a curated local bioinformatics script library organized by analysis purpose.
It is intended to be reused across disease, immunology, transcriptomics, single-cell, spatial, genetics, and drug-response projects.

## Quick Entry

- Full overview: `README_代码库总览_2026-05-13.md`
- Script index: `全部代码总索引_按用途分类_2026-05-13.md`
- Machine-readable script index: `全部代码总索引_按用途分类_2026-05-13.csv`
- Rename map: `脚本重命名对照表_2026-05-13.csv`
- GitHub supplement candidates: `GitHub_生信代码补充候选_2026-05-13.md`

## Module Layout

Scripts are grouped under `按分析用途分类/`:

1. Network pharmacology and target databases
2. GO/KEGG enrichment
3. GEO transcriptome preparation and DEG analysis
4. Machine-learning feature selection
5. Diagnostic model validation
6. Immune infiltration and immune visualization
7. Molecular docking visualization
8. Single-cell, spatial transcriptomics, and trajectory analysis
9. Mendelian randomization and GWAS processing
10. TWAS and single-cell eQTL weights
11. WGCNA co-expression network
12. TCGA survival analysis examples
13. Transcription-factor regulation and circular genome plots
14. Single-cell virtual perturbation databases
15. Drug perturbation and repurposing
16. Spatial communication and cell fate
17. Advanced result figures and closed-loop visualization

## Reuse Notes

- Treat scripts as reusable templates, not one-click workflows.
- Check each script's required input paths, package dependencies, and phenotype labels before running it in a new project.
- Keep generated results, caches, and large omics files outside this repository unless they are small examples.
- Prefer adding project-specific wrappers in the target project and referencing this library as the upstream code source.
