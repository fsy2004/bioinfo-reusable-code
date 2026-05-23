# Reusable Bioinformatics Code Library

This repository is a curated local bioinformatics script library organized by analysis purpose.
It is intended to be reused across disease, immunology, transcriptomics, single-cell, spatial, genetics, drug-response, pharmacovigilance, and multi-omics projects.

## AI / New Project Entry

For future AI-assisted topic selection or project setup, read this first:

- [AI_新课题优先阅读指南_2026-05-24.md](<C:/Users/fsy/Desktop/原始代码/AI_新课题优先阅读指南_2026-05-24.md>)
- [方法指南_2026-05-24/生信分析方法与算法指南_2026-05-24.md](<C:/Users/fsy/Desktop/原始代码/方法指南_2026-05-24/生信分析方法与算法指南_2026-05-24.md>)
- [方法指南_2026-05-24/代码方法覆盖矩阵_2026-05-24.csv](<C:/Users/fsy/Desktop/原始代码/方法指南_2026-05-24/代码方法覆盖矩阵_2026-05-24.csv>)

High-frequency workflows:

- [workflow_geo_diagnosis.md](<C:/Users/fsy/Desktop/原始代码/方法指南_2026-05-24/workflows/workflow_geo_diagnosis.md>)
- [workflow_singlecell_mechanism.md](<C:/Users/fsy/Desktop/原始代码/方法指南_2026-05-24/workflows/workflow_singlecell_mechanism.md>)
- [workflow_mr_coloc_twas.md](<C:/Users/fsy/Desktop/原始代码/方法指南_2026-05-24/workflows/workflow_mr_coloc_twas.md>)

## Original Index Files

- Full overview: [README_代码库总览_2026-05-13.md](<C:/Users/fsy/Desktop/原始代码/README_代码库总览_2026-05-13.md>)
- Script index: [全部代码总索引_按用途分类_2026-05-13.md](<C:/Users/fsy/Desktop/原始代码/全部代码总索引_按用途分类_2026-05-13.md>)
- Machine-readable script index: [全部代码总索引_按用途分类_2026-05-13.csv](<C:/Users/fsy/Desktop/原始代码/全部代码总索引_按用途分类_2026-05-13.csv>)
- Rename map: [脚本重命名对照表_2026-05-13.csv](<C:/Users/fsy/Desktop/原始代码/脚本重命名对照表_2026-05-13.csv>)
- GitHub supplement candidates: [GitHub_生信代码补充候选_2026-05-13.md](<C:/Users/fsy/Desktop/原始代码/GitHub_生信代码补充候选_2026-05-13.md>)

## Module Layout

Scripts are grouped under [按分析用途分类](<C:/Users/fsy/Desktop/原始代码/按分析用途分类>):

1. Network pharmacology and target databases
2. GO/KEGG enrichment
3. GEO transcriptome preparation and DEG analysis
4. Machine-learning feature selection
5. Diagnostic model validation
6. Immune infiltration and immune visualization
7. Molecular docking and MD result visualization
8. Single-cell, spatial transcriptomics, and trajectory analysis
9. Mendelian randomization and GWAS processing
10. TWAS and single-cell eQTL weights
11. WGCNA co-expression network
12. TCGA survival analysis examples
13. Transcription-factor regulation, SCENIC wrappers, and circular genome plots
14. Single-cell virtual perturbation and perturbation databases
15. Drug perturbation, repurposing, and FAERS pharmacovigilance
16. Spatial communication, cell fate, NicheNet, and spatial niche wrappers
17. Advanced result figures and closed-loop visualization
18. External method source code waiting for selected integration
19. Multi-omics integration and subtype/pattern templates
20. Mutation, CNV, methylation, proteomics, and metabolomics templates

## Reuse Notes

- Treat scripts as reusable templates, not one-click workflows.
- Check each script's required input paths, package dependencies, and phenotype labels before running it in a new project.
- Keep generated results, caches, and large omics files outside this repository unless they are small examples.
- Prefer project-specific wrappers in the target project and reference this library as the upstream code source.
- Do not vendor entire third-party source trees into the main Git history unless license, size, and maintenance strategy are clear.
