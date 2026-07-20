# 06 · Bulk 组学 — Bulk omics

本域共 19 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 差异表达 — Differential expression

- [008_geo_expression_matrix](01_differential_expression/008_geo_expression_matrix) — Annotate GEO probes → gene-level matrix
- [009_geo_sample_grouping](01_differential_expression/009_geo_sample_grouping) — Normalize matrix + attach group labels
- [010_geo_deg_volcano_heatmap_pca](01_differential_expression/010_geo_deg_volcano_heatmap_pca) — limma two-group DE with three figures
- [056_geo_multicohort_batch_correction](01_differential_expression/056_geo_multicohort_batch_correction) — Merge multi-cohort data + remove batch effect

## 富集分析 — Enrichment

- [007_go_kegg_enrichment](02_enrichment/007_go_kegg_enrichment) — GO/KEGG over-representation for a gene list
- [546_enrichplot_emap_cnet_tree](02_enrichment/546_enrichplot_emap_cnet_tree) — enrichGO + advanced plots replacing the plain bar (cnet/emap/tree)
- [549_goplot_chord_enrichment](02_enrichment/549_goplot_chord_enrichment) — GOplot gene×pathway chord / circle / membership heatmap

## 共表达网络(WGCNA 家族) — Co-expression networks

- [054_wgcna_coexpression](03_coexpression_networks/054_wgcna_coexpression) — Bulk WGCNA co-expression + module–trait
- [504_hdwgcna_single_cell](03_coexpression_networks/504_hdwgcna_single_cell) — hdWGCNA single-cell co-expression (metacell)
- [538_netrep_module_preservation](03_coexpression_networks/538_netrep_module_preservation) — NetRep permutation test of WGCNA module preservation across cohorts
- [539_smccnet_multiomics_network](03_coexpression_networks/539_smccnet_multiomics_network) — SmCCNet trait-driven sparse multi-omics network vs unsupervised baseline
- [540_cwgcna_causal_module](03_coexpression_networks/540_cwgcna_causal_module) — CWGCNA causal-direction (mediation) inference on WGCNA modules vs correlation baseline

## 多组学整合与分型 — Multi-omics integration

- [083_mofa_diablo_multiomics.R](04_multiomics_integration/083_mofa_diablo_multiomics.R) — MOFA2 / mixOmics DIABLO multi-omics latent integration
- [084_nmf_consensus_clustering](04_multiomics_integration/084_nmf_consensus_clustering) — NMF rank selection + consensus clustering subtyping

## 突变/甲基化/蛋白/代谢 — Mutation, methylation, proteome, metabolome

- [522_mutation_maftools_pipeline.R](05_mutation_methylation_proteome/522_mutation_maftools_pipeline.R) — Somatic mutation summary template
- [523_methylation_minfi_champ_pipeline.R](05_mutation_methylation_proteome/523_methylation_minfi_champ_pipeline.R) — Methylation differential analysis template
- [524_proteomics_limma_msstats_pipeline.R](05_mutation_methylation_proteome/524_proteomics_limma_msstats_pipeline.R) — Proteomics differential analysis template
- [525_metabolomics_metaboanalystR_pipeline.R](05_mutation_methylation_proteome/525_metabolomics_metaboanalystR_pipeline.R) — Metabolomics differential analysis template
- [526_cnv_gistic_or_cnvkit_pipeline.md](05_mutation_methylation_proteome/526_cnv_gistic_or_cnvkit_pipeline.md) — CNV analysis entry note (GISTIC2/CNVkit/inferCNV)
