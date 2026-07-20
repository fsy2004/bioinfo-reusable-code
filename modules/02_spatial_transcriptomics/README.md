# 02 · 空间转录组 — Spatial transcriptomics

本域共 21 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 上游与细胞分割 — Pipeline & segmentation

- [027_spatial_seurat_auto.R](01_pipeline_segmentation/027_spatial_seurat_auto.R) — Spatial transcriptomics read/cluster/visualize
- [050_spatial_cluster_annot_trajectory.R](01_pipeline_segmentation/050_spatial_cluster_annot_trajectory.R) — Spatial cluster annotation + monocle3 pseudotime
- [573_proseg_cell_segmentation](01_pipeline_segmentation/573_proseg_cell_segmentation) — (用途待补)

## 空间域、空间可变基因与空间统计 — Domains, SVG & spatial statistics

- [505_spatial_advanced](02_domains_svg_stats/505_spatial_advanced) — Spatial advanced: RCTD deconv + NMF niche + interface degree
- [541_banksy_spatial_domains](02_domains_svg_stats/541_banksy_spatial_domains) — BANKSY neighbor-augmented spatial-domain segmentation vs non-spatial baseline
- [542_nnsvg_spatial_svg](02_domains_svg_stats/542_nnsvg_spatial_svg) — nnSVG spatially-variable genes (NNGP) vs non-spatial HVG baseline
- [543_squidpy_spatial_statistics](02_domains_svg_stats/543_squidpy_spatial_statistics) — squidpy spatial stats (Moran / nhood enrichment / co-occurrence / Ripley)
- [574_stair_spatial_integration](02_domains_svg_stats/574_stair_spatial_integration) — (用途待补)
- [575_scale_spatial_method](02_domains_svg_stats/575_scale_spatial_method) — (用途待补)

## 解卷积与单细胞映射 — Deconvolution & mapping

- [074_tangram_sc_to_spatial.py](03_deconvolution_mapping/074_tangram_sc_to_spatial.py) — Tangram single-cell → spatial mapping
- [080_cell2location_squidpy_niche.py](03_deconvolution_mapping/080_cell2location_squidpy_niche.py) — cell2location abundance + Squidpy neighborhood (GPU)
- [545_spotlight_deconvolution](03_deconvolution_mapping/545_spotlight_deconvolution) — SPOTlight (NMF+NNLS) spot cell-type deconvolution vs known-mix baseline

## 切片配准与三维重建 — Slice alignment & 3D

- [544_paste2_slice_alignment](04_alignment_3d/544_paste2_slice_alignment) — PASTE optimal-transport spatial slice alignment / 3D stacking

## 细胞通讯 — Cell-cell communication

- [051_scrna_cellchat.R](05_cell_communication/051_scrna_cellchat.R) — CellChat cell-communication network
- [073_commot_spatial_communication.py](05_cell_communication/073_commot_spatial_communication.py) — COMMOT spatial ligand–receptor communication
- [077_nichenet_ligand_target.R](05_cell_communication/077_nichenet_ligand_target.R) — NicheNet ligand activity + ligand–target links
- [509_communication_functional_loop](05_cell_communication/509_communication_functional_loop) — Communication functional loop: ligand→UCell→enrich→Venn
- [531_liana_consensus_cci](05_cell_communication/531_liana_consensus_cci) — LIANA+ rank-aggregate consensus cell-cell communication (6 methods)
- [576_cellnest_spatial_ccc](05_cell_communication/576_cellnest_spatial_ccc) — (用途待补)
- [577_spider_spatial_ccc](05_cell_communication/577_spider_spatial_ccc) — (用途待补)

## 空间多组学 — Spatial multi-omics

- [521_spatialglue_multiomics](06_spatial_multiomics/521_spatialglue_multiomics) — SpatialGlue spatial multi-omics domains (GNN; baseline local)
