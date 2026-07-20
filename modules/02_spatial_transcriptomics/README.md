# 02 · 空间转录组 — Spatial transcriptomics

本域共 24 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 上游与细胞分割 — Pipeline & segmentation

- [027_spatial_seurat_auto.R](01_pipeline_segmentation/027_spatial_seurat_auto.R) — Spatial transcriptomics read/cluster/visualize
- [050_spatial_cluster_annot_trajectory.R](01_pipeline_segmentation/050_spatial_cluster_annot_trajectory.R) — Spatial cluster annotation + monocle3 pseudotime
- [573_proseg_cell_segmentation](01_pipeline_segmentation/573_proseg_cell_segmentation) — 成像空间转录组转录本点云的细胞分割:自带「最近核外扩」可跑基线(半径扫描 + recall/precision/ambient-leak/ARI 评分),Proseg 本体为 Rust CLI 的守卫式命令行封装。

## 空间域、空间可变基因与空间统计 — Domains, SVG & spatial statistics

- [505_spatial_advanced](02_domains_svg_stats/505_spatial_advanced) — Spatial advanced: RCTD deconv + NMF niche + interface degree
- [541_banksy_spatial_domains](02_domains_svg_stats/541_banksy_spatial_domains) — BANKSY neighbor-augmented spatial-domain segmentation vs non-spatial baseline
- [542_nnsvg_spatial_svg](02_domains_svg_stats/542_nnsvg_spatial_svg) — nnSVG spatially-variable genes (NNGP) vs non-spatial HVG baseline
- [543_squidpy_spatial_statistics](02_domains_svg_stats/543_squidpy_spatial_statistics) — squidpy spatial stats (Moran / nhood enrichment / co-occurrence / Ripley)
- [574_stair_spatial_integration](02_domains_svg_stats/574_stair_spatial_integration) — 多切片空间转录组整合模块：本机可跑的三级整合阶梯基线（PCA / ComBat+PCA / 空间平滑+ComBat+PCA，scIB 式双轴评分），外加签名逐行核对自上游源码的 STAIR HGAT 守卫式封装
- [575_scale_spatial_method](02_domains_svg_stats/575_scale_spatial_method) — 空间组学多尺度空间域识别:空间平滑×Leiden分辨率网格上的跨种子稳定性搜索(朴素基线)+ SCALE 上游守卫式封装,含无空间信息对照。

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
- [576_cellnest_spatial_ccc](05_cell_communication/576_cellnest_spatial_ccc) — 空间转录组细胞通讯:空间受限 LR 共表达乘积基线(本机 CPU 即跑)+ CellNEST(GATv2 图注意力)守卫式封装,输出与 CellNEST 官方 9 列 schema 对齐
- [577_spider_spatial_ccc](05_cell_communication/577_spider_spatial_ccc) — 在相邻 spot 之间构建 interface(容量约束 Delaunay 图),对配体-受体对做 interface 层面 Moran's I 置换检验找空间可变互作(SVI),并与不建 interface 的 spot-level 共表达基线并排对照;官方 spider-st 包为守卫式可选路径。

## 空间多组学 — Spatial multi-omics

- [521_spatialglue_multiomics](06_spatial_multiomics/521_spatialglue_multiomics) — SpatialGlue spatial multi-omics domains (GNN; baseline local)
- [578_spatialex_omics_translation](06_spatial_multiomics/578_spatialex_omics_translation) — 用 H&E 形态学做锚,把一张切片测到的组学 panel 跨切片翻译到另一张切片(SpatialEx/SpatialEx+ 的 panel 对角整合),自带 Ridge+空间平滑的可跑基线与均值地板对照
- [579_simo_spatial_multiomics](06_spatial_multiomics/579_simo_spatial_multiomics) — 579 · SIMO —— 把无空间坐标的单细胞多组学(scRNA + 非转录组模态)通过最优传输概率性映射到空间转录组切片。模块并排跑三条路线:A 贪心相关性朴素基线(地板对照)、B 自写 POT fused-Gromov-Wasserstein 传输参照、C SIMO 正牌路线(守卫式,未装 simo-omics 即优雅退出打印安装命令,不静默降级)。合成数据自带 layer ground truth,四项外部指标(层准确率/中位位移/spot 占用率/单 spot 最大堆叠)全部落盘。

## 空间基础模型 — Spatial foundation models

- [580_novae_spatial_fm](07_foundation_models/580_novae_spatial_fm) — 580 · Novae — 多切片空间转录组的空间域/niche 划分与跨切片可迁移性评估。默认跑三级朴素基线阶梯(expression-only PCA+KMeans / 空间 kNN niche 平滑 / niche 平滑+逐切片 z-score),并用复刻自上游 novae/monitor/eval.py 的 FIDE / JSD / heuristic 三项指标打分,另加 ARI 与诊断列 ARI_celltype。Novae 本体走守卫式封装(--run-novae),包未装、基因名非真实 symbol 或拿不到 HuggingFace 权重时如实返回 skipped/failed,不伪造结果。
