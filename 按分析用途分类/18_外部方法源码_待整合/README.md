# 18 — External method sources (manifest only, not vendored)

Whole third-party source trees are not committed to this repository, because they bloat git (~19k files) and duplicate code that already lives upstream. This folder keeps only a provenance manifest of the external methods collected locally. Fetch the actual source from the upstream URL or pip when a method is needed.

## Manifest files
- `外部方法源码_pip来源清单_2026-05-21.csv` — pip / install provenance
- `外部方法源码清单_2026-05-21.partial.json` — source list (names, origins)
- `DocumentsGitHub移动结果_2026-05-24.csv`, `手动下载源码移动结果_2026-05-22.csv` — local move logs

## Key upstream methods (see manifest CSVs for exact URLs)
Virtual perturbation: GEARS, CellOracle, scPerturb, scTenifoldKnk · Drug: chemCPA, scDrug ·
Trajectory/fate: CellRank, scVelo, Slingshot, tradeSeq, Palantir, CytoTRACE2 ·
Spatial: cell2location, RCTD/spacexr, BayesSpace, SpaGCN, Squidpy, stLearn, COMMOT, Tangram ·
Communication: NicheNet, MultiNicheNet · Regulation: SCENIC, decoupler ·
Multi-omics/subtype: MOFA2, mixOmics, NMF · Causal: TwoSampleMR, coloc, MVMR ·
Docking/MD: AutoDock Vina, GROMACS, gmx_MMPBSA, MDAnalysis · PV: FAERS pvda.

Runnable wrappers that call these (installed packages, not vendored source) live in
modules `14_单细胞虚拟扰动_扰动数据库`, `15_药物扰动_药物重定位`, `16_空间通讯_细胞命运`.
