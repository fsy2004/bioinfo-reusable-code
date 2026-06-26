# 18 — External method sources (manifest only, not vendored)

Whole third-party source trees are not committed to this repository, because they bloat git (~19k files) and duplicate code that already lives upstream. This folder records the provenance of external methods collected locally; fetch the actual source from the upstream URL or pip when a method is needed. Any sub-folder here is git-ignored (`modules/18_external_sources/*/`) and kept as a local reference only.

## Local references (git-ignored, not committed)
- `14_ai_scientific_figures/` — vendored AutoFigure-Edit (ICLR'26): turns a method description into editable SVG schematics. Kept locally for reference; the upstream repo + paper are bundled but not tracked.

## Key upstream methods
Virtual perturbation: GEARS, CellOracle, scPerturb, scTenifoldKnk · Drug: chemCPA, scDrug ·
Trajectory/fate: CellRank, scVelo, Slingshot, tradeSeq, Palantir, CytoTRACE2 ·
Spatial: cell2location, RCTD/spacexr, BayesSpace, SpaGCN, Squidpy, stLearn, COMMOT, Tangram ·
Communication: NicheNet, MultiNicheNet · Regulation: SCENIC, decoupler ·
Multi-omics/subtype: MOFA2, mixOmics, NMF · Causal: TwoSampleMR, coloc, MVMR ·
Docking/MD: AutoDock Vina, GROMACS, gmx_MMPBSA, MDAnalysis · PV: FAERS pvda.

Runnable wrappers that call these (installed packages, not vendored source) live in
modules `14_singlecell_perturbation`, `15_drug_perturbation`, `16_spatial_communication`.
