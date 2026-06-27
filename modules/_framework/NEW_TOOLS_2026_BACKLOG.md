# New-tools-2026 backlog — turnkey-module candidates

Scan date **2026-06-26**. Source: a 7-line gap analysis of the 2026-H1 tool survey
(`bioinfo-DL-library/analysis-tools-2026/`) cross-checked **module-by-module against
[`../CATALOG.md`](../CATALOG.md)** to honor the reuse-not-rewrite rule. Every candidate
below ① has a **verified installable R/Python package** (install command web-checked on
the scan date), ② has **no equivalent module** already in the catalog, and ③ produces
**bar-free figures** (lollipop/dot/dumbbell/violin/raincloud/heatmap/network/forest/
scatter/contour), per the repo figure rule.

**30 distinct candidates**, suggested numbers **527–556** (next free after 526).
Status legend identical to CATALOG: ✅ turnkey-local · 🟡 core/baseline local + real
package on server. None need GPU. **Build order priority** is marked ⭐ (highest-value,
CRAN/Bioc, easiest to verify locally first).

> Two tools surfaced on two lines each and are **merged to one module**: ProLIF
> (netpharm + docking_md → 547) and the free-energy-landscape figure (netpharm FEL +
> docking_md bio3d FEL/PCA/DCCM → 548, using the richer `bio3d` route).

---

## 08 · Single-cell / spatial / trajectory

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 527 ⭐ | sccomp + voomCLR composition test | R | ✅ | `BiocManager::install("sccomp")`; `remotes::install_github("koenvandenberge/voomCLR")` | No single-cell **differential-abundance/composition** module (021=bulk fractions, 084=subtyping). Operationalizes 7th iron rule (anti-pseudoreplication). | composition box + per-sample lines, credible-effect lollipop, proportion raincloud | sccomp PNAS 2023; voomCLR Bioinformatics 2026 (btaf637) |
| 528 ⭐ | miloR neighborhood DA | R | ✅ | `BiocManager::install("miloR")` | KNN-neighborhood DA without discrete clusters; the signature **DA beeswarm** has no producing module. | DA beeswarm (lineage-ordered, color=logFC), neighborhood network, neighborhood volcano | miloR Nat Biotechnol 2022 |
| 530 | CopyKAT scRNA CNV/aneuploidy | R | 🟡 | `remotes::install_github("navinlabcode/copykat")` | inferCNV **deprecated**; no scRNA-CNV module (526=bulk GISTIC/CNVkit). RAM-heavy on real data → synthetic-core local, full run on server. | CNV heatmap (cells×bins, aneuploid/diploid split), UMAP by ploidy, subclone-size lollipop | CopyKAT Nat Biotechnol 2021 |
| 541 ⭐ | BANKSY spatial domain ID | R/Python | ✅ | `BiocManager::install("Banksy")` / `pip install banksy_py` | No neighbor-augmented spatial-**domain** ID (027/050 = read/cluster only). 2026 non-DL CPU anchor. | domain-segmentation spatial scatter, neighbor-feature UMAP, domain-marker violin | Nat Genet 2024 (s41588-024-01664-3) |
| 542 ⭐ | Castl + SPARK-X spatially variable genes | R | ✅ | `remotes::install_github("TheY11/Castl")`; baseline `remotes::install_github("xzhoulab/SPARK")`; `BiocManager::install("nnSVG")` | **Zero SVG module** in catalog though SVG is core spatial analysis. Castl=2026 ensemble, SPARK-X=honest CPU baseline. | SVG rank lollipop/dot (Moran's I / aggregated p), top-SVG spatial feature-map, method dumbbell, FDR scatter | Castl Brief Bioinform 2026-03; SPARK-X Genome Biol 2021 |
| 543 ⭐ | squidpy spatial statistics | Python | ✅ | `pip install squidpy` | No spatial-stats module (510=metabolism). Standard Moran/Geary/nhood-enrichment/co-occurrence/Ripley toolkit. | nhood-enrichment z-score heatmap, Moran's I lollipop/scatter, co-occurrence curve, Ripley's L curve | Nat Methods 2021 (s41592-021-01358-2) |
| 544 | PASTE2 / STalign slice alignment + 3D | Python | ✅ | `pip install paste-bio`; `pip install git+https://github.com/JEFworks-Lab/STalign` | No multi-slice **registration/3D-reconstruction** module. CPU baselines from 2026 Nat Comput Sci registration benchmark. | before/after aligned-spot scatter, 3D stacked reconstruction, displacement-field arrows | PASTE2 Nat Methods 2022; STalign Nat Commun 2023 |

## 03 · Transcriptomics & DE

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 529 ⭐ | muscat pseudobulk multi-sample DS | R | ✅ | `BiocManager::install("muscat")` | 2026 **gold-standard pseudobulk** DS absent (008-010=bulk limma only). Sample-level aggregation → no I-type inflation; directly enacts anti-pseudoreplication rule. | pseudobulk MDS/PCA by sample, pathway-guided volcano per cell type, row-z DE heatmap, top-DS lollipop | muscat Nat Commun 2020 |

## 09 · Mendelian randomization (all summary-data, fully local, no OpenGWAS API)

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 533 ⭐ | MRcare (CARE/RIVW) | R | ✅ | `remotes::install_github("ChongWuLab/MRcare")` | No **winner's-curse-corrected** estimator (032/519=IVW/Egger/median/PRESSO only). Most substantive 2026 robust-MR advance. | per-method estimate lollipop, CARE-vs-RIVW-vs-IVW forest, selection-bias-corrected slope scatter | JASA 2026 (PMID 41869282) |
| 534 ⭐ | MendelianRandomization mr_mvcML (MVMR-cML-DP) | R | ✅ | `install.packages("MendelianRandomization")` (v0.10.0) | Upgrades self-coded MVMR (079) to constrained-ML, robust to correlated+uncorrelated pleiotropy. | MVMR-cML direct-effect forest, MVMR-cML-vs-IVW dumbbell, exposure×outcome β heatmap | CRAN v0.10.0; MVMR-cML method |
| 535 | MRBEEX (cis-MRBEE) | R | ✅ | `remotes::install_github("harryyiheyang/MRBEEX")` | Fine-map-anchored weak-IV + pleiotropy-robust **multivariable cis-MR** (drug-target workhorse), absent (075=MR+coloc, 079=basic MVMR). | cis multivariable effect lollipop, MRBEE-vs-naive forest, fine-mapped-weight scatter | Brief Bioinform 2025 (bbaf250) |
| 536 | MR-link-2 (mrlink2) | Python | ✅ | `pip install git+https://github.com/adriaan-vd-graaf/mrlink2` | **Single-region cis-MR** (where classic MR Type-I error explodes); jointly estimates causal+pleiotropy. CPU numpy/scipy. | per-region causal-vs-pleiotropy scatter, region-wise lollipop, vs-single-region-IVW forest | Nat Commun 2025 (PMID 40610416) |
| 537 | SharePro_coloc | Python | ✅ | `pip install git+https://github.com/zhwm/SharePro_coloc` | Effect-group variational coloc (sensitivity method, beats coloc+SuSiE under multiple causal variants); 075 uses single-causal coloc. | locuscompare-style scatter (effect-group colored), per-group PP heatmap, PP lollipop | Bioinformatics 2024 (btae295) |

## 11 · WGCNA co-expression

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 538 ⭐ | NetRep module-preservation test | R | ✅ | `install.packages("NetRep")` | No **cross-dataset module preservation/reproducibility** permutation test; pairs with 054/504 for external validation. | preservation scatter (Zsummary/medianRank vs size + thresholds), permutation null violin, per-module p lollipop | NetRep Cell Systems 2016; CRAN 1.2.10 |
| 540 | CWGCNA module causal direction | R | 🟡 | `devtools::install_github("yuabrahamliu/CWGCNA")` | Module→trait vs trait→module **causal/mediation** inside WGCNA (answers correlation≠causation); 508/499 do MR/SEM on scores, not modules. | causal-direction forest/dumbbell, mediation path diagram, topology-feature heatmap | NAR Genom Bioinform 2024-06 (PMID 38666214) |

## 19 · Multi-omics integration

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 539 ⭐ | SmCCNet trait-driven sparse-CCA subnetwork | R | ✅ | `install.packages("SmCCNet")` | Phenotype-driven **sparse-CCA cross-omics subnetwork** — distinct from 083 (MOFA/DIABLO latent-factor). | trait-specific cross-omics network (igraph/ggraph), selected-subnetwork adjacency heatmap, top canonical-weight lollipop | SmCCNet 2.0 BMC Bioinformatics 2024; CRAN 2.0.7 |

## 16 · Spatial communication

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 531 ⭐ | LIANA+ consensus cell-communication | Python | ✅ | `pip install liana` | Catalog has single-method comm only (051 CellChat, 077 NicheNet). LIANA+ = consensus meta-framework (CellPhoneDB/CellChat/NATMI/Connectome/SCSignalR + rank-aggregate). | consensus L-R dotplot (specificity×magnitude), cross-condition tile heatmap, comm network, source-target chord | liana Nat Cell Biol 2024 (PMC11392821) |

## 02 · Enrichment

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 546 ⭐ | enrichplot emapplot + treeplot | R | ✅ | `BiocManager::install(c("enrichplot","clusterProfiler","ggtangle"))` | 007 stops at dotplot/network; **emapplot** (pathway-module map) & **treeplot** (term tree) are the flagship reverse-bar enrichment figures, missing. | emapplot (pathway-module network), treeplot (clustered term tree), circular cnetplot | enrichplot Bioc; ggtangle ≥0.0.5 |
| 549 | GOplot GOCircle / GOCluster | R | ✅ | `install.packages("GOplot")` | GOCircle's z-score + per-gene logFC circle is a distinct enrichment figure no module produces (515=generic chord). | GOCircle (term ring + logFC scatter + z-score), GOCluster (dendro+circle), GOChord | GOplot CRAN 1.0.2 |

## 07 · Molecular docking & dynamics (analysis-side figures; engines stay in 086)

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 547 ⭐ | ProLIF interaction fingerprint *(merged netpharm+docking)* | Python | ✅ | `pip install prolif` | Decodes docking/MD poses into residue-level H-bond/π/hydrophobic fingerprint; 022/086 cannot. Replaces contact-count bars. | interaction-fingerprint barcode/heatmap (residue×type), per-residue contact lollipop, pose-similarity heatmap, ligand-network | ProLIF J Cheminform 2021; v2.2.0 |
| 548 ⭐ | bio3d FEL / PCA / DCCM conformational landscape *(merged)* | R | ✅ | `install.packages("bio3d")` | 086 emits only RMSD/RMSF/Rg/SASA; the **FEL contour / PCA / DCCM** top-journal conformational figures are absent. Pure Boltzmann-inverted 2D histogram (no GROMACS at render). | FEL 2D contour (PC1×PC2, −kT lnP), PCA scatter + porcupine, DCCM cross-correlation heatmap | bio3d Bioinformatics; CRAN 2.4-5 |
| 556 ⭐ | PoseBusters pose validity panel | Python | ✅ | `pip install posebusters` | Encodes the docking IRON RULE (DL pose generator + PoseBusters gatekeeper, report PB-valid %). No physical-validity figure exists. | PB-valid pass/fail tick heatmap (pose×check), per-check pass-rate lollipop, DL-vs-physics validity dumbbell | PoseBusters Chem Sci 2024 (D3SC04185A) |

## 17 · Advanced figures

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 532 | SCpubr publication-grade wrapper | R | ✅ | `install.packages("SCpubr")` | Single-call publication-ready single-cell figures (color-blind-safe). Standardizes aesthetic; partial overlap 046/17. | UMAP, dotplot, enrichment, communication, CNV, alluvial — one call | SCpubr Bioinform Adv 2026 (vbag151) |
| 545 | scatterpie / SPOTlight spatial scatterpie | R | ✅ | `install.packages("scatterpie")`; `BiocManager::install("SPOTlight")` | The standard bar-replacement deconvolution figure (per-spot pie at xy) is absent (505 niche-map ≠ scatterpie). | scatterpie (per-spot proportion pie on coords), scatterbar/PieGlyph variant | scatterpie CRAN 0.2.6 |

## 05 · Diagnostic models / 12 · Prognosis

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 550 ⭐ | TabPFN-2.5 tabular foundation model | Python | ✅ | `pip install tabpfn` | Only tool with hard **"beats default XGBoost"** evidence; CPU sub-second on ≲1000 samples (DL-strategy route A). Must ship LASSO/XGBoost honest baseline + in-fold prefilter (anti-leakage). | overlaid ROC+PR, calibration curve, SHAP/permutation lollipop, confusion heatmap (3-way vs LASSO/XGB) | TabPFN v2 Nature 2025; 2.5 arXiv:2511.08667 |
| 551 | aorsf accelerated oblique RSF | R | ✅ | `install.packages("aorsf")` | Oblique-split RSF, more accurate + interpretable than standard RSF (045/059/496); 057=Cox. Must compare regularized Cox. | importance lollipop/beeswarm, risk-stratified KM, time-dependent AUC | aorsf; Sci Rep 2025 |
| 552 | survex + SurvSHAP(t) survival explainability | R/Python | ✅ | `install.packages("survex")`; `pip install survshap` | Extends classification-only SHAP (052) to **time-dependent** survival explanation (RSF/aorsf downstream). | SurvSHAP(t) time-dependent contribution curves, SurvLIME local, importance lollipop | SurvSHAP(t) Knowl-Based Syst 2023 |
| 553 ⭐ | riskRegression honest survival evaluation | R | ✅ | `install.packages("riskRegression")` | Aligns with NeurIPS 2025 "Stop Chasing the C-index": adds **IBS + D-calibration + time-AUC** (016=binary DCA, 057=single-model timeROC). | time-dependent AUC/Brier curves, calibration curve, IBS dot/line, multi-model comparison | riskRegression CRAN 2026.03 |

## 04 · ML feature selection

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 554 | RobustRankAggreg consensus / stability FS | R | ✅ | `install.packages("RobustRankAggreg")` | Probabilistic **rank-aggregation** consensus (gives p) + cross-resample stability (Jaccard/Kuncheva) — distinct from pure intersection/vote (015/035/502). | multi-filter rank dumbbell/slopegraph, stability Jaccard heatmap, consensus lollipop | RRA; Metabolites 2025/26, Allergy 2026 |

## 23 · Uncertainty / conformal *(new category)*

| # | Module | Lang | St | Install | Fills gap vs catalog | Figures (bar-free) | Pub |
|---|--------|------|----|---------|----------------------|--------------------|-----|
| 555 | Conformal-prediction UQ | R/Python | ✅ | `pip install mapie` / `pip install crepes`; R `install.packages("probably")` | **Whole-library UQ blank.** Model-agnostic statistically-valid coverage for diagnostic/prognostic signatures (low-cost 二区 add). Needs separate calibration set. | prediction-set-size dot, coverage-vs-target calibration scatter, per-sample interval dumbbell | Front Bioinform 2025 (conformal genomic ML) |

---

## Explicitly NOT built (and why) — so they are not silently lost

- **Julia-only**: SGCRNA (genuine 2026-H1 method, but no R/Python package → cannot be a module here; the scale-free-free idea is reachable via SGCP which is already folded into 054).
- **GPU/DL-heavy** (belong to `bioinfo-DL-library/02 §E`, not turnkey-local): CELLama, STARS, Cellist, CellNEST, Nicheformer/SToFM, Uni-Dock/UniDock-Pro, Gnina, DiffDock-L/SigmaDock, Boltz-2/Boltzina, Deep Docking, BindFlow/OpenFE/FEP-ABFE, PLUMED/OPES.
- **No installable package / paper-only**: FusioMR, TLMR, LMM-MEC, life-course-pathway MR, common+rare-variant MR, InferPloidy, SpatialDG, st-Xprop, SPHENIC, Spa3D, STCS, SpaNiche, TGCN, InterVelo/Cell2fate/BayVel/etc preprints.
- **Web databases / GUI / web tools (not code packages)**: TCMSP, BATMAN-TCM 2.0, HERB 2.0, SwissTargetPrediction, STRING/STITCH, SymMap/TCMBank/NPASS, Open Targets/GeneCards/OMIM/DisGeNET (druggability already in 493), KEGG/Reactome/WikiPathways, SRplot, Cytoscape+cytoHubba+MCODE, Discovery Studio/LigPlot+/PyMOL, PLIP (subsumed by ProLIF).
- **Benchmark / review / guidance papers (cite-only, no package)**: all the FM/velocity/integration/CNV/registration benchmarks, Nat Methods 2025 perturbation critique, Stop-Chasing-C-index, Nat Methods 2024 leakage guide, SHAP critiques, WFCMS guideline.
- **Already covered** (reuse-not-rewrite): CellChat v2→051, CellRank2→072, COMMOT→073, NicheNet/MultiNicheNet→077, scVI/scANVI/Symphony/Azimuth→506, Geneformer/FMs→507, decoupler→076, Slingshot/Palantir→082/087, scMetabolism→510, sc-TWMR/FUSION→036-042, coloc→075, two-step/SEM mediation→508/499, IOBR→492, Mime/combo→496/059, SHAP(classification)→052, LASSO/SVM-RFE/RF→012/013/014, DCA/calibration/nomogram→016, LODO/meta→503, RCTD/NMF/MISTy→505, BayesPrism→520, SpatialGlue→521, oncoPredict/beyondcell→518, PyWGCNA/CEMiTool/GWENA/MEGENA/SGCP→054, pySCENIC→081, engines(Vina/GROMACS/gmx_MMPBSA)→086, ΔG heatmap→022, UpSet/Sankey/chord/network figures→003-011/498/515/047/007.

## Build sequence (recommended)

⭐ = highest-value + CRAN/Bioc + easiest local verification → build first:
527 sccomp, 528 miloR, 529 muscat, 531 LIANA+, 533 MRcare, 534 MVMR-cML, 538 NetRep,
539 SmCCNet, 541 BANKSY, 542 SVG, 543 squidpy-stats, 546 emapplot, 547 ProLIF,
548 bio3d-FEL, 550 TabPFN, 553 riskRegression, 556 PoseBusters.
Then the rest. **Each new module must ship**: synthetic example data + fixed seed +
vector figures + README + an **honest-baseline contrast** (the repo's "two cleavers"
rule: e.g. TabPFN vs LASSO/XGBoost; in-fold FS to avoid leakage) + a CATALOG line +
category-README line + a SERVER_DEPENDENCIES entry if 🟡 + clean `qc_lint`.

*Origin: 2026-06-26 user task — fold the non-DL seven-line 2026 survey into the reusable
code library "if it has a code tool". Gap analysis vs CATALOG done with a 7-agent
workflow; install commands web-verified that day.*
