# Server dependency inventory — full-version provisioning

Many turnkey modules run locally with **synthetic data + a re-implemented method core
or a simple baseline**, so they work without heavy/GPU packages. For real analyses on
the **server**, install the full package listed here. Status legend:

| Tag | Meaning |
|-----|---------|
| ✅ **real, local** | uses the real package; installed on this machine; runs locally |
| 🟡 **core/baseline local** | turnkey re-implements the method core or ships an honest baseline; runs locally; **install the real package on the server for the full version** |
| 🔴 **needs install / GPU** | requires a heavy or GPU package not installed here; script is ready, marked needs-server |

Audited 2026-06-25 by direct `find_spec` / `installed.packages()` probe.
Install proxy/mirror tips: see [[user_github_cli]] — GitHub via `https://gh-proxy.org/`,
pip via `-i https://pypi.tuna.tsinghua.edu.cn/simple` (international PyPI drops mid-download
through the proxy), `git push` via `http://127.0.0.1:7892`.

---

## This-session modules (506-521)

| # | Module | Real package | Status | Turnkey behavior | Server install |
|---|--------|--------------|:------:|------------------|----------------|
| 506 | scVI/scANVI integration | `scvi-tools` | ✅ | real model, trains on CPU | `pip install scvi-tools` (GPU optional, faster) |
| 507 | Geneformer in-silico | `geneformer`,`transformers`,`datasets` | 🔴 | PCA/Leiden baseline runs; FM path marked needs-GPU | `pip install transformers datasets`; clone Geneformer + download pretrained model (GPU) |
| 508 | Two-step mediation MR | base R | ✅ | real (closed-form IVW + Monte-Carlo) | — |
| 509 | Communication functional loop | `nichenetr` | 🟡 | NicheNet core re-implemented (synthetic prior) | `remotes::install_github("saeyslab/nichenetr")` + real `ligand_target_matrix` |
| 510 | scMetabolism pathway activity | `scMetabolism` | 🟡 | AUCell/UCell-style scoring core | `remotes::install_github("wu-yc/scMetabolism")` (or `UCell` + KEGG sets) |
| 511 | TF convergence (regulon×JASPAR×DepMap) | `RcisTarget`/pySCENIC, JASPAR, DepMap | 🟡 | 3-evidence convergence on synthetic inputs | pySCENIC (mod 081) + `JASPAR2024`/`TFBSTools` + DepMap CRISPR table |
| 512-516 | Advanced figures (raincloud/ridge/dumbbell/chord/composite) | base ggplot/matplotlib | ✅ | real | — |
| 517 | VECTOR differentiation direction | `Vector` (jumphone) | 🟡 | VECTOR core re-implemented (expression breadth → vector field) | `remotes::install_github("jumphone/Vector")`; potency from CytoTRACE2 (mod 082) |
| 518 | beyondcell drug response | `beyondcell` | 🟡 | BCS core re-implemented (UCell up−down) | `BiocManager::install("beyondcell")` + PSc/SSc signatures |
| 519 | Local MR pipeline (zero API) | `TwoSampleMR`,`ieugwasr`,`MRPRESSO`,`plinkbinr` | 🟡 | runs; **local LD clumping skipped without a 1000G `--bfile`** | `MRPRESSO`, `ieugwasr`, `plinkbinr` + 1000G EUR plink reference |
| 520 | BayesPrism deconvolution | `BayesPrism` | ✅ | real, runs locally | `remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")` |
| 521 | SpatialGlue multi-omics | `SpatialGlue`,`torch-geometric` | 🔴 | concat-PCA baseline runs; GNN needs GPU | `pip install SpatialGlue torch torch-geometric` (GPU) |

---

## Pre-existing heavy modules (install on server as needed)

| # | Module | Real package | Status | Server install |
|---|--------|--------------|:------:|----------------|
| 047 | RcisTarget TF motif network | `RcisTarget` | 🔴 | `BiocManager::install("RcisTarget")` + motif rankings |
| 068 | GEARS combo perturbation | `gears` | 🔴 | `pip install cell-gears` (GPU) |
| 069 | CellOracle GRN perturbation | `celloracle` | 🔴 | `pip install celloracle` |
| 070 | chemCPA | `chemCPA` | 🔴 | clone chemCPA repo (GPU) |
| 071 | scDrug | `scDrug` deps | 🔴 | per scDrug repo |
| 072 | CellRank fate drivers | `cellrank`,`scvelo` | 🔴 | `pip install cellrank scvelo` |
| 073 | COMMOT spatial communication | `commot` | 🔴 | `pip install commot` |
| 074 | Tangram sc→spatial | `tangram-sc` | 🔴 | `pip install tangram-sc` (GPU optional) |
| 076 | decoupler TF/pathway activity | `decoupler` | ✅ | installed |
| 077 | NicheNet ligand-target | `nichenetr` | 🔴 | `remotes::install_github("saeyslab/nichenetr")` |
| 080 | cell2location niche | `cell2location`,`squidpy` | 🔴 | `pip install cell2location squidpy` (GPU) |
| 081 | pySCENIC regulon | `pyscenic` | 🔴 | `pip install pyscenic` + cisTarget DBs |
| 082 | Trajectory (Palantir/Slingshot/tradeSeq/CytoTRACE2) | `palantir`; `slingshot`,`tradeSeq`,`CytoTRACE2` | 🔴 | `pip install palantir`; `BiocManager::install(c("slingshot","tradeSeq"))`; CytoTRACE2 from GitHub |
| 083 | MOFA/DIABLO multi-omics | `MOFA2`,`mixOmics` | 🔴 | `BiocManager::install(c("MOFA2","mixOmics"))` |
| 085 | Squidiff diffusion perturbation | `squidiff` | 🔴 | per repo (GPU) |
| 494 | GenKI VGAE-KO | `GenKI` | 🔴 | `pip install GenKI` (GPU) |

Installed and ready locally (no action): `Seurat`, `CellChat`, `monocle3`, `WGCNA`,
`hdWGCNA`, `Scissor`, `UCell`, `ConsensusClusterPlus`, `NMF`, `TwoSampleMR`,
`MendelianRandomization`, `coloc`, `Boruta`, `ComplexHeatmap`, `circlize`, `ggalluvial`,
`survival`, `timeROC`, `rms`, `glmnet`, `igraph`; Python `torch`, `scvi-tools`, `lightning`,
`scanpy`, `anndata`, `decoupler`, `pyro`.

---

## Server provisioning quick-start

```bash
# Python (GPU box) — use a domestic mirror if pip stalls: -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scvi-tools transformers datasets cell2location squidpy cellrank scvelo \
            commot tangram-sc palantir celloracle pyscenic cell-gears \
            SpatialGlue torch torch-geometric
```
```r
# R — via gh-proxy if GitHub is blocked: remotes::install_github("https://gh-proxy.org/https://github.com/<u>/<repo>")
BiocManager::install(c("RcisTarget","beyondcell","slingshot","tradeSeq","MOFA2","mixOmics","UCell"))
remotes::install_github(c("saeyslab/nichenetr","wu-yc/scMetabolism","jumphone/Vector",
                          "Danko-Lab/BayesPrism/BayesPrism","cytotrace2/CytoTRACE2"))
install.packages(c("MRPRESSO","ieugwasr","plinkbinr"))
```

GPU-heavy (Geneformer, SpatialGlue, GEARS, CellOracle, cell2location, chemCPA, Squidiff,
GenKI) belong on the AutoDL GPU box; the local turnkey versions are baselines/cores for
development and figure design only.
