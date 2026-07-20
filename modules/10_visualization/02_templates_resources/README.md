# 10.02 — 模板与外部资源（清单式，不 vendored）

第三方完整源码树不提交进本仓库：它们会把 git 撑大（约 19k 文件）、且重复了上游已有的代码。
这里只记录来源；真正要用某个方法时，从上游 URL 或 pip / CRAN 取。

## 本地引用（git-ignored，不提交）

- `ai_scientific_figures/` — vendored AutoFigure-Edit（ICLR'26）：把方法描述转成可编辑的 SVG 示意图。
  本地留作参考，上游仓库与论文一并下载但不跟踪（`.gitignore`）。

## 已提交的内容

- `templates/` — 出图模板。
- `advanced_figure_tools.csv` · `download_advanced_figure_tools.ps1` · `literature_download_links_for_fdm.txt`
  — 高级图工具清单与批量下载脚本。

## 主要上游方法（按域索引到本库的可跑封装）

| 域 | 上游方法 | 本库封装位置 |
|---|---|---|
| 虚拟扰动 | GEARS, CellOracle, scPerturb, scTenifoldKnk, Geneformer, RegVelo | `03_virtual_perturbation/01_insilico_knockout` |
| 基因调控网络 | SCENIC / pySCENIC, RcisTarget, decoupler | `03_virtual_perturbation/02_grn_inference` |
| 药物扰动 | chemCPA, scDrug, beyondcell | `03_virtual_perturbation/04_drug_perturbation` |
| 轨迹 / 命运 | CellRank, scVelo, scTour, Slingshot, tradeSeq, Palantir, CytoTRACE2 | `01_single_cell/06_trajectory_velocity` |
| 空间 | cell2location, RCTD/spacexr, SPOTlight, BANKSY, nnSVG, Squidpy, COMMOT, Tangram, PASTE2 | `02_spatial_transcriptomics/*` |
| 细胞通讯 | CellChat, NicheNet, LIANA | `02_spatial_transcriptomics/05_cell_communication` |
| 多组学 / 分型 | MOFA2, mixOmics, NMF | `06_bulk_omics/04_multiomics_integration` |
| 因果推断 | TwoSampleMR, coloc, MVMR, MRBEE, SharePro | `04_causal_inference_genetics/*` |
| 对接 / 动力学 | AutoDock Vina, smina, GROMACS, gmx_MMPBSA, MDAnalysis, ProLIF | `08_structure_drug_design/*` |
| 药物警戒 | FAERS pvda | `07_clinical_translational/04_pharmacovigilance` |
