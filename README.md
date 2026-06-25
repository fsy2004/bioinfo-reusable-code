# Reusable Bioinformatics Code Library

A set of self-contained R and Python modules for common bioinformatics analyses.
Each module ships with a small example dataset and runs from the command line to
produce vector, journal-style figures. Replace the example with your own data to
reuse it in a project.

- ~28 runnable modules across 21 analysis categories
- A shared plotting framework (`_framework/`) for consistent figure style
- A project scaffold, quality checklist, and static linter for reproducible pipelines
- Tested with R 4.4 and Python 3.12

## Repository layout

```
按分析用途分类/              analysis modules, grouped by purpose
├── _framework/             shared toolkit (themes, palettes, scaffold, linter)
├── 01_..21_..              module categories (catalog below)
└── <NNN_module>/           one folder per module:
    ├── <NNN_module>.R|.py   main script (runs on example_data/ by default)
    ├── README.md            input spec, method, outputs
    ├── example_data/        small synthetic input
    └── assets/              committed preview figures
```

Run-time outputs (`results/`, `figures/`) are git-ignored.

## Quick start

```bash
git clone https://github.com/fsy2004/bioinfo-reusable-code.git
cd bioinfo-reusable-code/按分析用途分类

# run a module on its bundled example data
Rscript 03_GEO转录组整理与差异分析/010_GEO差异分析_火山热图PCA/010_*.R

# run on your own data
Rscript 03_GEO转录组整理与差异分析/010_GEO差异分析_火山热图PCA/010_*.R \
        --input your_matrix.csv --outdir results/run1
```

Each module folder documents its exact input format, method, and outputs.

## Example outputs

Rendered directly from the bundled example data:

| Differential expression | Single-cell clustering | Mendelian randomization |
|:---:|:---:|:---:|
| ![volcano](按分析用途分类/03_GEO转录组整理与差异分析/010_GEO差异分析_火山热图PCA/assets/DEG_volcano.png) | ![umap](按分析用途分类/08_单细胞_空间转录组_细胞轨迹/046_单细胞发表级图/assets/UMAP_clusters.png) | ![mr](按分析用途分类/09_孟德尔随机化_GWAS处理/032_MR_TwoSampleMR分析/assets/MR_scatter.png) |

## Module catalog

| #  | Category | Modules | Typical output |
|----|----------|---------|----------------|
| 01 | Network pharmacology & target databases | 001–006, 011 | Venn, UpSet, target tables |
| 02 | GO / KEGG enrichment | 007 | dot/bar plots, pathway graph |
| 03 | Transcriptomics (GEO) & differential expression | 008–010, 056 | volcano, heatmap, PCA, batch correction |
| 04 | Machine-learning feature selection | 012–015, 034, 035, 052, 059 | LASSO, RF, SVM-RFE, SHAP, AUC heatmap |
| 05 | Diagnostic models & validation | 016, 063 | ROC, calibration, DCA, nomogram |
| 06 | Immune infiltration | 017–021, 492 | composition, boxplot, correlation |
| 07 | Molecular docking & dynamics | 022, 086 | binding-energy bubble, MD metrics |
| 08 | Single-cell / spatial / trajectory | 023–027, 044–051, 058, 062, 082 | UMAP, dot plot, marker heatmap |
| 09 | Mendelian randomization & GWAS | 028–033, 043, 055, 075, 079 | MR scatter, forest, funnel, leave-one-out |
| 10 | TWAS (single-cell eQTL weights) | 036–042 | weight tables |
| 11 | WGCNA co-expression | 054 | soft-threshold, module-trait heatmap |
| 12 | TCGA prognosis (reference only) | 048, 057, 060 | KM, time-dependent ROC, risk plot |
| 13 | Transcription-factor regulation / circos | 047, 053, 081 | chromosome circos, regulon network |
| 14 | Single-cell in-silico perturbation | 067–069, 085, 494, 495 | gene-knockout effects |
| 15 | Drug perturbation / repurposing | 070, 071, 078 | pharmacovigilance signals |
| 16 | Spatial communication / cell fate | 072–074, 076, 077, 080 | CellRank, niche maps |
| 17 | Advanced result figures | 498 | alluvial / Sankey |
| 18 | External method sources | manifest only | — |
| 19 | Multi-omics integration & subtyping | 083, 084 | MOFA, consensus clustering |
| 20 | Mutation / CNV / methylation / proteome / metabolome | 5 templates | oncoprint, volcano, heatmap |
| 21 | Disease burden (GBD / NHANES / CHARLS) | external / spec | — |

Categories 10, 14, 16 and parts of 07/12 require heavy or GPU-bound toolchains
(FUSION, GROMACS, deep-learning models); their scripts and dependency notes are
kept for reference rather than local one-command rendering.

## Framework (`_framework/`)

Shared by all modules so figures and I/O stay consistent:

- `theme_pub.R` / `pubstyle.py` — journal theme, discrete palettes (NPG/AAAS/Lancet/…),
  viridis for continuous scales, and `save_fig()` (vector PDF + 300 dpi PNG)
- `CONVENTIONS.md` — module layout, run conventions, figure rules
- `ANALYSIS_TEMPLATE/` — scaffold for a new multi-step project: central config
  (seed, relative paths, parameters), setup with checkpointed steps
  (`cache_step`), logged statistics, and an environment snapshot; R and Python versions
- `QUALITY_CHECKLIST.md` — pre-/in-/post-analysis checklist
- `qc_lint.py` — static checks for hard-coded paths, missing random seeds,
  non-vector figure exports, and missing environment snapshots

## Conventions

- Modules run on bundled example data with no edits; use `--input` / `--outdir` to switch.
- No absolute paths or `setwd()`; figures are exported as vector PDF + 300 dpi PNG.
- Reuse the framework instead of re-implementing themes or I/O.
- Figure text in English; analysis logic left intact when standardizing a module.

## License

Each module follows the license of the tools and methods it uses. Vendored
third-party code (e.g. category 18) keeps its original license — see the relevant
module README and upstream repository.
