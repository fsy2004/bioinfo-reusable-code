# 084 · NMF + ConsensusClusterPlus molecular subtyping

Unsupervised molecular subtyping from a feature-by-sample matrix using NMF rank selection, a consensus matrix, and a subtype heatmap.

| | |
|---|---|
| **Language / main dependencies** | R · `NMF` `ConsensusClusterPlus` `ComplexHeatmap` |
| **Purpose** | Partition samples into stable molecular subtypes without supervision |
| **Input** | `example_data/feature_matrix.csv` |
| **Output** | Subtype table and figures in `results/`; display figures in `assets/` |

## Input

CSV with feature names (genes / immune scores / pathway scores) in the first column and samples in the remaining columns. Values are non-negative.

## Method

`NMF` performs rank selection over k=kmin..kmax (cophenetic correlation) and extracts subtypes at the optimal k. `ConsensusClusterPlus` runs resampling-based consensus clustering and outputs a consensus matrix as robustness evidence, followed by a feature heatmap annotated by subtype.

Method citations: Brunet *et al.*, *PNAS* 2004 (NMF subtyping); Wilkerson & Hayes, *Bioinformatics* 2010 (ConsensusClusterPlus).

## Use cases

Tumor / disease molecular subtyping, immune co-infiltration subtyping, and pathway-activity subtyping, to identify patient subgroups with distinct molecular features.

## Features

- Runs directly on a matrix. `--k` sets the number of subtypes manually (default selects k by cophenetic correlation; manual confirmation against the rank-selection and consensus plots is recommended).
- Figures: NMF rank-selection curve, consensus matrix heatmap, and subtype-annotated feature heatmap.

## Outputs

| File | Type | Description |
|------|------|------|
| `assets/Consensus_matrix.png` | Consensus matrix | Subtyping robustness (diagonal blocks = subtypes) |
| `assets/NMF_rank_survey.png` | Curve | cophenetic vs k |
| `assets/Subtype_heatmap.png` | Heatmap | Feature pattern per subtype |

![consensus](assets/Consensus_matrix.png)

## Usage

```bash
Rscript 084_NMF_consensus_subtyping.R                              # 示例(自动选 k)
Rscript 084_NMF_consensus_subtyping.R --input data/mat.csv --k 3   # 指定 3 亚型
```

## Dependencies

```r
install.packages("NMF"); BiocManager::install(c("ConsensusClusterPlus","ComplexHeatmap"))
```
