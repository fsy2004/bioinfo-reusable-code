# 009 · GEO sample grouping and normalization

Normalizes an expression matrix and appends sample-group labels, producing `Sample_Type_Matrix.csv` for use by module 010.

| | |
|---|---|
| **Language / main dependency** | R · `limma` (optional `readxl`) |
| **Purpose** | Normalize the expression matrix and write sample-group labels |
| **Input** | `example_data/geneMatrix.csv` + `sample_group.csv` |
| **Output** | `results/Sample_Type_Matrix.csv` |

## Input

| File | Format | Required columns | Notes |
|------|------|------|------|
| `--expr` | csv | First column gene names + sample columns | Gene-level matrix produced by module 008 |
| `--group` | csv/xlsx | Column 1 sample name, column 2 type | Sample names must match the expr column names; types such as `con`/`tre` |

## Method

`limma::avereps` (average duplicate genes), quantile check to decide whether normalization is needed, `log2(x+1)`, `normalizeBetweenArrays` (between-array normalization), filter/sort samples by the group table, then append a `_type` suffix to each sample name.

## Purpose

Links modules 008 and 010: normalizes the gene matrix and applies group labels so that downstream differential analysis (010) can identify Control/Disease from the suffix.

## Notes

- Runs from `--expr`/`--group`; the group table accepts csv or xlsx.
- Automatic log2 decision: log transformation is applied based on the expression value distribution, avoiding double log.
- Output file name and format match the input expected by module 010.

## Outputs

No figures. `results/Sample_Type_Matrix.csv` (sample names carry a `_type` suffix) and `Sample_Summary.txt` (sample count per group).

## Usage

```bash
Rscript 009_GEO_sample_grouping.R                                       # 示例
Rscript 009_GEO_sample_grouping.R --expr geneMatrix.csv --group sample_group.csv
```

## Dependencies

```r
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install("limma")
# 若分组表为 xlsx: install.packages("readxl")
```
