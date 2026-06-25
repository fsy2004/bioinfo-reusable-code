# 008 · GEO expression matrix tidying (probe to gene)

Maps a GSE series matrix and GPL platform annotation to a gene-level expression matrix `geneMatrix.csv` for downstream use in modules 009/010.

| | |
|---|---|
| **Language / main dependency** | R · base (no third-party packages) |
| **Purpose** | Map and collapse probe-level GEO data into a gene-level matrix |
| **Input** | `example_data/` (GSE series matrix + GPL platform txt) |
| **Output** | `results/geneMatrix.csv` |

## Input

| File | Format | Description |
|------|--------|-------------|
| `GSE*_series_matrix.txt` | tab txt | Probe × sample expression table with an `ID_REF` header row (GEO standard download format) |
| `GPL*.txt` | tab txt | Platform annotation; one column holds the gene Symbol. The column index is set by `--symcol` (1-based, example = 2) |

The directory automatically detects `GSE*`/`GPL*` files. The Symbol column differs across platforms, so verify `--symcol`.

## Method

Locate the `ID_REF` start row and read the expression table, parse the GPL to build a probe-to-Symbol mapping (take the left side of `///`, discard non-words containing spaces), align with `merge`, then collapse multiple probes for the same gene with `aggregate` using the mean. The result is a gene-level matrix.

## Usage

This is the first preprocessing step for GEO microarray analysis. It produces a standard gene matrix that feeds module 009 (grouping/normalization) and then module 010 (differential analysis).

```bash
Rscript 008_GEO_expr_matrix_tidy.R                                  # example
Rscript 008_GEO_expr_matrix_tidy.R --gse GSExxx_series_matrix.txt --gpl GPLxxx.txt --symcol 11
```

## Outputs

No figures. `results/geneMatrix.csv`: first column `geneSymbol` followed by per-sample expression values.

## Dependencies

No additional installation required (base R).
