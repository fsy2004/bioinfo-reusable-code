# 004 · GeneCards Disease Target Extraction

Extract a deduplicated list of disease targets from a GeneCards export, optionally filtered by relevance score.

| | |
|---|---|
| **Language / Dependencies** | R · base |
| **Input** | `example_data/GeneCards_export.csv` |
| **Output** | `results/targets.csv` |

## Input

GeneCards export CSV. The `Gene Symbol` column and the `Relevance score` column are detected automatically.

## Method

Extract the gene column, optionally filter by `Relevance score >= --score-min` (typically 1 to 10), then deduplicate.

## Usage

```bash
Rscript 004_extract_targets.R --score-min 5
```

## Outputs

No figures. `results/targets.csv`.

The output can be merged with disease genes from other sources such as OMIM (module 005).

## Notes

The relevance score threshold is adjustable. This module shares the same engine as 001 and 002.
