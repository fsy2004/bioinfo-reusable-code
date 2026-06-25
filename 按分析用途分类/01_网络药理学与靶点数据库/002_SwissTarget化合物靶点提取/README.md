# 002 · SwissTargetPrediction compound target extraction

Extract a deduplicated target list from a SwissTargetPrediction export, with optional filtering by prediction probability.

| | |
|---|---|
| Language / main dependency | R · base |
| Input | `example_data/SwissTargetPrediction_export.csv` |
| Output | `results/targets.csv` |

## Input

SwissTargetPrediction export CSV. The `Gene` column and the `Probability` score column are detected automatically.

## Method

Extract the gene column, optionally filter by `Probability >= --score-min` (typically 0.1), then deduplicate.

## Usage

```bash
Rscript 002_extract_targets.R --score-min 0.1
```

## Outputs

No figures. `results/targets.csv`.

## Notes

The probability threshold is adjustable. Predicted compound targets can be merged with other sources such as CTD (see 003). Uses the same engine as 001 and 004.

## Dependencies

R, base.
