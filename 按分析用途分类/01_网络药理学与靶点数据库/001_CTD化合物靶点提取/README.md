# 001 · CTD compound target extraction

Extracts a deduplicated target gene list from a CTD database export file for downstream Venn or enrichment analysis.

| | |
|---|---|
| **Language / main dependency** | R · base |
| **Input** | `example_data/CTD_export.csv` |
| **Output** | `results/targets.csv` |

## Input

CTD export CSV. The `Gene Symbol` column (and the `Reference Count` score column) are detected automatically.

## Method

Read the export table, extract the gene column, optionally filter by score, deduplicate, and produce the target list.

## Usage

Network pharmacology first step: organize CTD compound-gene associations into a standard target list.

```bash
Rscript 001_extract_targets.R                                  # example
Rscript 001_extract_targets.R --input data/CTD_export.csv
```

## Outputs

No figures. `results/targets.csv` (deduplicated targets).

## Dependencies

R base only. The gene column and score column are detected automatically.
