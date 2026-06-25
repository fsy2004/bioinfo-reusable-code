# 035 · Multi-method combination intersection selection

Enumerate all k-method combinations of feature lists from multiple ML methods, rank them by intersection size, and report the best combination with an UpSet plot.

| | |
|---|---|
| **Language / main dependencies** | R · `theme_pub` + `UpSetR` |
| **Purpose** | Find the method combination with the largest intersection and its consensus features |
| **Input** | `example_data/method_sets/` (6 method feature lists) |
| **Output** | `results/` combination table, intersection, and figures · display figures in `assets/` |

## Input

The input is a directory (`--input`) containing several method feature lists (csv, with column `variable` or the gene in the first column). Files named `importanceGene.<method>.csv` or anything else; the method name is taken automatically.

## Method

Enumerate all `--pick` (default 5) method combinations, compute the intersection size of each combination, and sort. By default the combination with the largest intersection is selected (or specified by `--methods`), and its consensus features plus an UpSet plot are output.

## Use case

When there are many methods, objectively select a method combination that is both robust and retains enough features, avoiding subjective selection.

## Features

- Combination enumeration and ranking in one command; `--pick/--methods` configurable.
- Top combination intersection ranking (viridis) and UpSet plot for the selected combination.

## Outputs

| File | Figure type | Description |
|------|------|------|
| `assets/Combo_ranking.png` | Ranking | Intersection size of top combinations |
| `assets/Combo_UpSet.png` | UpSet | Feature intersection of the selected combination |
| `results/all_combinations.csv` · `selected_combo_intersection.txt` | Table | All combinations / selected intersection |

![ranking](assets/Combo_ranking.png)

## Usage

```bash
Rscript 035_method_combo_intersection.R                                   # 示例
Rscript 035_method_combo_intersection.R --input data/method_sets --pick 5 --methods "RF,Lasso,SVM,GBM,PLS"
```

## Dependencies

```r
install.packages("UpSetR")
```
