# 09 · Mendelian Randomization / GWAS Processing

GWAS instrument preparation, MR causal inference, and sensitivity, directionality, and mediation analysis. [032](032_mr_twosamplemr/) is the primary module, providing a self-contained MR implementation that does not depend on TwoSampleMR.

## Primary module

| Module | Purpose | Language | Output figures |
|------|------|------|--------|
| [032 MR analysis](032_mr_twosamplemr/) | IVW/Egger/WM + sensitivity | R | scatter, forest, funnel, leave-one-out |

## GWAS instrument processing (upstream helpers)

| Module | Function |
|------|------|
| 028 VCF significant SNP filtering, 029 LD clumping, 030 add EAF, 031 weak instrument (F) filtering | Produce harmonized instruments from raw GWAS (input to 032) |

## Advanced MR variants (kept for reference)

| Module | Method | Dependency |
|------|------|------|
| 033 MR fallback template | Basic MR | Same as 032 |
| 043 MR + directionality, 055 immune-cell bidirectional MR | Steiger directionality | TwoSampleMR |
| 075 MR + coloc causal evidence chain | colocalization | coloc, LocusZoom |
| 079 pQTL MVMR protein mediation | Multivariable MR | MVMR |
| 499 lavaan SEM mediation paths | Structural equation modeling | lavaan |

Module 032 follows the [unified framework conventions](../_framework/CONVENTIONS.md); its core IVW/Egger/weighted-median methods are a self-contained implementation and can be applied directly to the data in 033/043/055. The advanced variants (coloc/MVMR/SEM) depend on dedicated packages, and the original scripts are kept for reference.
