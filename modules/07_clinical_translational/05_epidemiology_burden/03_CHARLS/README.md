# 03_CHARLS

CHARLS-style longitudinal cohort analysis (wave description, cross-wave equating,
longitudinal trend, mixed model, incident-event survival).

| Module | Purpose | Language |
|--------|---------|----------|
| [529 CHARLS longitudinal cohort](529_charls_longitudinal_cohort/) | Table 1 + equipercentile equating + LMM + Cox/KM | R |

Turnkey on a synthetic multi-wave panel; grounded in
`../99_external_sources/charls_memory_equating/`. The `equate` package is
reimplemented from weighted-ECDF inversion (not installed); see the module README.
