# 02_NHANES

NHANES complex-survey analysis (survey-weighted descriptives, weighted
regression, weighted prevalence).

| Module | Purpose | Language |
|--------|---------|----------|
| [528 NHANES survey-weighted](528_nhanes_survey_weighted/) | svydesign → svymean/svyby/svyglm; weighted vs unweighted | R |

Turnkey on synthetic NHANES-shape data; survey calls grounded in the real
`../99_external_sources/nhanes/` vignettes. See the module README for the design
pitfalls (build design before subsetting; weights + strata + PSU).
