## Cleanup changelog (resolved 2026-06-26)

All previously-flagged duplicates and numbering inconsistencies have been actioned:

| Item | Action taken |
|------|--------------|
| `06/019_immune_infiltration_source_generic.R` | **Deleted** — was byte-identical to 017 (only the header `# 编号` line differed) |
| `06/020_immune_infiltration_scoring_generic.R` | **Deleted** — was byte-identical to 018 |
| `06/021_immune_infiltration_viz.R` (loose) | **Deleted** — stale pre-modularization copy of the `021_immune_infiltration_viz/` directory module |
| `08/082_palantir_branch_probability.py` | **Renumbered → 087** (resolves the 082 collision with the Slingshot/tradeSeq/CytoTRACE2 module) |
| `09/497_lavaan_sem_mediation_path.R` | **Renumbered → 499** (resolves the cross-category 497 collision with 12's scSurvival) |
| `18/14_ai_scientific_figures/` | **Renamed → `ai_scientific_figures/`** (dropped the stale `14_` category-number prefix) |
| `20/*` five templates | **Renumbered → 522–526** to match the repo-wide `NNN_` convention |
| `04/045_multimodalad_ml_models.R` | **Marked TEMPLATE** in its header — it `source()`s a project-specific `refer.ML.R` helper that is not bundled, so it is a reference, not a turnkey run. For a turnkey multi-method run use 034; for prognostic combos use 059/496 |

> 045, 059, 496 are heavy/upstream-derived scripts that sit in category 04 but are
> really modelling/prognostic work — kept here for provenance; see status marks.
> New module numbers continue at **561+**.
