# 01 · Network pharmacology and target databases

Extract targets from databases, intersect across multiple sources, and identify core action targets.

| Module | Purpose | Language | Output figures |
|------|------|------|--------|
| [001 CTD target extraction](001_CTD化合物靶点提取/) | CTD export to target list | R | none (produces list) |
| [002 SwissTarget extraction](002_SwissTarget化合物靶点提取/) | Predicted targets (probability filter) | R | none (produces list) |
| [004 GeneCards extraction](004_GeneCards疾病靶点提取/) | Disease targets (score filter) | R | none (produces list) |
| [003 CTD∩Swiss Venn](003_CTD_Swiss靶点并集Venn/) | Compound target intersection/union | R | Venn, bar |
| [005 OMIM∩GeneCards Venn](005_OMIM_GeneCards疾病靶点Venn/) | Disease target intersection/union | R | Venn, bar |
| [006 Disease×compound Venn](006_疾病化合物靶点交集Venn/) | Core action target intersection | R | Venn, bar |
| [011 DEG×drug target](011_差异基因_药物靶点交集/) | Multi-set intersection | R | Venn, UpSet, bar |
| 493 Druggability scoring | OpenTargets/DGIdb/ChEMBL scoring | Python | none (live API table) |

## Typical pipeline
```
001/002 compound targets ─┐
                    ├─▶ 003 merge ─┐
004 disease targets ───────┘             ├─▶ 006 disease×compound intersection = core targets ─▶ 007 enrichment / PPI
005 disease targets merge ─────────────────┘
```

003/005/006/011 share the same target-intersection engine (`venn_pub`, dependency-free Venn plus UpSetR). 001/002/004 share the extraction engine.

493: Python, depends on three live APIs (OpenTargets/DGIdb/ChEMBL) plus mygene, outputs a druggability scoring table (no figures), kept as a reference script (requires network access; the script header notes that installation must be confirmed manually).

All modules follow the [shared framework conventions](../_framework/CONVENTIONS.md).
