# 15 · Drug Perturbation / Drug Repositioning

Links single-cell states, candidate targets, and real-world drug signals to drug response prediction, repositioning, and pharmacovigilance.

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [078 FAERS pharmacovigilance](078_FAERS药物警戒信号挖掘/) | ROR/PRR/BCPNN/EBGM signals | R | ROR forest plot, signal heatmap | Runnable |
| 070 chemCPA drug perturbation expression prediction | Deep-learning perturbation profile prediction | Python | Scatter, heatmap | Heavy (GPU) |
| 071 scDrug single-cell drug response | Cluster-level drug sensitivity | Python | UMAP, heatmap | Heavy |

**078**: produces pharmacovigilance figures from the report table (follows the [unified framework conventions](../_framework/CONVENTIONS.md)).
**070/071**: deep-learning / external models (chemCPA, scDrug) requiring GPU and large models; not rendered locally. The original scripts are kept for reference (see per-script header dependencies and reproducible commands).
