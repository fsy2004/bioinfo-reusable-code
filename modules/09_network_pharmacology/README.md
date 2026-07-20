# 09 · 网络药理学 — Network pharmacology

本域共 8 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 靶点数据库提取 — Target databases

- [001_ctd_compound_targets](01_target_databases/001_ctd_compound_targets) — Extract & de-dup compound targets from a CTD export
- [002_swisstarget_compound_targets](01_target_databases/002_swisstarget_compound_targets) — Extract compound targets from a SwissTargetPrediction export
- [004_genecards_disease_targets](01_target_databases/004_genecards_disease_targets) — Extract disease targets from a GeneCards export

## 靶点交集与集合图 — Target intersection

- [003_ctd_swiss_target_union_venn](02_target_intersection/003_ctd_swiss_target_union_venn) — Union/intersection of CTD vs Swiss compound targets
- [005_omim_genecards_target_venn](02_target_intersection/005_omim_genecards_target_venn) — Union/intersection of OMIM vs GeneCards disease targets
- [006_disease_compound_target_venn](02_target_intersection/006_disease_compound_target_venn) — Disease ∩ compound targets → core targets
- [011_deg_drug_target_intersection](02_target_intersection/011_deg_drug_target_intersection) — Multi-set DEG ∩ drug ∩ disease target intersection

## 成药性评分 — Druggability

- [493_opentargets_dgidb_chembl_druggability.py](03_druggability/493_opentargets_dgidb_chembl_druggability.py) — Composite druggability score for a target set (live APIs)
