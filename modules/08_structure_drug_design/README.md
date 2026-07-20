# 08 · 结构与药物设计 — Structure & drug design

本域共 6 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 分子对接 — Molecular docking

- [022_docking_binding_energy_viz](01_docking/022_docking_binding_energy_viz) — Binding-energy heatmap + strongest-binding ranking
- [547_prolif_interaction_fingerprint](01_docking/547_prolif_interaction_fingerprint) — ProLIF protein-ligand interaction fingerprint + residue occupancy
- [556_posebusters_validity_panel](01_docking/556_posebusters_validity_panel) — PoseBusters physical-validity check panel for docking/AI poses

## 分子动力学 — Molecular dynamics

- [086_vina_gromacs_mmpbsa_mdanalysis_pipeline.py](02_md_simulation/086_vina_gromacs_mmpbsa_mdanalysis_pipeline.py) — Vina docking + GROMACS MD + MM-PBSA pipeline
- [548_bio3d_md_dccm_pca](02_md_simulation/548_bio3d_md_dccm_pca) — bio3d ensemble/MD: PCA + DCCM + RMSF with collectivity null

## 虚拟筛选与打分 — Virtual screening & scoring

- [596_scorch2_virtual_screening](03_virtual_screening/596_scorch2_virtual_screening) — SCORCH2 双视图共识 ML 重打分的守卫式封装 + 本机可跑的虚拟筛选富集评测骨架(EF1%/EF5%/BEDROC/AUROC,按靶点分层、GroupKFold 防泄漏)
