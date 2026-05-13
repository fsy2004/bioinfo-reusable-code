# 14_单细胞虚拟扰动_扰动数据库

本模块用于把现有单细胞流程扩展到 perturb-seq 数据、GRN 虚拟扰动和多基因组合扰动预测。

建议输入：
- Seurat RDS 或 AnnData h5ad
- perturbation/control 分组列
- 候选基因列表

建议输出：
- 扰动距离和 E-test 结果
- GEARS 预测表达矩阵
- CellOracle GRN 扰动向量和调控靶基因

