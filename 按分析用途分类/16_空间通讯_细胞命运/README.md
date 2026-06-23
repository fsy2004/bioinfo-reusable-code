# 16_空间通讯_细胞命运

> ⏭️ **状态**:本类别为空间通讯/细胞命运方法,多依赖 Python 空间组学栈(CellRank/COMMOT/Tangram/cell2location/Squidpy/decoupler)或 NicheNet 大型先验网络,**本地未渲染示例图**,保留原脚本作参考(各脚本头部含依赖与可复现命令)。配图规范见 [统一框架](../_framework/CONVENTIONS.md)。

---

本模块用于空间通讯、单细胞到空间映射、细胞命运推断、空间生态位和 TF/通路活性评分。

## 脚本

| 脚本 | 作用 |
|---|---|
| `072_CellRank_命运概率与驱动基因.py` | 基于 scVelo velocity graph 运行 CellRank，输出终末状态、命运概率和 driver genes。 |
| `073_COMMOT_空间细胞通讯.py` | 基于空间坐标和 ligand-receptor 数据库计算空间通讯 sender/receiver score。 |
| `074_Tangram_单细胞到空间映射.py` | 把 scRNA-seq 细胞映射到空间转录组 spot。 |
| `076_decoupler_TF通路活性评分.py` | 用 CollecTRI 或 PROGENy 推断 TF/通路活性。 |
| `077_NicheNet_配体靶基因通信推断.R` | 从 receiver DE genes 和 LR/target prior 推断 ligand activity 与 ligand-target 链条。 |
| `080_cell2location_Squidpy_空间生态位.py` | 整合 cell abundance 表并运行 Squidpy 空间邻域富集。 |

## 推荐输入

- AnnData h5ad 或 Seurat 转换结果。
- 空间坐标、组织切片信息。
- ligand-receptor 数据库、receiver 差异基因、背景表达基因。
- 细胞类型注释或 cell abundance 矩阵。

## 推荐输出

- CellRank fate probabilities 和 driver genes。
- COMMOT 空间通讯网络。
- Tangram 单细胞到空间映射结果。
- NicheNet ligand activity 和 ligand-target links。
- Squidpy neighborhood enrichment。
