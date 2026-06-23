# 14_单细胞虚拟扰动_扰动数据库

> ⏭️ **状态**:本类别为单细胞虚拟扰动方法,多依赖深度学习/图神经网络(GEARS/CellOracle/Squidiff/GenKI,需 GPU 与训练好的模型)或专用扰动数据库,**本地未渲染示例图**,保留原脚本作参考(各脚本头部含依赖与可复现命令)。067_scPerturb(R)与 495_bulkVGK(scTenifoldKnk)较轻,有合适数据时可运行。配图规范见 [统一框架](../_framework/CONVENTIONS.md)。

---

本模块用于把单细胞分析从“描述细胞状态”推进到“预测扰动后细胞会如何变化”。

## 脚本

| 脚本 | 作用 |
|---|---|
| `067_scPerturb_扰动数据Etest.R` | 对 perturbation/control 分组计算扰动距离和 E-test。 |
| `068_GEARS_单细胞组合扰动预测.py` | 调用 GEARS 预测单基因或多基因组合扰动表达响应。 |
| `069_CellOracle_GRN虚拟扰动.py` | 基于 CellOracle Oracle 对象做 GRN 虚拟敲低/敲除。 |
| `085_Squidiff_扩散模型单细胞扰动预测.py` | 调用外部 Squidiff/PerturbDiff 训练或采样脚本，记录可复现命令。 |

## 推荐输入

- Seurat RDS 或 AnnData h5ad。
- perturbation/control 分组列。
- 候选基因列表。
- 已训练模型或 GRN/Oracle 对象。

## 推荐输出

- 扰动距离、E-test 结果。
- GEARS 预测表达矩阵和运行摘要。
- CellOracle 扰动后状态转移结果。
- Squidiff/PerturbDiff 预测结果和运行日志。
