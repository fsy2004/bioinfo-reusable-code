# 15_药物扰动_药物重定位

本模块用于把单细胞状态、候选靶点和真实世界药物信号连接到药物响应预测、药物重定位和药物警戒。

## 脚本

| 脚本 | 作用 |
|---|---|
| `070_chemCPA_药物扰动表达预测.py` | 调用 chemCPA 预测单细胞药物扰动表达响应。 |
| `071_scDrug_单细胞药物响应预测.py` | 调用 scDrug 进行 cluster-level 药敏和治疗组合筛选。 |
| `078_FAERS_ROR_PRR_BCPNN_EBGM药物警戒.R` | 对 drug-event 表做 ROR、PRR、BCPNN/IC 和 EBGM proxy 信号挖掘。 |

## 推荐输入

- AnnData h5ad。
- 细胞类型或 cluster 标签。
- 候选药物、候选靶点或扰动设置。
- FAERS drug-event-case 原始表或 n11/n10/n01/n00 计数表。

## 推荐输出

- 药物扰动表达响应预测。
- cluster-level IC50/AUC、drug kill efficacy 和治疗组合报告。
- FAERS 不成比例报告信号表和阳性信号表。
