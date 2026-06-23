# 078 · FAERS 药物警戒信号挖掘

> 药物-不良事件报告 → 一条命令 → ROR/PRR/BCPNN/EBGM 四算法信号 + 森林图 + 信号热图。

| | |
|---|---|
| **语言 / 主依赖** | R · `ggplot2` |
| **一句话用途** | 不良事件不相称性分析,识别药物安全信号 |
| **输入** | `example_data/drug_event_cases.csv` |
| **输出** | `results/signals.csv` + `assets/` |

---

## ① 输入数据

二选一:① 原始报告行(列 `case_id, drug, event`);② 预计算四格计数(列 `drug, event, n11, n10, n01, n00`)。

## ② 方法 / 原理

构建药物×事件四格表 → 计算 **ROR**(报告比值比)、**PRR**(比例报告比)、**BCPNN-IC**(信息成分)、**EBGM**;以"ROR025>1 且 PRR≥2 且 IC025>0 且 n11≥3"为共识信号判定。

> 方法引用:Evans 2001(PRR);Bate 1998(BCPNN);药物警戒不相称性分析标准。

## ③ 用途

从 FAERS/JADER 等自发呈报数据库挖掘药物-不良反应信号,支持药物安全性研究与重定位风险评估。

## ④ 特点 / 亮点

- **Turnkey**:原始报告或计数均可;四算法一次算全。
- **顶刊图**:ROR 森林图(95%CI,信号着色)+ 药×事件信号热图(★标信号)。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/ROR_forest.png` | 森林图 | top 药-事件 ROR + CI |
| `assets/Signal_heatmap.png` | 热图 | log2(ROR),★=共识信号 |

![forest](assets/ROR_forest.png)

---

## 运行

```bash
Rscript 078_FAERS_pharmacovigilance.R                              # 示例
Rscript 078_FAERS_pharmacovigilance.R --input data/cases.csv
```

## 依赖安装

```r
install.packages("ggplot2")
```
