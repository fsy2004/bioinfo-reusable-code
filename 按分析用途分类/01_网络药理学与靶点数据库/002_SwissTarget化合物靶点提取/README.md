# 002 · SwissTargetPrediction 化合物靶点提取

> SwissTargetPrediction 导出 → 一条命令 → 按概率过滤的去重靶点列表。

| | |
|---|---|
| **语言 / 主依赖** | R · base |
| **输入** | `example_data/SwissTargetPrediction_export.csv` |
| **输出** | `results/targets.csv` |

## ① 输入数据
SwissTarget 导出 CSV;自动识别 `Gene` 列与 `Probability` 评分列。

## ② 方法 / 原理
提取基因列 →(可选)按 `Probability >= --score-min`(常用 0.1)过滤 → 去重。

## ③ 用途
获取化合物的预测靶点,与 CTD 等合并(→003)。

## ④ 特点 / 亮点
Turnkey;概率阈值可调;与 001/004 同引擎。

## ⑤ 输出结果
无图。`results/targets.csv`。

## 运行
```bash
Rscript 002_extract_targets.R --score-min 0.1
```
