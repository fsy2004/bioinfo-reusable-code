# 004 · GeneCards 疾病靶点提取

> GeneCards 导出 → 一条命令 → 按相关性评分过滤的去重疾病靶点列表。

| | |
|---|---|
| **语言 / 主依赖** | R · base |
| **输入** | `example_data/GeneCards_export.csv` |
| **输出** | `results/targets.csv` |

## ① 输入数据
GeneCards 导出 CSV;自动识别 `Gene Symbol` 列与 `Relevance score` 评分列。

## ② 方法 / 原理
提取基因列 →(可选)按 `Relevance score >= --score-min`(常用 1~10)过滤 → 去重。

## ③ 用途
获取疾病相关基因,与 OMIM 等合并(→005)。

## ④ 特点 / 亮点
Turnkey;相关性阈值可调;与 001/002 同引擎。

## ⑤ 输出结果
无图。`results/targets.csv`。

## 运行
```bash
Rscript 004_extract_targets.R --score-min 5
```
