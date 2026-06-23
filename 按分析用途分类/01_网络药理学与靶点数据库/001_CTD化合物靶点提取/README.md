# 001 · CTD 化合物靶点提取

> CTD 数据库导出文件 → 一条命令 → 去重靶点基因列表(供下游 Venn/富集)。

| | |
|---|---|
| **语言 / 主依赖** | R · base |
| **输入** | `example_data/CTD_export.csv` |
| **输出** | `results/targets.csv` |

## ① 输入数据
CTD 导出 CSV;自动识别 `Gene Symbol` 列(及评分列 `Reference Count`)。

## ② 方法 / 原理
读取导出表 → 提取基因列 → 可选按评分过滤 → 去重 → 靶点列表。

## ③ 用途
网络药理学第一步:把 CTD 化合物-基因关联整理成标准靶点列表。

## ④ 特点 / 亮点
Turnkey;自动识别基因列与评分列;零第三方依赖。

## ⑤ 输出结果
无图。`results/targets.csv`(去重靶点)。

## 运行
```bash
Rscript 001_extract_targets.R                                  # 示例
Rscript 001_extract_targets.R --input data/CTD_export.csv
```
