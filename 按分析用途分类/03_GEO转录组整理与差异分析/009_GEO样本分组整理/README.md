# 009 · GEO 样本分组整理 + 归一化

> 基因矩阵 + 分组表 → 一条命令 → 归一化并加分组后缀的矩阵 `Sample_Type_Matrix.csv`(直接进 010)。

| | |
|---|---|
| **语言 / 主依赖** | R · `limma`(可选 `readxl`) |
| **一句话用途** | 表达矩阵归一化 + 写入样本分组标签 |
| **输入** | `example_data/geneMatrix.csv` + `sample_group.csv` |
| **输出** | `results/Sample_Type_Matrix.csv` |

---

## ① 输入数据

| 文件 | 格式 | 必需列 | 说明 |
|------|------|------|------|
| `--expr` | csv | 首列基因名 + 样本列 | 008 产出的基因级矩阵 |
| `--group` | csv/xlsx | 第1列样本名、第2列类型 | 样本名须与 expr 列名对应;类型如 `con`/`tre` |

## ② 方法 / 原理

`limma::avereps`(重复基因平均)→ 分位数判断是否需要 → `log2(x+1)` → `normalizeBetweenArrays`(数组间归一化)→ 按分组表筛选/排序样本 → 样本名追加 `_类型` 后缀。

## ③ 用途

衔接 008 与 010:把基因矩阵标准化并打上分组标签,使下游差异分析(010)能自动按后缀识别 Control/Disease。

## ④ 特点 / 亮点

- **Turnkey**:`--expr`/`--group` 即跑;分组表 csv/xlsx 通吃。
- **自动判断 log2**:依据表达值分布自动决定是否对数化,避免重复 log。
- **链式衔接**:输出文件名/格式与 010 输入无缝对接。

## ⑤ 输出结果

无图。`results/Sample_Type_Matrix.csv`(样本名带 `_类型` 后缀)+ `Sample_Summary.txt`(各组样本数)。

---

## 运行

```bash
Rscript 009_GEO_sample_grouping.R                                       # 示例
Rscript 009_GEO_sample_grouping.R --expr geneMatrix.csv --group sample_group.csv
```

## 依赖安装

```r
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install("limma")
# 若分组表为 xlsx: install.packages("readxl")
```
