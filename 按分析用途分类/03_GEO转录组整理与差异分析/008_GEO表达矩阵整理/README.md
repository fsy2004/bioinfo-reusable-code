# 008 · GEO 表达矩阵整理(探针 → 基因)

> GSE 序列矩阵 + GPL 平台注释 → 一条命令 → 基因级表达矩阵 `geneMatrix.csv`(供 009/010 下游)。

| | |
|---|---|
| **语言 / 主依赖** | R · base(无第三方包) |
| **一句话用途** | 把探针级 GEO 数据映射并折叠为基因级矩阵 |
| **输入** | `example_data/`(GSE 序列矩阵 + GPL 平台 txt) |
| **输出** | `results/geneMatrix.csv` |

---

## ① 输入数据

| 文件 | 格式 | 说明 |
|------|------|------|
| `GSE*_series_matrix.txt` | tab txt | 含 `ID_REF` 表头行的探针 × 样本表达表(GEO 标准下载格式) |
| `GPL*.txt` | tab txt | 平台注释,某列为基因 Symbol;列号由 `--symcol` 指定(1 起始,示例=2) |

**约定**:目录内自动识别 `GSE*`/`GPL*`;不同平台 Symbol 所在列不同,务必核对 `--symcol`。

## ② 方法 / 原理

定位 `ID_REF` 起始行读表达表 → 解析 GPL 建 探针→Symbol 映射(取 `///` 左侧、丢弃含空格的非单词)→ `merge` 对齐 → 同一基因多探针 `aggregate` 取均值 → 基因级矩阵。

## ③ 用途

GEO 微阵列分析的第一步前处理;产出标准基因矩阵,接 009(分组/归一化)→ 010(差异分析)。

## ④ 特点 / 亮点

- **Turnkey**:目录放入 GSE/GPL → 一条命令;零第三方依赖。
- **稳健**:自动定位 ID_REF;多探针按基因均值合并。

## ⑤ 输出结果

无图。`results/geneMatrix.csv`:首列 `geneSymbol` + 各样本表达。

---

## 运行

```bash
Rscript 008_GEO_expr_matrix_tidy.R                                  # 示例
Rscript 008_GEO_expr_matrix_tidy.R --gse GSExxx_series_matrix.txt --gpl GPLxxx.txt --symcol 11
```

## 依赖安装

无需额外安装(base R)。
