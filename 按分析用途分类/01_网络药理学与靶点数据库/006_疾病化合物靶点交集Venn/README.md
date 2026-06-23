# 006 · 疾病 × 化合物靶点交集 Venn

> 疾病靶点 + 化合物靶点 → 一条命令 → 交集(=核心作用靶点)+ Venn。

| | |
|---|---|
| **语言 / 主依赖** | R · `theme_pub` + `UpSetR` |
| **输入** | `example_data/`(disease + compound 靶点 csv) |
| **输出** | `results/` 交集表 + `assets/` 图 |

## ① 输入数据
`--input` 目录,含疾病靶点与化合物靶点列表(csv;自动识别 `Gene`/首列)。

## ② 方法 / 原理
求两集交集 = **疾病-化合物共同靶点**(网络药理学核心作用靶点)→ `venn_pub` + 柱状图。

## ③ 用途
网络药理学关键一步:化合物靶点 ∩ 疾病靶点 = 化合物治疗该病的候选作用靶点,供下游 PPI/富集。

## ④ 特点 / 亮点
Turnkey;零依赖 Venn;交集基因直接输出供 PPI/富集(→007)。

## ⑤ 输出结果图
`assets/Target_Venn.png`(Venn)· `assets/Set_size_bar.png`

![Venn](assets/Target_Venn.png)

## 运行
```bash
Rscript 006_disease_compound_venn.R
```
依赖:`install.packages("UpSetR")`
