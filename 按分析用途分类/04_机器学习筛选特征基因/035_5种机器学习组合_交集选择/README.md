# 035 · 多方法组合交集选择

> 多种 ML 方法的特征列表 → 一条命令 → 枚举所有 k 方法组合的交集大小,选最优组合 + UpSet。

| | |
|---|---|
| **语言 / 主依赖** | R · `theme_pub` + `UpSetR` |
| **一句话用途** | 找出交集最大的方法组合 + 其一致特征 |
| **输入** | `example_data/method_sets/`(6 份方法特征列表) |
| **输出** | `results/` 组合表/交集+图 · 展示图见 `assets/` |

---

## ① 输入数据

**输入是一个目录**(`--input`),内含多份方法特征列表(csv,列名 `variable` 或首列=基因)。文件名 `importanceGene.<方法>.csv` 或任意,自动取方法名。

## ② 方法 / 原理

枚举所有 `--pick`(默认 5)方法组合 → 计算各组合交集大小 → 排序;默认选交集最大的组合(或 `--methods` 指定),输出其一致特征 + UpSet。

## ③ 用途

当方法很多时,客观选出"既稳健又保留足够特征"的方法组合,避免主观挑选。

## ④ 特点 / 亮点

- **Turnkey**:组合枚举 + 排行,一条命令;`--pick/--methods` 可调。
- **顶刊图**:top 组合交集排行榜(viridis)+ 选定组合 UpSet。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/Combo_ranking.png` | 排行榜 | top 组合交集大小 |
| `assets/Combo_UpSet.png` | UpSet | 选定组合特征交集 |
| `results/all_combinations.csv` · `selected_combo_intersection.txt` | 表 | 全组合 / 选定交集 |

![ranking](assets/Combo_ranking.png)

---

## 运行

```bash
Rscript 035_method_combo_intersection.R                                   # 示例
Rscript 035_method_combo_intersection.R --input data/method_sets --pick 5 --methods "RF,Lasso,SVM,GBM,PLS"
```

## 依赖安装

```r
install.packages("UpSetR")
```
