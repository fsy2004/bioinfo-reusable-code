# 013 · SVM-RFE 特征基因筛选

> 表达矩阵 + 候选基因 → 一条命令 → SVM 递归特征消除排序 + CV 准确率曲线选最优子集。

| | |
|---|---|
| **语言 / 主依赖** | R · `e1071` `ggplot2` |
| **一句话用途** | 用线性 SVM 权重递归排序、选最优特征数 |
| **输入** | `example_data/Sample_Type_Matrix.csv` + `candidate_genes.csv` |
| **输出** | `results/` 排名/子集+图 · 展示图见 `assets/` |

---

## ① 输入数据

同 [012](../012_LASSO特征基因筛选/):表达矩阵(样本名后缀分组)+ 可选候选基因。

## ② 方法 / 原理

线性核 SVM 拟合 → 以权重平方为特征评分,每轮剔除评分最低者(SVM-RFE 递归消除)得到完整排名 → 对 top-k 子集做 k 折交叉验证,准确率最高处为最优特征数。

> 方法引用:Guyon *et al.*, *Machine Learning* 2002(SVM-RFE)。

## ③ 用途

与 LASSO/RF 互补的特征选择;给出"最少多少基因即可达到最佳判别"的答案。

## ④ 特点 / 亮点

- **Turnkey**:零改动跑示例;`--maxk/--folds` 可调。
- **顶刊图**:CV 准确率-特征数曲线(标注最优 n)+ RFE 排名图(选中/未选中着色)。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/SVMRFE_CV_accuracy.png` | 曲线 | CV 准确率 vs 特征数,标注最优 |
| `assets/SVMRFE_top_rank.png` | 排名图 | top 特征 RFE 排名 |
| `results/SVMRFE_ranking.csv` · `SVMRFE_selected_genes.txt` | 表 | 全排名 / 最优子集 |

![CV](assets/SVMRFE_CV_accuracy.png)

---

## 运行

```bash
Rscript 013_SVM_RFE_feature_selection.R                                  # 示例
Rscript 013_SVM_RFE_feature_selection.R --input data/expr.csv --maxk 30 --folds 5
```

## 依赖安装

```r
install.packages(c("e1071","ggplot2"))
```
