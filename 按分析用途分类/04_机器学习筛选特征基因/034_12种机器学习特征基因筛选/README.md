# 034 · 多种机器学习方法比较 + 特征筛选

> 小候选基因集 → 一条命令 → 训练多种 ML 分类器,比较 ROC/AUC,并对各法 top 特征取交集(UpSet)。

| | |
|---|---|
| **语言 / 主依赖** | R · `caret` `pROC` `UpSetR` + 各算法包 |
| **一句话用途** | 多模型横向比较 + 一致性特征筛选 |
| **输入** | `example_data/Sample_Type_Matrix.csv`(小候选集) |
| **输出** | `results/` AUC表/交集 · 展示图见 `assets/` |

---

## ① 输入数据

表达矩阵 CSV(首列基因,样本名后缀 `*_con`/`*_tre`)。**建议输入为上游交集得到的小候选集**(数十个基因),避免高维下 LDA/PLS 奇异。

## ② 方法 / 原理

`caret` 统一接口训练 10 种方法(Lasso/ElasticNet/RF/SVM/LDA/GBM/NeuralNet/PLS/kNN/LogitBoost),`repeatedcv` 调参 → 测试集 `pROC` 算 AUC;`caret::varImp` 提取各法重要特征 → 取交集。**缺失算法包的方法自动跳过**,不中断。

> 注:原脚本用 DALEX 计算重要性,本重构改用 `caret::varImp`(免重依赖,语义等价)。

## ③ 用途

回答"哪种模型对该数据最优 + 哪些特征被多数模型一致选中"。一致性特征 = 最稳健的标志物候选。

## ④ 特点 / 亮点

- **Turnkey + 稳健**:依赖缺失/单个方法失败均自动跳过并提示。
- **顶刊图**:多模型 ROC 叠加(AUC 图例)+ viridis AUC 排行榜 + 特征交集 UpSet。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/ROC_overlay.png` | ROC 叠加 | 各模型测试集 ROC + AUC |
| `assets/AUC_leaderboard.png` | 排行榜 | 模型 AUC 排序条形图 |
| `assets/Feature_UpSet.png` | UpSet | 各法 top 特征交集 |
| `results/model_AUC.csv` · `intersect_genes.txt` | 表 | AUC / 一致特征 |

![ROC](assets/ROC_overlay.png)
![AUC](assets/AUC_leaderboard.png)

---

## 运行

```bash
Rscript 034_multiML_feature_selection.R                                  # 示例
Rscript 034_multiML_feature_selection.R --input data/signature.csv --topgene 10 --train 0.7
```

## 依赖安装

```r
install.packages(c("caret","pROC","UpSetR","glmnet","randomForest","kernlab",
                   "gbm","nnet","pls","kknn","caTools","MASS"))
```
