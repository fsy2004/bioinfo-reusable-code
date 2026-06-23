# 052 · SHAP 机器学习解释分析

> 表达矩阵 → 一条命令 → 训练多模型选最优 + kernelshap 解释 + 顶刊级 SHAP 图(重要性/蜂群/依赖/瀑布/力图)。

| | |
|---|---|
| **语言 / 主依赖** | R · `caret` `kernelshap` `shapviz` `pROC` |
| **一句话用途** | 用 SHAP 解释模型,定位并解读关键特征基因 |
| **输入** | `example_data/geneexp.csv` |
| **输出** | `results/` SHAP 表+图 · 展示图见 `assets/` |

---

## ① 输入数据

表达矩阵 CSV(首列基因,样本列名后缀分组,默认对照 `*_con`、实验 `*_tra`,可用 `--ctrl/--case` 改)。

## ② 方法 / 原理

`caret` 训练 RF/SVM/XGB → 按测试 AUC 选最优 → `kernelshap`(模型无关 Kernel SHAP)计算特征贡献 → `shapviz` 出图。

> 方法引用:Lundberg & Lee, *NeurIPS* 2017(SHAP)。

## ③ 用途

不仅筛选特征,还**解释**每个基因如何、向哪个方向影响预测(全局重要性 + 单样本归因),增强模型可信度。

## ④ 特点 / 亮点

- **Turnkey**:零改动跑示例;自动选最优模型再解释。
- **顶刊级 SHAP 全家桶**:重要性条形 · 蜂群 · 依赖 · 瀑布 · 力图,均 theme_pub 化(原"可爱配色"已升级期刊风)。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/SHAP_beeswarm.png` | 蜂群图 | 全局特征贡献分布(色=特征值) |
| `assets/SHAP_importance_bar.png` | 条形 | 平均 |SHAP| 重要性 |
| `assets/SHAP_dependence.png` | 依赖图 | top 基因 SHAP-表达关系 |
| `assets/SHAP_waterfall.png` · `SHAP_force.png` | 瀑布/力图 | 单样本归因 |
| `assets/Model_ROC.png` | ROC | 多模型比较 |

![beeswarm](assets/SHAP_beeswarm.png)

---

## 运行

```bash
Rscript 052_SHAP_interpretation.R                              # 示例
Rscript 052_SHAP_interpretation.R --input data/geneexp.csv --case _tre
```

## 依赖安装

```r
install.packages(c("caret","kernelshap","shapviz","pROC","randomForest","kernlab","xgboost"))
```
