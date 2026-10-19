# 607 · Missingness and confounding diagnostics

检索词：missing data, MICE, multiple imputation, propensity score, IPW, ATE, ATT, overlap weights, SMD, effective sample size, DAG, 混杂、缺失、多重插补、倾向评分、平衡诊断。

三种独立模式：`audit` 输出逐字段缺失；`mi_glm` 对研究者指定的变量做链式多重插补、逐插补数据拟合 GLM 后用 Rubin 规则汇总；`ps_diagnostics` 用基础 R 估计 Logistic PS 和 ATE/ATT/ATO 权重，输出重叠、平衡与有效样本量。代码为原创。**不按单因素 P 值选混杂，不从缺失率判定 MAR，不默认把倾向权重乘上调查权重。** 不下载数据、不安装依赖。

```powershell
Rscript 607_missingness_confounding.R --config config.example.json
```

CLI 需要 `jsonlite`；audit 和 PS 仅依赖基础 R；MI 另需 `mice,broom`。亦可 `source()` 后调用 `run_missingness_confounding(list(...))`。输入一人一行，`factor_variables` 显式列出分类字段。输出拒绝覆盖非空目录。真实输入和模型 RDS 应按原始数据权限保存。

两个模型模式均需要 `exposure,confounders,dag_assumptions`。`confounders` 应是暴露前混杂因素的预指定集，不能无差别纳入中介、碰撞变量或暴露后的药物适应证。DAG 注释是可追溯的设计假设，本入口不会证明其正确性。官方 `dagitty::adjustmentSets` 列为 `reference_only`，未在此自动生成调整集。

MI 模式额外需要 `outcome,imputation_variables,imputation_methods,m,maxit,seed,family`。`imputation_variables` 包含暴露、结局、混杂和合理辅助变量，排除 ID；每个含 NA 的变量必须显式指定方法。例如 `"imputation_methods":{"age":"pmm","smoking":"logreg"}` 只是格式演示，真实方法需匹配变量类型和模型结构。默认预测矩阵为所有插补变量互相预测，排除对角；可用 `predictor_matrix` CSV 显式覆盖，行列名及顺序须与插补变量一致。`family` 为 `binomial` 或 `gaussian`。保留 chainMean、loggedEvents 和完整 RDS，需检查链、共线性、兼容性、结局处理与 MNAR 敏感性；这里没有“MI 通过即可发表”的自动判定。复杂抽样 MI、纵向多层 MI、MI+PS 的联合方差方案均需另行设计。

PS 模式额外需要 `estimand=ATE/ATT/ATO` 以及 `positivity_bounds`（例如 `[0.01,0.99]` 只作格式演示）；阈值用于诊断，不截断、不修剪、不改目标人群。暴露要求 0/1，PS 模型变量不得缺失。ATE 权重为按暴露比例稳定的逆概率权重；ATT 和 ATO 使用对应目标人群公式。输出每一分类水平的标准化差异，前后使用同一个未加权合并 SD 作分母；报告两组 PS 范围、阈值外人数、最大权重和 ESS。SMD 是诊断而非因果识别保证，常数或组内方差不足时记 NA。

该入口只交付权重及诊断，**未拟合加权结局效应**，因此不会把把估计权重当已知造成的方差遗漏隐藏进效应报告。需要正式效应估计时，官方 `ngreifer/WeightIt` 的 `glm_weightit()` 与 `ngreifer/cobalt::bal.tab()` 作为 `reference_only` 路线，先核对估计量、权重估计的不确定性、聚类和 survey 目标。当前本机这两个包未安装，本次不安装。

来源、commit 与许可详见 `module.json`。人工合成接口检查只能证明代码路径，不能验证 MAR、无未测混杂、正值性或糖皮质激素对癌症的因果效应。
