# 605 · Survey cancer association models

检索词：NHANES, complex survey, svyglm, svycoxph, interaction, natural cubic spline, 患癌现况、复杂抽样。

用途：设计权重描述、癌症现患二元结局的 Logistic 关联、乘法尺度交互、预指定自然三次样条，以及有外部死亡随访时的加权 Cox。**横断面癌症史不能给出肿瘤发生率；癌症死亡不是肿瘤新发。** 本模块不下载数据、不猜测周期合并权重、不选择混杂因素、不安装依赖。

```powershell
Rscript 605_survey_cancer_models.R --config config.example.json
```

运行前修改示例路径和字段。依赖 `survey,survival`；JSON 命令行入口另外依赖 `jsonlite`。亦可 `source()` 后调用 `run_survey_cancer(list(...))`。所有依赖由研究者预先准备。输入每个 ID 一行；设计变量不许缺失，权重正且有限，结局显式为 0/1/NA。先在完整设计上 `svydesign()`，再对可选的 0/1 `domain` 列限定分析域。`lonely_psu` 必须明确选择；示例 `fail` 让研究者处理设计问题。`weight_provenance` 记录周期、抽样子集、权重及官方规则，示例文字须替换。

`covariates` 为已由研究问题和 DAG 确定的调整变量，示例不是通用调整集。`modifier` 为可选交互变量，输出整个交互项的 `regTermTest()`。自然样条可加入 `"spline":{"knots":[5,10],"boundary_knots":[0,20]}`，数值仅演示 JSON 格式，真实结点必须预指定并匹配量纲；输出基函数系数，不把各基函数指数化值当暴露水平对照的 OR，不自动宣称非线性显著。

死亡模型可加入 `"mortality":{"time":"followup_months","status":"death","endpoint":"all_cause_mortality"}` 或 `cancer_mortality`。须先核对链接适格人群、行政截尾、死亡结局定义和时间单位。癌症死亡的其他死因是否按原因别风险删失必须在上游明确编码；本入口没有 Fine–Gray 模型。时间不允许负值。缺失结局或协变量由各模型显式 complete-case 处理，`used_n` 逐模型导出。

输出系数/区间、现患率分组表、模型 RDS、描述和交互检验文本、配置、输入 MD5 与 `sessionInfo()`；已有非空结果目录会拒绝覆盖。指数化 Logistic 系数是现患 OR，Cox 系数是 HR。RDS 含参与者变量，真实分析时应按原始数据权限保存。

源 API：`survey::svydesign/svymean/svyglm/regTermTest/svycoxph`、`survival::Surv`、`splines::ns`。GitHub 镜像与官方源的身份、commit 和许可在 `module.json` 中分别标明。人工合成接口检查不构成 NHANES 实际分析或因果证据。
