# 606 · Longitudinal incident cancer

检索词：CHARLS, incident cancer, interval censoring, AFT, discrete time, cause-specific hazard, competing death, 队列、新发肿瘤、区间删失、竞争风险。

用途：在研究者已经定义癌症结局和时间轴后，用三种互斥入口处理不同随访结构。每个入口核对基线癌症，排除基线患癌和基线状态未知者并记录人数。**同一份数据应按事件时间的真实分辨率选择模型，不能把访谈波次发现癌症的日期当精确发病日期。** 不下载数据、不安装依赖。

```powershell
Rscript 606_longitudinal_incident_cancer.R --config config.example.json
```

依赖 `survival`，JSON CLI 另需 `jsonlite`；`source()` 后可传 list 给 `run_incident_cancer()`。必需配置为示例所示字段，协变量由研究设计确定。输出拒绝覆盖非空目录。所有模式使用 `status: 0=该时点未观察到癌症事件/截尾,1=新发癌症,2=竞争死亡`；禁止把失访当成未发癌的永久状态。缺失模型变量会停止，需要研究者提前确定 MI 或 complete-case 纳排。

| mode | 时间输入与单位 | 估计和边界 |
|---|---|---|
| `cause_specific_cox` | 每人一行 `entry_time,exit_time,status`，自共同时间原点计时 | `Surv(entry,exit,status==1)`；支持延迟入组，竞争死亡在其死亡时刻退出风险集；输出原因别癌症 HR 与 `cox.zph`。不输出累计发病概率或 Fine–Gray HR。 |
| `interval_aft` | 每人一行 `left_time,right_time,status`。癌症发生于两次评估之间；精确事件可左右相等；截尾 `right_time=NA`、`left_time=末次无癌随访时间` | `Surv(type='interval2')` + `survreg`。另需 `aft_distribution=weibull/lognormal/loglogistic`。输出时间比，不能解释成 HR；本入口不支持左截断或竞争死亡，存在 status=2 会停止。 |
| `discrete_time` | 每人每风险区间一行 `period_start,period_end,period,status`，同一人区间连续且不重叠，事件/死亡只在最后一行 | `glm(cloglog)` + 波次因子 + 受试者聚类 sandwich SE；只接受基线固定暴露/协变量。竞争死亡的最终区间无癌症事件，随后无风险行。不假设行相互独立，不恢复精确发病日期。 |

离散入口**强制同一 `period` 在不同人之间具有完全相同的起止边界**；各 period 之间的长度可以不同，其差异由 period 因子吸收。不规则的个体访谈区间会拒绝运行，应另行设计时间模型。首次自报癌症可能有回忆误差，既往不明不能默认无癌。AFT 的分布和截尾假设要另做适合性/敏感性评估；若有竞争死亡，不能先删除死者来绕过保护。

输出 `cohort_flow.csv`、`coefficients.csv`、模型 RDS、配置、输入 MD5、`sessionInfo()` 和解释边界。Cox 和 AFT 均不自动应用 CHARLS 的横断/纵向权重或失访权重；这些权重的目标人群和方差方案需额外设计。原因别 Cox、Fine–Gray 和累计发生函数回答不同问题；官方 `cmprsk::crr/cuminc` 仅列为 `reference_only`，本代码未实现其封装。药物适应证混杂、时变用药、反向因果、潜伏期和 immortal-time bias 需在上游研究方案中处理。

代码为原创适配器，使用 `therneau/survival` 官方 API，固定 commit/许可记录在 `module.json`。人工合成接口检查只检验代码路径，不构成 CHARLS 实际研究或科学结论。
