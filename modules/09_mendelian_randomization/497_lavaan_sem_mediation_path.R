# =============================================================================
# 497_lavaan_sem_mediation_path.R
# 用途    : 结构方程 / 路径中介模型——给"X→中介→结局"机制叙事一张带标准化β的通路图
#           (MR 的转录组中介搭档: MR 给遗传因果, SEM 给转录中介 → 双层证据)
# 来源    : lavaan (Rosseel 2012) + semPlot::semPaths   (CRAN)
# 补库依据 : 原 21 类无 SEM；论文2 (Biomedicines 2026, 糖尿病足) 用 lavaan 拟合
#           Transport→Inflammation→Healing 路径(标准化β=-0.92, CFI/SRMR/RMSEA)。
#           对共病主线 (CHIP→inflammaging→{癌,CVD}) 是最高叙事杠杆；同时补 SEM 路径图缺口。
# 依赖    : install.packages(c("lavaan","semPlot"))      # 先确认再装
# 输入    : df —— 含各节点"观测复合评分"(如各通路 z 标准化均值) 的 data.frame
# 输出    : 拟合指标(CFI/TLI/RMSEA/SRMR) + 标准化路径系数 + SEM 路径图(PDF)
# =============================================================================
library(lavaan); library(semPlot)

model <- '
  inflammation ~ a*exposure              # 暴露 -> 中介
  outcome      ~ b*inflammation + c*exposure   # 中介/暴露 -> 结局
  indirect := a*b                        # 中介(间接)效应
  total    := c + a*b                    # 总效应
'
# fit <- sem(model, data = df, estimator = "ML")
# summary(fit, standardized = TRUE, fit.measures = TRUE)
# parameterEstimates(fit, standardized = TRUE)
#
# semPaths(fit, what = "std", layout = "tree2", edge.label.cex = 1.1,
#          residuals = FALSE, nCharNodes = 0, sizeMan = 9)   # = 论文 Fig3B 路径图
#
# 备注: 观测复合评分当节点时为"路径分析"(SEM 的特例)；小样本 RMSEA 易偏高(论文亦如此),
#       优先看 CFI(>0.95) 与 SRMR(<0.08)。bootstrap 间接效应: sem(..., se="bootstrap").
# 与 cat09 关系: 这是"转录组中介"，与本目录 075/079 的"遗传(MR)中介"互补，共同支撑因果叙事。
