# =============================================================================
# 编号   : 553
# 脚本名 : 生存/风险模型的诚实评估 — 校准曲线 + 决策曲线(DCA)+ 时变 AUC/Brier
# 分类   : 12_tcga_prognosis
# 用途   : 给生存预后模型(Cox/风险评分)做"超越 C-index"的诚实评估:
#          ① 区分度  time-dependent AUC(riskRegression::Score)
#          ② 校准    calibration 曲线(预测风险 vs 观测发生率)
#          ③ 临床有用性  决策曲线 DCA 的净获益 net benefit(dcurves::dca)
#          ④ 整体预测  Brier / 积分 Brier 评分 IBS(越低越好)
# ★诚实基线(2026 范式):只报 C-index/AUC 的区分度是不够的——一个高 AUC 模型
#   可能严重校准失衡、或在任何阈值下都不如"全治/不治"(无临床净获益)。本模块
#   强制把【校准 + DCA 净获益 + IBS】与区分度并列,并内置 3 个竞争模型对照
#   (null / clinical / full),让"看着好看的指标"无处遁形。
# 依赖   : riskRegression · dcurves · survival · prodlim · ggplot2 · data.table
# 运行   : Rscript 553_riskregression_dca_calibration.R            # 合成示例
#          Rscript 553_riskregression_dca_calibration.R --input data/你的.csv --time_eval 5 --outdir results/run1
# 输入   : 一个 csv,每行一个样本,至少含:
#          time(随访时间,数值) · status(0=删失,1=事件) · 若干协变量列(数值)
#          默认用合成示例(synthetic, demo only),无需任何外部输入。
# =============================================================================

## ---- 0. 定位框架 + 载顶刊主题 --------------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(survival); library(prodlim); library(riskRegression)
  library(dcurves);  library(ggplot2); library(data.table)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input     = file.path(DDAT, "survival_demo.csv"),
  outdir    = file.path(SCRIPT_DIR, "results"),
  time_col  = "time", status_col = "status",
  time_eval = "5",                 # DCA / 校准评估时点(与随访单位一致)
  eval_grid = "1,2,3,4,5,6,8"))    # 时变 AUC/Brier 的评估时点序列
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
t_eval    <- as.numeric(args$time_eval)
eval_grid <- as.numeric(strsplit(args$eval_grid, ",")[[1]])

## ---- 1. 合成示例数据(synthetic, demo only) ----------------------------------
# 设计:3 类协变量,信号强度递增,使 null<clinical<full 的真区分度有梯度,
#       且让 "full" 含一个强但真实的预后因子 → 体现"加变量确实改善净获益"。
if (!file.exists(args$input)) {
  n <- 600
  age   <- round(rnorm(n, 60, 10))                       # 临床:弱信号
  stage <- sample(1:4, n, TRUE, prob = c(.3,.3,.25,.15)) # 临床:中信号
  gene1 <- rnorm(n); gene2 <- rnorm(n); gene3 <- rnorm(n)# 分子:gene1 强,gene2 中,gene3 噪声
  # 真线性预测子(决定真风险):分子主导,临床次之,gene3 无关
  lp <- 0.020*(age-60) + 0.45*(stage-2) + 0.85*gene1 + 0.40*gene2 + 0*gene3
  base_rate <- 0.06
  u <- runif(n); t_event <- -log(u) / (base_rate * exp(lp))   # 指数比例风险
  t_cens  <- rexp(n, rate = 1/12)                             # 独立删失
  time    <- pmin(t_event, t_cens)
  status  <- as.integer(t_event <= t_cens)
  dat <- data.frame(time = round(time, 3), status = status,
                    age = age, stage = stage,
                    gene1 = round(gene1, 4), gene2 = round(gene2, 4), gene3 = round(gene3, 4))
  write.csv(dat, args$input, row.names = FALSE)
  cat(sprintf("[gen] 合成生存数据 synthetic-demo: n=%d, 事件率=%.0f%%, 中位随访=%.1f\n",
              n, 100*mean(status), median(time)))
}

## ---- 2. 读数据 + 拟合 3 个竞争 Cox 模型 ---------------------------------------
cat("Step 1: 读数据 + 拟合 null / clinical / full 三个竞争模型...\n")
D <- read_table_smart(args$input)
tc <- args$time_col; sc <- args$status_col
stopifnot(tc %in% names(D), sc %in% names(D))
D <- D[is.finite(D[[tc]]) & D[[tc]] > 0 & D[[sc]] %in% c(0,1), ]
surv_lhs <- sprintf("Surv(%s, %s)", tc, sc)

# null=仅 age;clinical=age+stage;full=age+stage+gene1+gene2(gene3 故意排除=噪声)
f_null     <- as.formula(paste(surv_lhs, "~ age"))
f_clinical <- as.formula(paste(surv_lhs, "~ age + stage"))
f_full     <- as.formula(paste(surv_lhs, "~ age + stage + gene1 + gene2"))
m_null     <- coxph(f_null,     data = D, x = TRUE)
m_clinical <- coxph(f_clinical, data = D, x = TRUE)
m_full     <- coxph(f_full,     data = D, x = TRUE)
models <- list(`Age only` = m_null, Clinical = m_clinical, Full = m_full)
cat(sprintf("  样本=%d 事件=%d;模型: Age only / Clinical(+stage) / Full(+gene1+gene2)\n",
            nrow(D), sum(D[[sc]])))

## ---- 3. riskRegression::Score —— 时变 AUC + Brier + IBS(诚实区分度+整体预测)--
cat("Step 2: Score() 计算 time-dependent AUC / Brier / IBS...\n")
sc_obj <- Score(models,
                formula  = as.formula(paste("Hist(", tc, ",", sc, ") ~ 1")),
                data     = D, times = eval_grid, conf.int = TRUE,
                metrics  = c("auc","brier"), summary = c("ibs","risks"),
                plots    = c("calibration","roc"), null.model = TRUE)
auc_dt   <- as.data.table(sc_obj$AUC$score)
brier_dt <- as.data.table(sc_obj$Brier$score)
fwrite(auc_dt,   file.path(args$outdir, "timedep_AUC.csv"))
fwrite(brier_dt, file.path(args$outdir, "timedep_Brier_IBS.csv"))
# 取评估时点的汇总(IBS 取评估网格内累计到 t_eval 最近点)
ibs_at <- brier_dt[abs(times - t_eval) == min(abs(times - t_eval)),
                   .(model, IBS)][!duplicated(model)]
cat("  [诚实指标] 评估时点", t_eval, "附近 IBS(越低越好):\n")
print(ibs_at)

## ---- 4. 校准曲线数据(预测风险 vs 观测发生率) --------------------------------
cat("Step 3: 校准曲线(t =", t_eval, ")...\n")
pcal <- plotCalibration(sc_obj, times = t_eval, cens.method = "local",
                        method = "quantile", q = 10, plot = FALSE)
cal_df <- rbindlist(lapply(names(pcal$plotFrames), function(nm) {
  fr <- as.data.frame(pcal$plotFrames[[nm]]); data.table(model = nm, Pred = fr$Pred, Obs = fr$Obs)
}))
fwrite(cal_df, file.path(args$outdir, "calibration_points.csv"))

## ---- 5. 决策曲线 DCA —— 临床净获益(dcurves) ---------------------------------
cat("Step 4: 决策曲线 DCA(net benefit, t =", t_eval, ")...\n")
risk_df <- data.frame(time = D[[tc]], status = D[[sc]])
for (nm in names(models)) {
  risk_df[[nm]] <- as.numeric(predictRisk(models[[nm]], newdata = D, times = t_eval))
}
mod_terms <- paste(sprintf("`%s`", names(models)), collapse = " + ")
dca_obj <- dca(as.formula(paste0("Surv(time, status) ~ ", mod_terms)),
               data = risk_df, time = t_eval,
               thresholds = seq(0, 0.5, by = 0.01))
dca_dt <- as.data.table(dca_obj$dca)[, .(label, threshold, net_benefit)]
fwrite(dca_dt, file.path(args$outdir, "dca_net_benefit.csv"))

## =============================================================================
## 出图(全部顶刊风格;禁用平凡条形图 → 折线/带状/dumbbell/lollipop)
## =============================================================================
cat("Step 5: 出图(校准 / DCA / 时变 AUC / Brier-IBS)...\n")
PAL <- pal_pub(name = "npg")
cols <- setNames(PAL[seq_along(models)], names(models))

## (A) 校准曲线:对角线=完美校准;每模型一条带点折线 ---------------------------
p_cal <- ggplot(cal_df, aes(Pred, Obs, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = cols, name = "Model") +
  coord_equal(xlim = c(0, max(cal_df$Pred, cal_df$Obs)), ylim = c(0, max(cal_df$Pred, cal_df$Obs))) +
  labs(title = "Calibration curves",
       subtitle = sprintf("Predicted vs observed event risk at t = %g (quantile bins)", t_eval),
       x = "Predicted risk", y = "Observed event rate") +
  theme_pub()
save_fig(p_cal, file.path(ASSETS, "calibration_curve"), width = 6.2, height = 6.2)

## (B) 决策曲线 DCA:净获益 vs 阈概率;含 Treat All / Treat None 参照 -----------
dca_dt[, label := factor(label, levels = c("Treat All","Treat None", names(models)))]
ref_cols <- c(`Treat All` = "grey40", `Treat None` = "grey70")
dca_cols <- c(ref_cols, cols)
p_dca <- ggplot(dca_dt, aes(threshold, net_benefit, colour = label, linetype = label)) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = dca_cols, name = NULL) +
  scale_linetype_manual(values = c(`Treat All` = "dashed", `Treat None` = "dotted",
                                   setNames(rep("solid", length(models)), names(models))), name = NULL) +
  coord_cartesian(ylim = c(min(-0.02, min(dca_dt$net_benefit)), max(dca_dt$net_benefit) * 1.05)) +
  labs(title = "Decision curve analysis",
       subtitle = sprintf("Clinical net benefit across threshold probabilities at t = %g", t_eval),
       x = "Threshold probability", y = "Net benefit") +
  theme_pub()
save_fig(p_dca, file.path(ASSETS, "decision_curve"), width = 7, height = 5.2)

## (C) 时变 AUC:折线 + 95%CI 带 ----------------------------------------------
auc_plot <- auc_dt[model != "Null model"]
auc_plot[, model := factor(model, levels = names(models))]
p_auc <- ggplot(auc_plot, aes(times, AUC, colour = model, fill = model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = cols, name = "Model") +
  scale_fill_manual(values = cols, guide = "none") +
  labs(title = "Time-dependent discrimination (AUC)",
       subtitle = "Cumulative/dynamic AUC with 95% CI; 0.5 = no discrimination",
       x = "Time", y = "AUC(t)") +
  theme_pub()
save_fig(p_auc, file.path(ASSETS, "timedep_auc"), width = 6.8, height = 5)

## (D) 时变 Brier:折线(越低越好)+ IBS lollipop 内插对照 ----------------------
brier_plot <- brier_dt[, model := factor(model, levels = c("Null model", names(models)))]
p_brier <- ggplot(brier_plot, aes(times, Brier, colour = model)) +
  geom_line(linewidth = 0.85) + geom_point(size = 1.9) +
  scale_colour_manual(values = c(`Null model` = "grey55", cols), name = "Model") +
  labs(title = "Time-dependent prediction error (Brier)",
       subtitle = "Lower is better; Null model = Kaplan-Meier reference",
       x = "Time", y = "Brier score") +
  theme_pub()
save_fig(p_brier, file.path(ASSETS, "timedep_brier"), width = 6.8, height = 5)

## (E) IBS lollipop:跨模型整体预测误差对比(代替条形图)-----------------------
ibs_dt <- brier_dt[abs(times - t_eval) == min(abs(times - t_eval)),
                   .(model, IBS)][!duplicated(model)]
ibs_dt[, model := factor(model, levels = rev(c("Null model", names(models))))]
p_ibs <- ggplot(ibs_dt, aes(IBS, model, colour = model)) +
  geom_segment(aes(x = 0, xend = IBS, yend = model), linewidth = 0.9) +
  geom_point(size = 4.2) +
  geom_text(aes(label = sprintf("%.4f", IBS)), hjust = -0.3, size = 3.3, colour = "black") +
  scale_colour_manual(values = c(`Null model` = "grey55", cols), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Integrated Brier Score (overall prediction error)",
       subtitle = sprintf("IBS up to t = %g; lower = better calibrated + discriminating", t_eval),
       x = "IBS", y = NULL) +
  theme_pub()
save_fig(p_ibs, file.path(ASSETS, "ibs_lollipop"), width = 6.6, height = 3.6)

## ---- 6. 诚实评估汇总表(区分度 + 校准 + 净获益 + 整体)-------------------------
auc_te <- auc_dt[abs(times - t_eval) == min(abs(times - t_eval)) & model != "Null model",
                 .(model, AUC = round(AUC, 3))][!duplicated(model)]
# 模型在 0.2 阈值处的净获益(代表性阈值,临床有用性)
nb_te <- dca_dt[abs(threshold - 0.2) < 1e-9 & label %in% names(models),
                .(model = as.character(label), NetBenefit_at_0.2 = round(net_benefit, 4))]
summ <- Reduce(function(a,b) merge(a,b,by="model",all=TRUE),
               list(auc_te, ibs_dt[model %in% names(models), .(model=as.character(model), IBS=round(IBS,4))], nb_te))
fwrite(summ, file.path(args$outdir, "honest_eval_summary.csv"))
cat("\n[★诚实评估汇总] 区分度(AUC) + 整体(IBS) + 临床有用性(净获益@0.2):\n")
print(summ)
cat("\n完成。结果表见", normalizePath(args$outdir), ";展示图见 assets/\n")

## ---- 依赖版本快照(铁律 6) ---------------------------------------------------
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
