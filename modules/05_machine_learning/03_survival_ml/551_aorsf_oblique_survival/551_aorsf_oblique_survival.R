# =============================================================================
# 编号   : 551
# 脚本名 : 加速斜分裂随机生存森林 (aorsf, oblique RSF) — 生存预后 + 诚实基线
# 分类   : 12_tcga_prognosis
# 用途   : 用 aorsf 训练「斜分裂(oblique)随机生存森林」做生存预后建模:
#          每个分裂节点用变量的线性组合(而非单变量阈值)切分 → 比标准 RSF
#          更准、且可用 negation 重要性解释。全流程 train/test 划分 + 多时点评估。
# ★诚实基线 : 同一 train/test 上同时跑 CoxPH(正则线性基线)与 标准 RSF(轴对齐森林),
#             三者并列报 test C-index。大队列 ML 常仅与正则 Cox 持平 → 必报基线,
#             不只报 aorsf 好看的数字。基线对照实测见 results/cindex_comparison.csv。
# 依赖   : aorsf · survival · randomForestSRC · timeROC · ggplot2 (+ 框架 theme_pub.R)
# 运行   : Rscript 551_aorsf_oblique_survival.R                      # 合成示例(自动生成)
#          Rscript 551_aorsf_oblique_survival.R --input data/你的.csv --outdir results/run1
# 输入   : 一个 csv,含 time(随访时间,>0)、status(事件 1/删失 0)+ 任意数值/因子协变量。
#          列名 time/status 可用 --time_col / --status_col 覆盖。
# 注意   : 本环境 aorsf 多线程会段错误 → 全程 n_thread=1(稳定优先;真实大数据可调高自测)。
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(aorsf); library(survival); library(ggplot2)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input      = file.path(DDAT, "synthetic_survival.csv"),
  outdir     = file.path(SCRIPT_DIR, "results"),
  time_col   = "time",
  status_col = "status",
  n_tree     = 500,
  test_frac  = 0.30,
  horizons   = "365,730,1095"))   # time-dependent AUC 评估时点(天)
args$n_tree    <- as.integer(args$n_tree)
args$test_frac <- as.numeric(args$test_frac)
HORIZONS <- as.numeric(strsplit(as.character(args$horizons), ",")[[1]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

N_THREAD <- 1L   # 见头部注意:本环境多线程不稳

## ---- 0. 合成生存数据(若无输入) ------------------------------------------
# synthetic, for demo only:
#   12 个协变量,其中 5 个(x1..x5)真实关联生存,7 个(x6..x12)为噪声。
#   风险 = 线性主效应(x1..x5)+ 一个真实交互项(x1×x2)+ 一个非线性项(x3^2)。
#   交互/非线性结构正是斜分裂(oblique)随机森林相对 单变量 Cox / 轴对齐 RSF 的用武之地;
#   但效应温和 → 保持诚实(Cox 仍具竞争力),不人为夸大 aorsf 优势。删失 = 行政 + 随机。
if (!file.exists(args$input)) {
  set.seed(42); n <- 600; p_signal <- 5; p_noise <- 7
  X <- matrix(rnorm(n * (p_signal + p_noise)), n)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  beta <- c(0.8, -0.6, 0.5, -0.45, 0.4, rep(0, p_noise))   # 仅前 5 个真线性效应
  eta  <- as.vector(X %*% beta) +
          0.7 * X[, "x1"] * X[, "x2"] +                    # 真实交互(线性 Cox 默认捕捉不到)
          0.35 * (X[, "x3"]^2 - 1)                          # 真实非线性(中心化)
  lambda <- 0.00025 * exp(eta)                              # 个体风险率
  t_event <- rexp(n, rate = lambda)                         # 真事件时间
  t_cens  <- pmin(rexp(n, rate = 1 / 1500), 1825)           # 随机 + 5 年行政删失
  time   <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)
  df <- data.frame(time = round(time, 1), status = status, X)
  write.csv(df, args$input, row.names = FALSE)
  cat(sprintf("[gen] 合成生存数据: n=%d, 事件率=%.1f%%, 信号变量 x1-x5, 噪声 x6-x12 (synthetic demo only)\n",
              n, 100 * mean(status)))
}

## ---- 1. 读数据 + train/test 划分 -----------------------------------------
cat("Step 1: 读数据 + train/test 划分...\n")
dat <- read_table_smart(args$input)
stopifnot(args$time_col %in% names(dat), args$status_col %in% names(dat))
# 标准化为 time/status 命名,确保 status 为 0/1 整数
names(dat)[names(dat) == args$time_col]   <- "time"
names(dat)[names(dat) == args$status_col] <- "status"
dat <- dat[is.finite(dat$time) & dat$time > 0 & dat$status %in% c(0, 1), ]
dat$status <- as.integer(dat$status)
feat <- setdiff(names(dat), c("time", "status"))
# 字符列转因子(aorsf/Cox 都接受因子)
for (f in feat) if (is.character(dat[[f]])) dat[[f]] <- factor(dat[[f]])

set.seed(42)
idx_te <- sample(seq_len(nrow(dat)), floor(args$test_frac * nrow(dat)))
dtr <- dat[-idx_te, ]; dte <- dat[idx_te, ]
cat(sprintf("  n=%d (train=%d / test=%d), 特征=%d, 事件率 train=%.1f%% / test=%.1f%%\n",
            nrow(dat), nrow(dtr), nrow(dte), length(feat),
            100 * mean(dtr$status), 100 * mean(dte$status)))

# 评估时点裁剪到 test 随访范围内,避免外推
HORIZONS <- HORIZONS[HORIZONS < max(dte$time)]
if (!length(HORIZONS)) HORIZONS <- quantile(dte$time[dte$status == 1], c(.25, .5, .75))
cat("  评估时点 (天):", paste(round(HORIZONS), collapse = ", "), "\n")

# 反向 C-index:风险评分(越高越坏)与生存时间应负相关 → reverse=TRUE
cidx_risk <- function(time, status, risk)
  survival::concordance(survival::Surv(time, status) ~ risk, reverse = TRUE)$concordance

## ---- 2. aorsf 斜分裂随机生存森林 -----------------------------------------
cat("Step 2: 训练 aorsf 斜分裂随机生存森林 (oblique RSF)...\n")
fml <- stats::as.formula("Surv(time, status) ~ .")
fit_aorsf <- orsf(dtr, formula = fml, n_tree = args$n_tree, n_thread = N_THREAD,
                  importance = "negate")
risk_aorsf <- as.vector(predict(fit_aorsf, new_data = dte,
                                pred_type = "risk", pred_horizon = max(HORIZONS)))
c_aorsf <- cidx_risk(dte$time, dte$status, risk_aorsf)
cat(sprintf("  aorsf test C-index = %.3f (OOB train C = %.3f)\n",
            c_aorsf, tail(as.vector(fit_aorsf$eval_oobag$stat_values), 1)))

## ---- 3. ★诚实基线 A: CoxPH(正则线性基线) --------------------------------
cat("Step 3: 诚实基线 A — CoxPH...\n")
cox <- survival::coxph(fml, data = dtr, x = TRUE)
risk_cox <- as.vector(predict(cox, newdata = dte, type = "lp"))   # 线性预测子=风险评分
c_cox <- cidx_risk(dte$time, dte$status, risk_cox)
cat(sprintf("  CoxPH test C-index = %.3f\n", c_cox))

## ---- 4. ★诚实基线 B: 标准 RSF(轴对齐随机生存森林) -----------------------
cat("Step 4: 诚实基线 B — 标准 RSF (randomForestSRC)...\n")
has_rsf <- requireNamespace("randomForestSRC", quietly = TRUE)
if (has_rsf) {
  rsf <- randomForestSRC::rfsrc(fml, data = dtr, ntree = args$n_tree)
  risk_rsf <- predict(rsf, newdata = dte)$predicted   # 累积死亡率,越高越坏
  c_rsf <- cidx_risk(dte$time, dte$status, risk_rsf)
  cat(sprintf("  标准 RSF test C-index = %.3f\n", c_rsf))
} else {
  c_rsf <- NA_real_; risk_rsf <- NULL
  cat("  ⚠ 未装 randomForestSRC → 跳过轴对齐 RSF 基线(真实使用请安装以完整对照)。\n")
}

## ---- 诚实基线汇总表 -------------------------------------------------------
cmp <- data.frame(
  model    = c("aorsf (oblique RSF)", "CoxPH (linear baseline)", "RSF (axis-aligned)"),
  c_index  = c(c_aorsf, c_cox, c_rsf),
  stringsAsFactors = FALSE)
cmp <- cmp[!is.na(cmp$c_index), ]
write.csv(cmp, file.path(args$outdir, "cindex_comparison.csv"), row.names = FALSE)
cat("  —— 诚实基线 test C-index ——\n")
for (i in seq_len(nrow(cmp)))
  cat(sprintf("    %-26s %.3f\n", cmp$model[i], cmp$c_index[i]))
gain <- c_aorsf - c_cox
cat(sprintf("  aorsf − Cox = %+.3f  (%s)\n", gain,
            if (gain > 0.02) "斜森林明显占优" else if (gain > 0) "略优,与正则 Cox 基本持平" else "未超线性基线 → 该数据线性已够"))

## ---- 5. time-dependent AUC(各模型 × 各时点) ----------------------------
cat("Step 5: time-dependent AUC (timeROC)...\n")
suppressWarnings(suppressMessages(library(timeROC)))
mk_roc <- function(risk) timeROC::timeROC(T = dte$time, delta = dte$status,
                                          marker = risk, cause = 1, times = HORIZONS, iid = FALSE)
roc_list <- list("aorsf" = mk_roc(risk_aorsf), "CoxPH" = mk_roc(risk_cox))
if (has_rsf) roc_list[["RSF"]] <- mk_roc(risk_rsf)
auc_df <- do.call(rbind, lapply(names(roc_list), function(m)
  data.frame(model = m, time = HORIZONS, AUC = as.vector(roc_list[[m]]$AUC))))
write.csv(auc_df, file.path(args$outdir, "time_dependent_auc.csv"), row.names = FALSE)

## ---- 6. 图 1: time-dependent AUC 曲线(折线+点,非条形) -------------------
cat("Step 6: 图1 time-dependent AUC...\n")
auc_df$model <- factor(auc_df$model, levels = intersect(c("aorsf","CoxPH","RSF"), unique(auc_df$model)))
p_auc <- ggplot(auc_df, aes(time, AUC, colour = model, group = model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.6) +
  scale_color_pub("nejm") +
  scale_x_continuous(breaks = round(HORIZONS)) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(title = "Time-dependent AUC", subtitle = "aorsf vs honest baselines (test set)",
       x = "Follow-up time (days)", y = "Time-dependent AUC", colour = "Model") +
  theme_pub(base_size = 12)
save_fig(p_auc, file.path(ASSETS, "fig1_time_dependent_auc"), width = 6.4, height = 4.6)

## ---- 7. 图 2: negation 变量重要性 lollipop(非条形图) --------------------
cat("Step 7: 图2 negation 变量重要性 lollipop...\n")
vi <- orsf_vi_negate(fit_aorsf)
vi_df <- data.frame(variable = names(vi), importance = as.vector(vi))
vi_df <- vi_df[order(vi_df$importance), ]
vi_df$variable <- factor(vi_df$variable, levels = vi_df$variable)
vi_df$signal <- ifelse(vi_df$variable %in% paste0("x", 1:5), "True signal", "Noise / other")
p_vi <- ggplot(vi_df, aes(importance, variable)) +
  geom_segment(aes(x = 0, xend = importance, y = variable, yend = variable,
                   colour = signal), linewidth = 0.7) +
  geom_point(aes(colour = signal), size = 3.4) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c("True signal" = "#BC3C29", "Noise / other" = "#6F99AD")) +
  labs(title = "Negation variable importance (aorsf)",
       subtitle = "Drop in C-index when a predictor is negated; higher = more important",
       x = "Negation importance", y = NULL, colour = NULL) +
  theme_pub(base_size = 12)
save_fig(p_vi, file.path(ASSETS, "fig2_negation_importance_lollipop"), width = 6.6, height = 5)

## ---- 8. 图 3: 风险分层 KM 曲线(aorsf 风险三分位) ------------------------
cat("Step 8: 图3 风险分层 KM 曲线...\n")
grp <- cut(risk_aorsf, breaks = quantile(risk_aorsf, c(0, 1/3, 2/3, 1)),
           include.lowest = TRUE, labels = c("Low", "Intermediate", "High"))
km_dat <- data.frame(time = dte$time, status = dte$status, group = grp)
sf <- survival::survfit(survival::Surv(time, status) ~ group, data = km_dat)
sd <- survival::survdiff(survival::Surv(time, status) ~ group, data = km_dat)
logrank_p <- 1 - stats::pchisq(sd$chisq, df = length(sd$n) - 1)
# 手绘 KM(精修,非 base plot.survfit)
km_steps <- do.call(rbind, lapply(seq_along(sf$strata), function(i) {
  rng <- if (i == 1) 1:sf$strata[1] else (cumsum(sf$strata)[i-1] + 1):cumsum(sf$strata)[i]
  g   <- sub("group=", "", names(sf$strata)[i])
  data.frame(time = c(0, sf$time[rng]), surv = c(1, sf$surv[rng]), group = g)
}))
km_steps$group <- factor(km_steps$group, levels = c("Low", "Intermediate", "High"))
p_km <- ggplot(km_steps, aes(time, surv, colour = group)) +
  geom_step(linewidth = 1) +
  scale_color_pub("lancet") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Risk-stratified survival (aorsf tertiles)",
       subtitle = sprintf("Test set · log-rank p = %.2e", logrank_p),
       x = "Follow-up time (days)", y = "Survival probability", colour = "Risk group") +
  theme_pub(base_size = 12)
save_fig(p_km, file.path(ASSETS, "fig3_risk_stratified_km"), width = 6.4, height = 4.8)

## ---- 收尾 -----------------------------------------------------------------
cat("\n完成。\n")
cat("  结果表 :", normalizePath(args$outdir), "\n")
cat("    - cindex_comparison.csv  (★诚实基线对照)\n")
cat("    - time_dependent_auc.csv\n")
cat("  展示图 :", normalizePath(ASSETS), "(fig1/fig2/fig3 · PDF+PNG)\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖版本快照
