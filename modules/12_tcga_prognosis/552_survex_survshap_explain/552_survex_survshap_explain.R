# =============================================================================
# 编号   : 552
# 脚本名 : 生存模型可解释 — SurvSHAP(t) 时间依赖 SHAP + survex 解释 (RSF/Cox)
# 分类   : 12_tcga_prognosis
# 用途   : 对生存预后模型(Cox 比例风险 / 随机生存森林 RSF)做时间依赖可解释:
#          ① survex::explain() 统一封装模型 → ② model_survshap() 算 SurvSHAP(t)
#          时间依赖 SHAP(每个特征对生存函数的贡献随时间 t 变化的曲线)
#          ③ predict_parts(type="survshap"/"survlime") 出单样本 BD/SurvLIME profile
#          ④ model_performance() 报时间依赖 C-index / Brier / AUC。
# ★诚实基线(可解释基线,非统计显著性): 用 model_parts() 的「时间积分全局排列重要性」
#          (每特征一个标量)作对照,展示它把随时间变化的贡献压成一个数 → 会掩盖
#          「某特征早期重要、晚期失效」之类的动态;SurvSHAP(t) 把同一信息展开成
#          时间曲线,这正是它相对单一全局重要性的增量价值。脚本实测两者并并排出图。
# 依赖   : survex(>=1.2) · survival · ranger · ggplot2 (+ 框架 theme_pub.R)
# 运行   : Rscript 552_survex_survshap_explain.R                       # 合成示例,零改动即跑
#          Rscript 552_survex_survshap_explain.R --input data/你的.csv --time_col OS.time \
#                  --event_col OS --model rsf --outdir results/run1
# 输入   : 一张生存数据表(csv/tsv):每行一个样本,含 时间列 + 事件列(0/1) + ≥1 个数值特征列。
#          列名通过 --time_col / --event_col 指定;其余数值列自动作为特征(可 --features a,b,c 限定)。
#          合成示例数据写到 example_data/synthetic_survival.csv(synthetic, for demo only)。
# =============================================================================

## ---- 定位并加载顶刊绘图框架 ----------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(survex); library(survival); library(ggplot2)
  has_ranger <- requireNamespace("ranger", quietly = TRUE)
  if (has_ranger) library(ranger)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input     = file.path(DDAT, "synthetic_survival.csv"),
  time_col  = "time",
  event_col = "status",
  features  = "",                    # 逗号分隔限定特征列;留空=自动取除时间/事件外的数值列
  model     = "cox",                 # cox | rsf
  n_shap    = "60",                  # 算 SurvSHAP(t) 的样本数(背景计算量,越大越慢)
  outdir    = file.path(SCRIPT_DIR, "results")))
N_SHAP <- as.integer(args$n_shap)
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 0. 合成生存示例(若无输入文件) -------------------------------------
# 5 个特征,刻意制造「时间依赖效应」以凸显 SurvSHAP(t) 价值:
#   age      : 全程稳定有害(线性进入风险)
#   stage    : 早期主导风险(早期事件由它驱动),晚期减弱  → 时间依赖!
#   biomarker: 保护性(高表达降风险),中后期才显现
#   noise1/2 : 纯噪声,期望贡献≈0(作阴性对照,检验解释不虚报)
if (!file.exists(args$input)) {
  n <- 300
  age       <- rnorm(n, 0, 1)
  stage     <- rnorm(n, 0, 1)
  biomarker <- rnorm(n, 0, 1)
  noise1    <- rnorm(n, 0, 1)
  noise2    <- rnorm(n, 0, 1)
  # 基础线性预测子(决定整体风险高低):age 有害、biomarker 保护
  lp_base <- 0.85 * age - 0.70 * biomarker
  # stage 的危害「早期强」: 用一个早发的指数事件时钟体现 —— stage 高者其潜在事件时间被大幅压缩
  u  <- runif(n)
  t_base  <- rexp(n, rate = exp(lp_base) * 0.5)          # age/biomarker 驱动的基线事件时间
  t_stage <- rexp(n, rate = exp(1.1 * stage) * 1.5)      # stage 驱动的、整体更早的事件时间
  time <- pmin(t_base, t_stage)                          # 取最早者 → 早期事件更多由 stage 决定
  # 行政删失(约 30%)
  cens <- quantile(time, 0.70)
  status <- as.integer(time <= cens); time <- pmin(time, cens)
  syn <- data.frame(time = round(time, 4), status = status,
                    age = round(age, 4), stage = round(stage, 4),
                    biomarker = round(biomarker, 4),
                    noise1 = round(noise1, 4), noise2 = round(noise2, 4))
  writeLines("# synthetic survival data, for demo only (552_survex_survshap_explain)",
             con = file.path(DDAT, "_README_synthetic.txt"))
  write.csv(syn, args$input, row.names = FALSE)
  cat(sprintf("[gen] 合成生存数据 %d 样本 (%d 事件) → %s\n",
              n, sum(status), basename(args$input)))
}

## ---- 1. 读数据 + 组织 X / Surv(y) ----------------------------------------
cat("Step 1: 读生存数据 → 构造特征矩阵 X 与 Surv(y)...\n")
dat <- read_table_smart(args$input)
stopifnot(args$time_col %in% names(dat), args$event_col %in% names(dat))
ytime <- as.numeric(dat[[args$time_col]]); yevent <- as.integer(dat[[args$event_col]])
if (nzchar(args$features)) {
  feats <- trimws(strsplit(args$features, ",")[[1]])
} else {
  num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  feats <- setdiff(num_cols, c(args$time_col, args$event_col))
}
X <- dat[, feats, drop = FALSE]
keep <- is.finite(ytime) & ytime > 0 & yevent %in% c(0L, 1L) & stats::complete.cases(X)
X <- X[keep, , drop = FALSE]; ytime <- ytime[keep]; yevent <- yevent[keep]
ysurv <- survival::Surv(ytime, yevent)
cat(sprintf("  样本=%d  事件=%d (%.0f%%)  特征(%d): %s\n",
            nrow(X), sum(yevent), 100*mean(yevent), length(feats), paste(feats, collapse=", ")))

## ---- 2. 拟合生存模型(Cox 或 RSF) ---------------------------------------
cat(sprintf("Step 2: 拟合 %s 生存模型...\n", toupper(args$model)))
mdf <- cbind(X, .time = ytime, .status = yevent)
fml <- stats::as.formula(paste0("Surv(.time, .status) ~ ", paste(feats, collapse = " + ")))
if (args$model == "rsf") {
  if (!has_ranger) stop("--model rsf 需要 ranger 包: install.packages('ranger')")
  fit <- ranger::ranger(fml, data = mdf, num.trees = 300,
                        respect.unordered.factors = TRUE, seed = 42)
} else {
  fit <- survival::coxph(fml, data = mdf, x = TRUE, model = TRUE)
}

## ---- 3. survex::explain() 统一封装 + 时间依赖性能 ------------------------
cat("Step 3: survex::explain() 封装 → 时间依赖性能 (C-index / Brier / AUC)...\n")
expl <- survex::explain(fit, data = X, y = ysurv,
                        label = toupper(args$model), verbose = FALSE)
perf <- survex::model_performance(expl)
cidx <- as.numeric(perf$result[["C-index"]])
ibs  <- as.numeric(perf$result[["Integrated Brier score"]])
cat(sprintf("  模型 C-index=%.3f · Integrated Brier=%.3f · 评估时间点=%d\n",
            cidx, ibs, length(expl$times)))

## ---- 4. ★诚实基线:全局排列重要性(时间积分,每特征一个标量) ----------
# model_parts = 时间依赖排列重要性;对时间维做积分即得「单一全局重要性」——
# 这是 SurvSHAP(t) 要超越的对象。下面同时保留它的时间分辨版本以便对比。
cat("Step 4: [诚实基线] 全局排列重要性 model_parts() ...\n")
mparts <- survex::model_parts(expl, N = min(80, nrow(X)))
mp_df  <- mparts$result
mp_times <- mp_df[["_times_"]]
# 时间积分(梯形法)→ 每特征一个全局重要性标量
trap_int <- function(tt, vv) sum(diff(tt) * (head(vv, -1) + tail(vv, -1)) / 2) / (max(tt) - min(tt))
global_imp <- sapply(feats, function(f) trap_int(mp_times, mp_df[[f]]))
global_imp_df <- data.frame(feature = names(global_imp),
                            importance = as.numeric(global_imp))
global_imp_df <- global_imp_df[order(global_imp_df$importance), ]
write.csv(global_imp_df, file.path(args$outdir, "baseline_global_permutation_importance.csv"),
          row.names = FALSE)
cat("  全局(时间积分)重要性排序: ",
    paste(sprintf("%s=%.4f", rev(global_imp_df$feature), rev(global_imp_df$importance)),
          collapse = " > "), "\n")

## ---- 5. SurvSHAP(t):时间依赖 SHAP(核心增量) ---------------------------
cat(sprintf("Step 5: SurvSHAP(t) model_survshap() (n=%d 样本)...\n", min(N_SHAP, nrow(X))))
idx <- sample(seq_len(nrow(X)), min(N_SHAP, nrow(X)))
ss <- survex::model_survshap(expl, new_observation = X[idx, , drop = FALSE])
shap_times <- ss$eval_times
# 跨样本平均 |SHAP|(t):每特征一条「时间依赖重要性」曲线
mean_abs_shap <- Reduce(`+`, lapply(ss$result, function(d) abs(d[feats]))) / length(ss$result)
shap_curve <- data.frame(time = rep(shap_times, length(feats)),
                         feature = rep(feats, each = length(shap_times)),
                         mean_abs_shap = unlist(mean_abs_shap[feats], use.names = FALSE))
write.csv(shap_curve, file.path(args$outdir, "survshap_t_mean_abs.csv"), row.names = FALSE)
# 时间积分 SurvSHAP 重要性(用于与基线对比 & 排序)
shap_global <- sapply(feats, function(f) {
  v <- shap_curve$mean_abs_shap[shap_curve$feature == f]; trap_int(shap_times, v) })
shap_global_df <- data.frame(feature = names(shap_global),
                             importance = as.numeric(shap_global))
shap_global_df <- shap_global_df[order(shap_global_df$importance), ]
write.csv(shap_global_df, file.path(args$outdir, "survshap_time_integrated_importance.csv"),
          row.names = FALSE)

## ---- 6. 单样本 profile:SurvSHAP(t) + SurvLIME ---------------------------
cat("Step 6: 单样本 predict_parts (survshap + survlime)...\n")
foc <- idx[1]                                  # 取一个焦点样本
pp_shap <- survex::predict_parts(expl, new_observation = X[foc, , drop = FALSE],
                                 type = "survshap")
single_shap <- data.frame(time = rep(pp_shap$eval_times, length(feats)),
                          feature = rep(feats, each = length(pp_shap$eval_times)),
                          shap = unlist(pp_shap$result[feats], use.names = FALSE))
write.csv(single_shap, file.path(args$outdir, "single_obs_survshap_t.csv"), row.names = FALSE)
pp_lime <- tryCatch(
  survex::predict_parts(expl, new_observation = X[foc, , drop = FALSE], type = "survlime"),
  error = function(e) { cat("  ⚠ SurvLIME 失败:", conditionMessage(e), "\n"); NULL })
lime_df <- NULL
if (!is.null(pp_lime)) {
  lime_coef <- as.numeric(pp_lime$result[1, feats])
  lime_df <- data.frame(feature = feats, coef = lime_coef)
  write.csv(lime_df, file.path(args$outdir, "single_obs_survlime_coef.csv"), row.names = FALSE)
}

## ===========================================================================
## 出图(全部 lollipop / 时间曲线 / dumbbell / heatmap;无平凡条形图)
## ===========================================================================
cat("Step 7: 顶刊级出图(SurvSHAP(t) 曲线 / 基线对比 / 单样本 profile / 热图)...\n")
pal <- pal_pub(length(feats), "npg"); names(pal) <- feats

## (1) SurvSHAP(t) 时间依赖解释曲线 —— 核心图,展示贡献随时间变化
p1 <- ggplot(shap_curve, aes(time, mean_abs_shap, colour = feature)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pal, name = "Feature") +
  labs(title = "SurvSHAP(t): time-dependent feature importance",
       subtitle = sprintf("Mean |SHAP| over %d samples · %s model (C-index=%.2f)",
                          length(idx), toupper(args$model), cidx),
       x = "Survival time t", y = "Mean |SHAP(t)|") +
  theme_pub(base_size = 12)
save_fig(p1, file.path(ASSETS, "01_survshap_t_curves"), width = 7.2, height = 5)

## (2) ★诚实基线对比:全局(时间积分)重要性 vs SurvSHAP —— dumbbell
##     两法各自 min-max 归一到 [0,1] 再并排,凸显「单一标量 vs 时间展开」之别
norm01 <- function(v) { r <- range(v); if (diff(r) == 0) rep(0.5, length(v)) else (v - r[1]) / diff(r) }
cmp <- merge(
  data.frame(feature = global_imp_df$feature, perm = norm01(global_imp_df$importance)),
  data.frame(feature = shap_global_df$feature, shap = norm01(shap_global_df$importance)),
  by = "feature")
cmp <- cmp[order(cmp$shap), ]; cmp$feature <- factor(cmp$feature, levels = cmp$feature)
cmp_long <- rbind(
  data.frame(feature = cmp$feature, method = "Permutation (global, time-integrated)", value = cmp$perm),
  data.frame(feature = cmp$feature, method = "SurvSHAP (time-integrated)", value = cmp$shap))
p2 <- ggplot() +
  geom_segment(data = cmp, aes(x = perm, xend = shap, y = feature, yend = feature),
               colour = "grey70", linewidth = 1.1) +
  geom_point(data = cmp_long, aes(value, feature, colour = method), size = 4) +
  scale_colour_manual(values = c("Permutation (global, time-integrated)" = "#3C5488",
                                 "SurvSHAP (time-integrated)" = "#E64B35"), name = NULL) +
  labs(title = "Honest baseline: single global score vs time-resolved SHAP",
       subtitle = "Permutation collapses time into one number; SurvSHAP(t) keeps the time axis (see fig.1)",
       x = "Normalized importance [0,1]", y = NULL) +
  theme_pub(base_size = 12, legend = "bottom")
save_fig(p2, file.path(ASSETS, "02_baseline_vs_survshap_dumbbell"), width = 7.4, height = 4.6)

## (3) 全局排列重要性 lollipop(基线本体的精修版,替代默认条形图)
gi <- global_imp_df; gi$feature <- factor(gi$feature, levels = gi$feature)
p3 <- ggplot(gi, aes(importance, feature)) +
  geom_segment(aes(x = 0, xend = importance, y = feature, yend = feature),
               colour = "grey75", linewidth = 0.9) +
  geom_point(aes(colour = feature), size = 5) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(title = "Baseline: global permutation importance (time-integrated)",
       subtitle = "A single scalar per feature — the view SurvSHAP(t) improves on",
       x = "Time-integrated permutation loss increase", y = NULL) +
  theme_pub(base_size = 12)
save_fig(p3, file.path(ASSETS, "03_baseline_permutation_lollipop"), width = 6.6, height = 4.2)

## (4) 单样本 SurvSHAP(t) profile(焦点样本,带正负贡献)
p4 <- ggplot(single_shap, aes(time, shap, colour = feature)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pal, name = "Feature") +
  labs(title = "Single-observation SurvSHAP(t) profile",
       subtitle = sprintf("Focus sample #%d · signed contribution to survival over time", foc),
       x = "Survival time t", y = "SHAP(t)  (+ raises risk / - lowers)") +
  theme_pub(base_size = 12)
save_fig(p4, file.path(ASSETS, "04_single_obs_survshap_profile"), width = 7.2, height = 5)

## (5) SurvSHAP(t) 热图:特征 × 时间(跨样本平均 |SHAP|),viridis 连续量
hm <- shap_curve
hm$feature <- factor(hm$feature, levels = rev(shap_global_df$feature))
p5 <- ggplot(hm, aes(time, feature, fill = mean_abs_shap)) +
  geom_tile() +
  scale_fill_cont(option = "D", name = "Mean |SHAP(t)|") +
  labs(title = "SurvSHAP(t) heatmap: feature x survival time",
       subtitle = "Rows ordered by time-integrated importance; bright = strong effect at that time",
       x = "Survival time t", y = NULL) +
  theme_pub(base_size = 12)
save_fig(p5, file.path(ASSETS, "05_survshap_t_heatmap"), width = 7.4, height = 4.4)

## (6) 单样本 SurvLIME 系数 lollipop(局部线性代理解释)
if (!is.null(lime_df)) {
  ld <- lime_df[order(lime_df$coef), ]; ld$feature <- factor(ld$feature, levels = ld$feature)
  ld$dir <- ifelse(ld$coef >= 0, "raises risk", "lowers risk")
  p6 <- ggplot(ld, aes(coef, feature)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = coef, y = feature, yend = feature, colour = dir),
                 linewidth = 1) +
    geom_point(aes(colour = dir), size = 4.5) +
    scale_colour_manual(values = c("raises risk" = "#E64B35", "lowers risk" = "#4DBBD5"),
                        name = NULL) +
    labs(title = "Single-observation SurvLIME explanation",
         subtitle = sprintf("Local linear surrogate coefficients · focus sample #%d", foc),
         x = "SurvLIME coefficient", y = NULL) +
    theme_pub(base_size = 12, legend = "bottom")
  save_fig(p6, file.path(ASSETS, "06_single_obs_survlime_lollipop"), width = 6.6, height = 4.2)
}

## ---- 收尾:依赖快照(铁律6) ---------------------------------------------
cat(sprintf("\n完成。模型=%s · C-index=%.3f\n", toupper(args$model), cidx))
cat("结果表:", normalizePath(args$outdir), "\n图(PDF+PNG):", normalizePath(ASSETS), "\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
