# =============================================================================
# 编号   : 535
# 脚本名 : 偏差校正估计方程 MR(MRBEE)— 去测量误差偏倚 + 抗水平多效性
# 分类   : 09_mendelian_randomization
# 用途   : 用 MRBEE(Bias-corrected Estimating Equation MR, Nat Commun 2024)做
#          单/多变量 MR:在估计方程中显式扣除工具变量测量误差(exposure GWAS 的
#          beta 本身带抽样误差)带来的偏倚,并用 IMRP 迭代多效性剔除(delta≠0 的
#          SNP 判为水平多效性离群)抗 pleiotropy。尤其适合 cis 区(工具少、常偏弱)
#          与多变量(MVMR)场景。
# ★诚实基线 : 同一套合成数据上并行跑 naive IVW(不校正测量误差偏倚)。在弱工具下
#          IVW 因 regression dilution 被系统性衰减(attenuation toward null),而
#          MRBEE 把估计拉回真值附近 —— 图中直接对照,不只报 MRBEE 的好看数字。
# 依赖   : MRBEE(必需,真包实跑;exported: MRBEE.IMRP / MRBEE.IMRP.UV / errorCov)
#          ggplot2(经 theme_pub.R)。★MRBEEX(多变量正则化扩展)装不上 → 本模块
#          只用基础 MRBEE 的 MV 函数 MRBEE.IMRP,不依赖 MRBEEX。
# 运行   : Rscript 535_mrbee_cis_mr.R                       # 零改动跑合成示例
#          Rscript 535_mrbee_cis_mr.R --uv exposures_uv.csv --mv exposures_mv.csv --outdir results/run1
# 输入   : 见 README ①。两张 summary 表:
#          (A) 单变量多暴露表 uv: 每行 = 1 个 (exposure, SNP) 的 GWAS summary
#              列 = exposure, SNP, bx, bxse, by, byse  (by/byse = 同一 outcome)
#          (B) 多变量表 mv: 每行 = 1 个 SNP,列 = SNP, bx_<E1>, bxse_<E1>, ...,
#              by, byse  (宽表,多列暴露)
# =============================================================================

## ---- 定位框架 theme_pub.R 并加载 -------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(MRBEE); library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  uv     = file.path(DDAT, "exposures_uv.csv"),   # 单变量多暴露长表
  mv     = file.path(DDAT, "exposures_mv.csv"),   # 多变量宽表
  outdir = file.path(SCRIPT_DIR, "results"),
  rho_xy = 0.05))                                 # 暴露/结局 GWAS 样本重叠引入的测量误差相关(无重叠=0)
args$rho_xy <- as.numeric(args$rho_xy)
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. 合成示例数据(synthetic, for demo only)
#    设计要点(让诚实基线有意义):
#    - 工具偏弱:exposure beta 自带可观测的抽样误差(bxse 相对 bx 不小)→
#      naive IVW 会被 regression dilution 系统性衰减;MRBEE 校正测量误差应拉回真值。
#    - 注入水平多效性:部分 SNP 的 by 叠加与 bx 无关的额外效应 → IMRP 应把它们
#      的 delta 标为非 0(离群),从而抗多效性。
# =============================================================================
gen_uv_exposure <- function(name, true_b, n, frac_plei = 0.12,
                            bx_sd = 0.06, bxse_rng = c(0.030, 0.045),
                            byse_rng = c(0.012, 0.020), plei_sd = 0.06) {
  bx_true <- rnorm(n, 0, bx_sd)
  bxse <- runif(n, bxse_rng[1], bxse_rng[2])
  bx   <- bx_true + rnorm(n, 0, 1) * bxse                  # 观测 = 真值 + 测量误差(关键)
  byse <- runif(n, byse_rng[1], byse_rng[2])
  by_true <- true_b * bx_true
  plei <- rep(0, n); idx <- sample(n, max(1, round(frac_plei * n)))
  plei[idx] <- rnorm(length(idx), 0, plei_sd)              # 水平多效性:与 bx 无关的额外效应
  by <- by_true + plei + rnorm(n, 0, 1) * byse
  data.frame(exposure = name, SNP = sprintf("%s_rs%05d", name, seq_len(n)),
             bx = bx, bxse = bxse, by = by, byse = byse,
             true_beta = true_b, is_pleio = plei != 0)
}

if (!file.exists(args$uv) || !file.exists(args$mv)) {
  cat("[gen] 未找到输入 → 生成合成 MR summary(synthetic demo only)\n")
  # (A) 单变量:5 个暴露,真因果各异(含 1 个零效应作阴性对照)
  uv_spec <- list(ExpA = 0.50, ExpB = 0.35, ExpC = -0.40, ExpD = 0.20, ExpE = 0.00)
  uv <- do.call(rbind, Map(function(nm, b) gen_uv_exposure(nm, b, n = 70),
                           names(uv_spec), unlist(uv_spec)))
  write.csv(uv[, c("exposure","SNP","bx","bxse","by","byse")], args$uv, row.names = FALSE)
  write.csv(uv[, c("exposure","SNP","true_beta","is_pleio")],
            file.path(DDAT, "exposures_uv_truth.csv"), row.names = FALSE)

  # (B) 多变量:3 个相关暴露共享一组 SNP(真效应 0.45 / -0.25 / 0.10)
  nmv <- 130; mv_true <- c(MV1 = 0.45, MV2 = -0.25, MV3 = 0.10)
  L <- matrix(c(1,0.3,0.2, 0.3,1,0.25, 0.2,0.25,1), 3, 3)        # 暴露间相关
  Z <- matrix(rnorm(nmv * 3), nmv, 3) %*% chol(L)
  bX_true <- Z * 0.06
  bXse <- matrix(runif(nmv * 3, 0.030, 0.045), nmv, 3)
  bX   <- bX_true + matrix(rnorm(nmv * 3), nmv, 3) * bXse        # 观测 + 测量误差
  byse <- runif(nmv, 0.012, 0.020)
  plei <- rep(0, nmv); pidx <- sample(nmv, round(0.12 * nmv))
  plei[pidx] <- rnorm(length(pidx), 0, 0.06)
  by <- as.numeric(bX_true %*% mv_true) + plei + rnorm(nmv, 0, 1) * byse
  mvdf <- data.frame(SNP = sprintf("MV_rs%05d", seq_len(nmv)))
  for (j in 1:3) { mvdf[[paste0("bx_", names(mv_true)[j])]] <- bX[, j]
                   mvdf[[paste0("bxse_", names(mv_true)[j])]] <- bXse[, j] }
  mvdf$by <- by; mvdf$byse <- byse
  write.csv(mvdf, args$mv, row.names = FALSE)
  attr_mv_true <- mv_true
}

# =============================================================================
# 1. naive IVW(诚实基线;不校正测量误差偏倚)
#    IVW = 过原点的逆方差加权回归(weights = 1/byse^2)。SE 取一阶 delta 法。
# =============================================================================
naive_ivw <- function(bx, by, byse) {
  w <- 1 / byse^2
  beta <- sum(w * bx * by) / sum(w * bx^2)
  se   <- sqrt(1 / sum(w * bx^2))               # 固定效应 IVW 标准误
  c(beta = beta, se = se)
}

# =============================================================================
# 2. 单变量 MRBEE(逐暴露)+ naive IVW 对照
#    MRBEE.IMRP.UV(by,bx,byse,bxse,Rxy): Rxy = 2x2 测量误差相关阵 [exposure, outcome]
#    返回: theta(因果估计) / vartheta(方差) / delta(每 SNP 多效性,≠0=离群)
# =============================================================================
cat("Step 1: 单变量 MRBEE vs naive IVW(逐暴露)...\n")
uv <- read.csv(args$uv, stringsAsFactors = FALSE)
truth_file <- file.path(DDAT, "exposures_uv_truth.csv")
truth <- if (file.exists(truth_file)) read.csv(truth_file, stringsAsFactors = FALSE) else NULL

Rxy_uv <- matrix(c(1, args$rho_xy, args$rho_xy, 1), 2, 2)
exps <- unique(uv$exposure)
uv_res <- do.call(rbind, lapply(exps, function(e) {
  d <- uv[uv$exposure == e, ]
  fit <- MRBEE.IMRP.UV(by = d$by, bx = d$bx, byse = d$byse, bxse = d$bxse, Rxy = Rxy_uv)
  iv  <- naive_ivw(d$bx, d$by, d$byse)
  tb  <- if (!is.null(truth)) truth$true_beta[truth$exposure == e][1] else NA_real_
  data.frame(exposure = e, n_snp = nrow(d),
             mrbee_beta = as.numeric(fit$theta), mrbee_se = sqrt(as.numeric(fit$vartheta)),
             ivw_beta = iv["beta"], ivw_se = iv["se"],
             n_pleio_flag = sum(fit$delta != 0), true_beta = tb,
             row.names = NULL)
}))
uv_res$mrbee_p <- 2 * pnorm(-abs(uv_res$mrbee_beta / uv_res$mrbee_se))
uv_res$ivw_p   <- 2 * pnorm(-abs(uv_res$ivw_beta   / uv_res$ivw_se))
write.csv(uv_res, file.path(args$outdir, "MRBEE_uv_results.csv"), row.names = FALSE)
if (all(!is.na(uv_res$true_beta))) {
  cat(sprintf("  平均绝对误差 |est-true|: MRBEE=%.3f  naive_IVW=%.3f  (MRBEE 应更小)\n",
      mean(abs(uv_res$mrbee_beta - uv_res$true_beta)),
      mean(abs(uv_res$ivw_beta   - uv_res$true_beta))))
}

# =============================================================================
# 3. 多变量 MRBEE(MVMR)+ 每暴露 naive IVW 对照
#    MRBEE.IMRP(by,bX,byse,bXse,Rxy): Rxy = (k+1)x(k+1),顺序 [exposures..., outcome]
#    返回: theta(向量) / covtheta(协方差阵) / delta(每 SNP 多效性)
# =============================================================================
cat("Step 2: 多变量 MRBEE(MVMR)...\n")
mv <- read.csv(args$mv, stringsAsFactors = FALSE, check.names = FALSE)
bx_cols  <- grep("^bx_",  names(mv), value = TRUE)
bxse_cols <- grep("^bxse_", names(mv), value = TRUE)
mv_names <- sub("^bx_", "", bx_cols)
bX   <- as.matrix(mv[, bx_cols, drop = FALSE]);  colnames(bX) <- mv_names
bXse <- as.matrix(mv[, bxse_cols, drop = FALSE])
k <- length(mv_names)
Rxy_mv <- diag(k + 1); Rxy_mv[k + 1, 1:k] <- Rxy_mv[1:k, k + 1] <- args$rho_xy
fit_mv <- MRBEE.IMRP(by = mv$by, bX = bX, byse = mv$byse, bXse = bXse, Rxy = Rxy_mv)
mv_true <- if (exists("attr_mv_true")) attr_mv_true else setNames(rep(NA_real_, k), mv_names)
mv_se <- sqrt(diag(fit_mv$covtheta))
# 每暴露 naive 单变量 IVW(忽略其它暴露 + 测量误差),作为对照
mv_ivw <- vapply(seq_len(k), function(j) naive_ivw(bX[, j], mv$by, mv$byse)["beta"], numeric(1))
mv_res <- data.frame(exposure = mv_names,
  mrbee_beta = as.numeric(fit_mv$theta), mrbee_se = mv_se,
  mrbee_p = 2 * pnorm(-abs(as.numeric(fit_mv$theta) / mv_se)),
  ivw_beta = mv_ivw, true_beta = as.numeric(mv_true[mv_names]),
  n_pleio_flag = sum(fit_mv$delta != 0), n_snp = nrow(mv))
write.csv(mv_res, file.path(args$outdir, "MRBEE_mv_results.csv"), row.names = FALSE)
cat(sprintf("  MVMR theta = %s ; IMRP 标记多效性 SNP = %d / %d\n",
    paste(sprintf("%s=%.3f", mv_names, fit_mv$theta), collapse = ", "),
    sum(fit_mv$delta != 0), nrow(mv)))

# =============================================================================
# 4. 顶刊级图(禁止平凡条形;用 lollipop / forest / 散点)
# =============================================================================
cat("Step 3: 出图(lollipop / forest / 工具效应散点)...\n")
pal <- pal_pub(name = "npg")
COL_MRBEE <- pal[1]; COL_IVW <- pal[4]

## 图1. 因果效应 lollipop(单变量多暴露,按 MRBEE 估计排序)----------------------
d1 <- uv_res[order(uv_res$mrbee_beta), ]
d1$exposure <- factor(d1$exposure, levels = d1$exposure)
d1$sig <- ifelse(d1$mrbee_p < 0.05, "p < 0.05", "n.s.")
p1 <- ggplot(d1, aes(x = mrbee_beta, y = exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = 0, xend = mrbee_beta, yend = exposure), colour = "grey70", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = mrbee_beta - 1.96 * mrbee_se, xmax = mrbee_beta + 1.96 * mrbee_se),
                 height = 0, colour = "grey45", linewidth = 0.6) +
  geom_point(aes(colour = sig), size = 4.2) +
  { if (all(!is.na(d1$true_beta)))
      geom_point(aes(x = true_beta), shape = 124, size = 4, colour = "black") } +
  scale_colour_manual(values = c("p < 0.05" = COL_MRBEE, "n.s." = "grey60"), name = NULL) +
  labs(title = "MRBEE causal effects across exposures",
       subtitle = "Lollipop = MRBEE theta +/- 95% CI; black tick = true effect (synthetic)",
       x = "Causal effect (theta)", y = NULL) +
  theme_pub(base_size = 12)
save_fig(p1, file.path(ASSETS, "fig1_mrbee_lollipop"), width = 7, height = 4.6)

## 图2. MRBEE vs naive IVW forest(同一暴露两估计并列,展示去衰减)---------------
fr <- rbind(
  data.frame(exposure = uv_res$exposure, method = "MRBEE",
             beta = uv_res$mrbee_beta, se = uv_res$mrbee_se, true_beta = uv_res$true_beta),
  data.frame(exposure = uv_res$exposure, method = "naive IVW",
             beta = uv_res$ivw_beta, se = uv_res$ivw_se, true_beta = uv_res$true_beta))
ord <- uv_res$exposure[order(uv_res$mrbee_beta)]
fr$exposure <- factor(fr$exposure, levels = ord)
fr$method <- factor(fr$method, levels = c("naive IVW", "MRBEE"))
p2 <- ggplot(fr, aes(x = beta, y = exposure, colour = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  { if (all(!is.na(fr$true_beta)))
      geom_point(aes(x = true_beta), shape = 124, size = 5, colour = "black",
                 position = position_dodge(width = 0.6), show.legend = FALSE) } +
  geom_errorbarh(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se),
                 height = 0, linewidth = 0.7, position = position_dodge(width = 0.6)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("MRBEE" = COL_MRBEE, "naive IVW" = COL_IVW), name = NULL) +
  labs(title = "Bias correction: MRBEE vs naive IVW",
       subtitle = "Weak instruments + measurement error -> IVW attenuated toward null; MRBEE restored",
       x = "Causal effect estimate", y = NULL) +
  theme_pub(base_size = 12)
save_fig(p2, file.path(ASSETS, "fig2_mrbee_vs_ivw_forest"), width = 7.4, height = 4.8)

## 图3. 工具效应散点(以最强暴露为例:bx vs by,标注 IMRP 多效性离群 + 两条斜率)---
e_demo <- uv_res$exposure[which.max(abs(uv_res$mrbee_beta))]
dd <- uv[uv$exposure == e_demo, ]
fit_demo <- MRBEE.IMRP.UV(by = dd$by, bx = dd$bx, byse = dd$byse, bxse = dd$bxse, Rxy = Rxy_uv)
dd$pleio <- ifelse(fit_demo$delta != 0, "pleiotropy (removed)", "valid instrument")
iv_demo <- naive_ivw(dd$bx, dd$by, dd$byse)["beta"]
p3 <- ggplot(dd, aes(x = bx, y = by)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_errorbar(aes(ymin = by - byse, ymax = by + byse), colour = "grey80", linewidth = 0.3, width = 0) +
  geom_errorbarh(aes(xmin = bx - bxse, xmax = bx + bxse), colour = "grey80", linewidth = 0.3, height = 0) +
  geom_point(aes(colour = pleio), size = 2.6, alpha = 0.9) +
  geom_abline(aes(slope = as.numeric(fit_demo$theta), intercept = 0, linetype = "MRBEE"),
              colour = COL_MRBEE, linewidth = 1) +
  geom_abline(aes(slope = iv_demo, intercept = 0, linetype = "naive IVW"),
              colour = COL_IVW, linewidth = 1) +
  scale_colour_manual(values = c("valid instrument" = "grey45",
                                 "pleiotropy (removed)" = pal[8]), name = NULL) +
  scale_linetype_manual(values = c("MRBEE" = "solid", "naive IVW" = "longdash"), name = "slope") +
  labs(title = sprintf("Instrument effects: %s -> outcome", e_demo),
       subtitle = "Each point = one SNP; red = IMRP-flagged pleiotropy; lines = MRBEE vs IVW slope",
       x = "SNP-exposure effect (bx)", y = "SNP-outcome effect (by)") +
  theme_pub(base_size = 12)
save_fig(p3, file.path(ASSETS, "fig3_instrument_scatter"), width = 7, height = 5.4)

cat("完成。结果表见", normalizePath(args$outdir), ";图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖版本快照(铁律6)
