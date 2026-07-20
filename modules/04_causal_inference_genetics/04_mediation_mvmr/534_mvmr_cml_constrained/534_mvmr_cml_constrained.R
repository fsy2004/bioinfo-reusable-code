# =============================================================================
# 编号   : 534
# 脚本名 : 多变量约束极大似然 MVMR(MVMR-cML / MVcML-DP)
# 分类   : 09_mendelian_randomization
# 用途   : 多个暴露 → 1 个结局的多变量 MR,用约束极大似然(constrained maximum
#          likelihood)同时抗【相关 + 非相关】水平多效性:以 BIC 自动选出"无效工具
#          变量"(invalid IVs)并将其多效性效应作为自由参数估计,得到去偏的各暴露
#          直接因果效应(direct effect);DP=data perturbation 给出稳健标准误/置信区间。
# ★诚实基线 : 内置对照 MVMR-IVW(mr_mvivw,不抗多效性)。合成数据故意注入【方向性
#             水平多效性】使 IVW 估计被系统性放大;脚本实测对比 cML 是否更接近真值——
#             不只报 cML 的好看指标,而是把"基线被带偏 vs cML 去偏"摆在同一张图上。
# 依赖   : MendelianRandomization(>=0.10.0, 提供 mr_mvinput/mr_mvcML/mr_mvivw)·
#          ggplot2 · 框架 theme_pub.R
# 运行   : Rscript 534_mvmr_cml_constrained.R                 # 合成示例(脚本内生成)
#          Rscript 534_mvmr_cml_constrained.R --input my_mvmr.csv --n 80000 --outdir results/run1
# 输入   : 一个长表 csv,每行 = 1 个工具 SNP × 1 个暴露的关联,列:
#          SNP, exposure, bx, bxse, outcome, by, byse   (见 README ①;脚本会 reshape 成矩阵)
# =============================================================================

## ---- 框架定位 + 依赖 --------------------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(MendelianRandomization); library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input   = file.path(DDAT, "mvmr_summary.csv"),
  outdir  = file.path(SCRIPT_DIR, "results"),
  n       = 80000,          # 结局 GWAS 样本量(cML 需要)
  num_pert = 100,           # data-perturbation 次数(>=100 推荐;调小更快)
  kmax    = 20))            # BIC 搜索的最大无效工具数(K_vec = 0:kmax)
for (k in c("n","num_pert","kmax")) args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 0. 合成 MVMR summary(若无输入文件) -----------------------------------
# synthetic, for demo only —— 3 暴露(LDL/HDL/TG)→ 主结局 CAD;另造 2 个次结局供热图。
# 设计要点(让诚实基线可被检验):
#   * 真直接效应 theta_true = (LDL=+0.40, HDL=-0.25, TG=0)      —— TG 为零效应(阴性对照)
#   * 在前 ~24/80 个 SNP 注入【方向性(全为正)水平多效性】→ 系统性抬高 by
#     → MVMR-IVW 把它误当作因果 → 估计被放大(尤其 LDL);MVMR-cML 应识别这些 SNP
#       为 invalid 并去偏,估计回到真值附近。
EXPO <- c("LDL", "HDL", "TG"); THETA_TRUE <- c(LDL = 0.40, HDL = -0.25, TG = 0.0)
make_mvmr <- function(p = 80, theta = THETA_TRUE, n_invalid = 24, outcome = "CAD",
                      plei_mean = 0.06, seed = 1) {
  set.seed(seed); k <- length(theta)
  # 先抽数值矩阵(决定 cML 的数值稳定性),再抽 rsid 标签 —— 顺序固定保证可复现
  bx   <- matrix(abs(rnorm(p * k, 0, 0.10)), p, k, dimnames = list(NULL, names(theta)))
  bxse <- matrix(runif(p * k, 0.008, 0.018), p, k, dimnames = list(NULL, names(theta)))
  plei <- numeric(p); inv <- sample.int(p, n_invalid)
  plei[inv] <- abs(rnorm(n_invalid, plei_mean, 0.03))     # 方向性(正)多效性
  by   <- as.vector(bx %*% theta) + plei + rnorm(p, 0, 0.008)
  byse <- runif(p, 0.008, 0.018)
  snp  <- sprintf("rs%07d", sample(1e6:9e6, p))
  rownames(bx) <- rownames(bxse) <- snp
  long <- do.call(rbind, lapply(seq_len(k), function(j) data.frame(
    SNP = snp, exposure = names(theta)[j], bx = bx[, j], bxse = bxse[, j],
    outcome = outcome, by = by, byse = byse, stringsAsFactors = FALSE)))
  long
}
if (!file.exists(args$input)) {
  # 主结局 CAD(强多效性,IVW 应明显被带偏)+ 两个次结局供"多暴露×多结局"热图
  d_cad <- make_mvmr(p = 80, theta = THETA_TRUE,             n_invalid = 24, outcome = "CAD",    plei_mean = 0.06, seed = 1)
  d_t2d <- make_mvmr(p = 80, theta = c(LDL=0.10,HDL=-0.05,TG=0.30), n_invalid = 12, outcome = "T2D",    plei_mean = 0.03, seed = 2)
  d_str <- make_mvmr(p = 80, theta = c(LDL=0.25,HDL=-0.30,TG=0.05), n_invalid = 18, outcome = "Stroke", plei_mean = 0.05, seed = 3)
  write.csv(rbind(d_cad, d_t2d, d_str), args$input, row.names = FALSE)
  cat("[gen] 合成 MVMR summary: 3 暴露 × 3 结局, 每结局 80 SNP (synthetic, demo only)\n")
}

## ---- 1. 读入 + reshape 成矩阵 ----------------------------------------------
cat("Step 1: 读 MVMR summary → reshape 成 (SNP × 暴露) 矩阵...\n")
raw <- read.csv(args$input, stringsAsFactors = FALSE)
stopifnot(all(c("SNP","exposure","bx","bxse","outcome","by","byse") %in% names(raw)))
PRIMARY <- raw$outcome[1]                      # 第一个出现的结局 = 主结局(forest/dumbbell 对象)
expo_levels <- unique(raw$exposure)
cat(sprintf("  暴露 = %s ; 结局 = %s ; 主结局 = %s\n",
            paste(expo_levels, collapse = "/"), paste(unique(raw$outcome), collapse = "/"), PRIMARY))

# 把一个 outcome 的长表转成 mr_mvinput(bx/bxse 矩阵 + by/byse 向量)
build_input <- function(df) {
  snps <- unique(df$SNP)
  # ★每个暴露须在【该暴露子表内部】用 match 取值:df$bx[match(snps, df$SNP[...])]
  #   会拿子表内位置去索全表 df$bx → 取错行(三列全塌成 LDL 块)。改为先取子表再 match。
  bx   <- sapply(expo_levels, function(e) { s <- df[df$exposure == e, ]; s$bx[match(snps, s$SNP)] })
  bxse <- sapply(expo_levels, function(e) { s <- df[df$exposure == e, ]; s$bxse[match(snps, s$SNP)] })
  # by/byse 与暴露无关(同一结局),取任一暴露子表即可
  sub  <- df[df$exposure == expo_levels[1], ]
  by   <- sub$by[match(snps, sub$SNP)]; byse <- sub$byse[match(snps, sub$SNP)]
  mr_mvinput(bx = as.matrix(bx), bxse = as.matrix(bxse), by = by, byse = byse,
             exposure = expo_levels, outcome = df$outcome[1], snps = snps)
}

## ---- 2. 主结局:MVMR-cML(DP) + 诚实基线 MVMR-IVW ---------------------------
cat("Step 2: 主结局 MVMR-cML-DP(抗多效性) vs MVMR-IVW(诚实基线,不抗多效性)...\n")
inp <- build_input(raw[raw$outcome == PRIMARY, ])

cml <- mr_mvcML(inp, n = args$n, DP = TRUE, num_pert = args$num_pert,
                K_vec = 0:args$kmax, seed = 42)
ivw <- mr_mvivw(inp)

cmp <- data.frame(
  exposure = expo_levels,
  theta_true = THETA_TRUE[expo_levels],
  cML_est = cml@Estimate, cML_lo = cml@CILower, cML_hi = cml@CIUpper, cML_p = cml@Pvalue,
  IVW_est = ivw@Estimate, IVW_lo = ivw@CILower, IVW_hi = ivw@CIUpper, IVW_p = ivw@Pvalue,
  row.names = NULL)
cmp$cML_absbias <- abs(cmp$cML_est - cmp$theta_true)
cmp$IVW_absbias <- abs(cmp$IVW_est - cmp$theta_true)
write.csv(cmp, file.path(args$outdir, "mvmr_cML_vs_IVW.csv"), row.names = FALSE)

n_inv <- length(cml@BIC_invalid)
# DP 模式下 @K_hat 是【每次 perturbation 各一个】的向量 → 取中位数作代表;
# @SNPs 本身就是"参与计算的 SNP 数"(标量),勿再 length()。
k_hat_typ <- stats::median(cml@K_hat)
writeLines(c(sprintf("主结局: %s", PRIMARY),
             sprintf("BIC 选出无效工具数 K_hat(每次扰动中位数) = %g | DP 收敛次数 = %d/%d",
                     k_hat_typ, cml@eff_DP_B, args$num_pert),
             sprintf("最终标记无效工具 = %d / %d SNP", n_inv, cml@SNPs),
             sprintf("无效工具索引: %s", paste(sort(cml@BIC_invalid), collapse=", "))),
           file.path(args$outdir, "mvmr_cML_invalid_IVs.txt"))

cat("  —— 直接效应估计(真值在括号内)——\n")
for (i in seq_len(nrow(cmp))) cat(sprintf(
  "  %-4s  真%+.2f | cML %+.3f (bias %.3f, p=%.1e) | IVW %+.3f (bias %.3f, p=%.1e)\n",
  cmp$exposure[i], cmp$theta_true[i], cmp$cML_est[i], cmp$cML_absbias[i], cmp$cML_p[i],
  cmp$IVW_est[i], cmp$IVW_absbias[i], cmp$IVW_p[i]))
cat(sprintf("  cML 识别 invalid IV = %d 个 | 平均 |bias|: cML %.3f vs IVW %.3f → %s\n",
  n_inv, mean(cmp$cML_absbias), mean(cmp$IVW_absbias),
  ifelse(mean(cmp$cML_absbias) < mean(cmp$IVW_absbias), "cML 更接近真值(基线被多效性带偏)", "本次基线未明显偏")))

## ---- 3. 多暴露 × 多结局:每结局都跑 cML(供热图) ---------------------------
cat("Step 3: 对每个结局跑 MVMR-cML(多暴露×多结局 beta 矩阵)...\n")
beta_list <- lapply(unique(raw$outcome), function(oc) {
  r <- mr_mvcML(build_input(raw[raw$outcome == oc, ]), n = args$n, DP = TRUE,
                num_pert = args$num_pert, K_vec = 0:args$kmax, seed = 42)
  data.frame(outcome = oc, exposure = expo_levels, beta = r@Estimate, pval = r@Pvalue)
})
beta_df <- do.call(rbind, beta_list)
write.csv(beta_df, file.path(args$outdir, "mvmr_cML_beta_matrix.csv"), row.names = FALSE)

## ============================================================================
## 4. 顶刊级图(全部非平凡条形:forest / dumbbell / heatmap)
## ============================================================================
cat("Step 4: 出图(forest / dumbbell / heatmap)...\n")
COL <- pal_pub(name = "npg")  # [1]红 [2]蓝 [3]绿 ...

## 4a. cML 各暴露直接效应 forest(点 + 95% CI;参考真值竖线)----------------
cmp$lab <- factor(cmp$exposure, levels = rev(cmp$exposure))
cmp$sig <- ifelse(cmp$cML_p < 0.05, "p < 0.05", "n.s.")
p_forest <- ggplot(cmp, aes(x = cML_est, y = lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  geom_errorbar(aes(xmin = cML_lo, xmax = cML_hi), orientation = "y", width = 0.18, linewidth = 0.7, colour = "grey25") +
  geom_point(aes(colour = sig), size = 4.2) +
  geom_point(aes(x = theta_true), shape = 124, size = 6, colour = COL[3], stroke = 1.2) +  # 真值刻度
  scale_colour_manual(values = c("p < 0.05" = COL[1], "n.s." = "grey60"), name = NULL) +
  labs(title = sprintf("MVMR-cML direct effects on %s", PRIMARY),
       subtitle = "Point = cML estimate; bar = 95% CI (data-perturbation); green tick = ground truth",
       x = "Direct causal effect (per-SD)", y = NULL) +
  theme_pub(legend = "top")
save_fig(p_forest, file.path(ASSETS, "fig1_mvcml_forest"), width = 7, height = 3.6)

## 4b. cML vs IVW dumbbell(哑铃图:同一暴露两法估计差距 = 多效性偏倚)------
dmb <- cmp[, c("exposure","cML_est","IVW_est","theta_true")]
dmb$exposure <- factor(dmb$exposure, levels = rev(cmp$exposure))
p_dumb <- ggplot(dmb, aes(y = exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80", linewidth = 0.4) +
  geom_segment(aes(x = IVW_est, xend = cML_est, yend = exposure),
               colour = "grey60", linewidth = 1.1) +
  geom_point(aes(x = IVW_est, colour = "MVMR-IVW (baseline)"), size = 4.3) +
  geom_point(aes(x = cML_est, colour = "MVMR-cML"), size = 4.3) +
  geom_point(aes(x = theta_true), shape = 124, size = 6, colour = COL[3], stroke = 1.2) +
  scale_colour_manual(values = c("MVMR-IVW (baseline)" = COL[2], "MVMR-cML" = COL[1]), name = NULL) +
  labs(title = sprintf("Pleiotropy bias: cML corrects IVW on %s", PRIMARY),
       subtitle = "Grey link = shift from IVW to cML; green tick = ground truth (cML lands closer)",
       x = "Direct causal effect (per-SD)", y = NULL) +
  theme_pub(legend = "top")
save_fig(p_dumb, file.path(ASSETS, "fig2_cml_vs_ivw_dumbbell"), width = 7, height = 3.6)

## 4c. 多暴露 × 多结局 cML beta 热图(发散 RdBu;显著者描边 + 数值)--------
beta_df$exposure <- factor(beta_df$exposure, levels = expo_levels)
beta_df$outcome  <- factor(beta_df$outcome,  levels = unique(raw$outcome))
beta_df$sig <- beta_df$pval < 0.05
p_heat <- ggplot(beta_df, aes(x = outcome, y = exposure, fill = beta)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_tile(data = subset(beta_df, sig), colour = "black", linewidth = 0.9, fill = NA) +  # 显著描黑边
  geom_text(aes(label = sprintf("%.2f", beta)), size = 3.6, colour = "grey10") +
  scale_fill_diverge(midpoint = 0, name = "cML\nbeta") +
  labs(title = "MVMR-cML direct effects: exposures x outcomes",
       subtitle = "Black border = p < 0.05; colour = signed direct effect (RdBu)",
       x = "Outcome", y = "Exposure") +
  coord_equal() +
  theme_pub()
save_fig(p_heat, file.path(ASSETS, "fig3_exposure_outcome_heatmap"), width = 6.2, height = 3.8)

cat("完成。结果表见", normalizePath(args$outdir), ";图见 assets/\n")

## ---- 依赖版本快照(铁律6)---------------------------------------------------
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
