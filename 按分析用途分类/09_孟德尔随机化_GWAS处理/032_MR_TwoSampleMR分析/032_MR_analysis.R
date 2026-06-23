# =============================================================================
# 编号       : 032
# 脚本名     : 孟德尔随机化(MR)分析 (turnkey + 顶刊图)
# 分类       : 09_孟德尔随机化_GWAS处理
# 用途       : 从 harmonized 工具变量数据做 MR(IVW / Egger / 加权中位数),输出散点、
#              森林、漏斗、留一法图与结果表。自包含实现,不依赖 TwoSampleMR/LD 服务。
# 方法/包    : base R(IVW/Egger/WM)+ ggplot2;主题 theme_pub.R
# 结果图     : MR_scatter;MR_forest;MR_funnel;MR_leaveoneout
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 032_MR_analysis.R
# 运行(自己): Rscript 032_MR_analysis.R --input data/harmonized.csv
# 输入规格 : harmonized CSV,需含列 SNP, beta.exposure, se.exposure, beta.outcome, se.outcome
#            (即 TwoSampleMR::harmonise_data 的输出格式)。
# 整理日期 : 2026-06-23(turnkey 重构;核心 MR 自包含实现,替代重型 TwoSampleMR 依赖)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages(library(ggplot2)))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "harmonized_data.csv"), outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/3: 读取 harmonized 数据 + MR 估计...\n")
d <- read_table_smart(args$input)
bx <- d$beta.exposure; by <- d$beta.outcome; sx <- d$se.exposure; sy <- d$se.outcome
w <- 1 / sy^2
# IVW(固定效应,过原点加权回归)
b_ivw <- sum(w * bx * by) / sum(w * bx^2); se_ivw <- sqrt(1 / sum(w * bx^2))
# MR-Egger(带截距加权回归)
eg <- summary(lm(by ~ bx, weights = w)); b_egger <- eg$coefficients[2, 1]; se_egger <- eg$coefficients[2, 2]
intercept <- eg$coefficients[1, 1]; p_pleio <- eg$coefficients[1, 4]
# 加权中位数
wald <- by / bx; vw <- (sy^2) / (bx^2)
o <- order(wald); ws <- (1 / vw)[o] / sum(1 / vw); cw <- cumsum(ws); ww <- wald[o]
b_wm <- approx(cw, ww, xout = 0.5, ties = "ordered")$y
methods <- data.frame(Method = c("IVW", "MR-Egger", "Weighted median"),
                      b = c(b_ivw, b_egger, b_wm), se = c(se_ivw, se_egger, NA))
methods$OR <- exp(methods$b); methods$p <- 2 * pnorm(-abs(methods$b / methods$se))
write.csv(methods, file.path(args$outdir, "MR_estimates.csv"), row.names = FALSE)
write.csv(data.frame(test = "Egger intercept", intercept = intercept, p = p_pleio), file.path(args$outdir, "MR_pleiotropy.csv"), row.names = FALSE)
cat(sprintf("  IVW beta=%.3f (p=%.1e) · Egger beta=%.3f · WM=%.3f · 多效性 p=%.2f\n", b_ivw, methods$p[1], b_egger, b_wm, p_pleio))

cat("Step 2/3: 散点图 + 漏斗图...\n")
dd <- data.frame(bx, by, sx, sy)
p_sc <- ggplot(dd, aes(bx, by)) +
  geom_errorbar(aes(ymin = by - sy, ymax = by + sy), colour = "grey80", width = 0) +
  geom_errorbarh(aes(xmin = bx - sx, xmax = bx + sx), colour = "grey80", height = 0) +
  geom_point(size = 2, colour = "#3C5488") +
  geom_abline(aes(slope = b_ivw, intercept = 0, colour = "IVW"), linewidth = 0.9) +
  geom_abline(aes(slope = b_egger, intercept = intercept, colour = "MR-Egger"), linewidth = 0.9) +
  geom_abline(aes(slope = b_wm, intercept = 0, colour = "Weighted median"), linewidth = 0.9, linetype = "dashed") +
  scale_colour_manual(values = pal_pub(3, "npg"), name = "Method") +
  labs(title = "MR scatter", x = "SNP effect on exposure", y = "SNP effect on outcome") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.02, .98), legend.justification = c(0, 1))
save_fig(p_sc, file.path(ASSETS, "MR_scatter"), 6.2, 5.5); save_fig(p_sc, file.path(args$outdir, "MR_scatter"), 6.2, 5.5)

dd$wald <- by / bx; dd$inv_se <- abs(bx) / sy
p_fn <- ggplot(dd, aes(wald, inv_se)) + geom_point(size = 2, colour = "#3C5488") +
  geom_vline(xintercept = b_ivw, colour = "#E64B35", linewidth = 0.8) +
  labs(title = "Funnel plot", x = "Wald ratio (per SNP)", y = "1 / SE (instrument strength)") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_fn, file.path(ASSETS, "MR_funnel"), 6, 5); save_fig(p_fn, file.path(args$outdir, "MR_funnel"), 6, 5)

cat("Step 3/3: 森林图 + 留一法...\n")
# 单 SNP Wald + 合并估计
snp <- data.frame(label = d$SNP, b = wald, se = sqrt(vw))
snp <- snp[order(snp$b), ]
comb <- data.frame(label = c("IVW", "MR-Egger"), b = c(b_ivw, b_egger), se = c(se_ivw, se_egger))
fp <- rbind(data.frame(snp, type = "SNP"), data.frame(comb, type = "Combined"))
fp$label <- factor(fp$label, levels = c(rev(snp$label), "MR-Egger", "IVW"))
p_fo <- ggplot(fp, aes(b, label, colour = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.3) +
  geom_point(aes(size = type)) +
  scale_colour_manual(values = c(SNP = "grey50", Combined = "#E64B35"), guide = "none") +
  scale_size_manual(values = c(SNP = 1.6, Combined = 3), guide = "none") +
  labs(title = "Forest plot (per-SNP & combined)", x = "MR effect (beta)", y = NULL) +
  theme_pub(base_size = 10, border = TRUE) + theme(axis.text.y = element_text(size = 7))
save_fig(p_fo, file.path(ASSETS, "MR_forest"), 6, 6.5); save_fig(p_fo, file.path(args$outdir, "MR_forest"), 6, 6.5)

# 留一法
loo <- data.frame(label = d$SNP, b = sapply(seq_len(nrow(d)), function(i) {
  ww <- w[-i]; sum(ww * bx[-i] * by[-i]) / sum(ww * bx[-i]^2) }),
  se = sapply(seq_len(nrow(d)), function(i) { ww <- w[-i]; sqrt(1 / sum(ww * bx[-i]^2)) }))
loo <- rbind(loo, data.frame(label = "All", b = b_ivw, se = se_ivw))
loo$label <- factor(loo$label, levels = rev(loo$label)); loo$all <- loo$label == "All"
p_loo <- ggplot(loo, aes(b, label, colour = all)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.3) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#E64B35"), guide = "none") +
  labs(title = "Leave-one-out", x = "IVW beta (excluding SNP)", y = NULL) +
  theme_pub(base_size = 10, border = TRUE) + theme(axis.text.y = element_text(size = 7))
save_fig(p_loo, file.path(ASSETS, "MR_leaveoneout"), 6, 6.5); save_fig(p_loo, file.path(args$outdir, "MR_leaveoneout"), 6, 6.5)
cat("完成。MR 图/表见", normalizePath(args$outdir), "\n")
