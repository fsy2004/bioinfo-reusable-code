# =============================================================================
# 编号       : 057
# 脚本名     : TCGA 预后风险模型可视化 (turnkey + 顶刊图)
# 分类       : 12_TCGA_肿瘤预后生存_仅参考
# 用途       : 基于风险评分文件绘制风险分布、生存状态、基因热图、KM 生存曲线、
#              1/3/5 年时间依赖 ROC —— 预后签名的标准五件套。
# 方法/包    : survival + survminer + timeROC + ComplexHeatmap;主题 theme_pub.R
# 结果图     : Risk_distribution;Survival_status;Risk_heatmap;KM_curve;timeROC
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 057_prognostic_risk_model.R
# 运行(自己): Rscript 057_prognostic_risk_model.R --input data/risk.csv
# 输入规格 : CSV,必含列 futime(天), fustat(0/1), riskScore, risk(low/high);
#            其余数值列视为风险基因表达。
# 整理日期 : 2026-06-23(turnkey 重构;base/ggplot 图统一 theme_pub)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(survival); library(survminer); library(timeROC); library(ggplot2); library(ComplexHeatmap); library(circlize) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "risk.csv"), outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)
COL <- c(low = "#0072B5", high = "#BC3C29")

cat("Step 1/5: 读取风险数据...\n")
rd <- read_table_smart(args$input)
rd$futime <- rd$futime / 365
rd <- rd[order(rd$riskScore), ]; rd$rank <- seq_len(nrow(rd))
rd$risk <- factor(rd$risk, levels = c("low", "high"))
n_low <- sum(rd$risk == "low"); med <- median(rd$riskScore)
cat("  ", nrow(rd), "患者 · low", n_low, "/ high", nrow(rd) - n_low, "\n")

# ---- Step 2. 风险评分分布 ----
cat("Step 2/5: 风险分布 + 生存状态...\n")
p_risk <- ggplot(rd, aes(rank, riskScore, colour = risk)) + geom_point(size = 1.6, alpha = .85) +
  scale_colour_manual(values = COL, name = "Risk") +
  geom_vline(xintercept = n_low, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = med, linetype = "dashed", colour = "grey50") +
  labs(title = "Risk score distribution", x = "Patients (ranked)", y = "Risk score") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_risk, file.path(ASSETS, "Risk_distribution"), 9, 3.5); save_fig(p_risk, file.path(args$outdir, "Risk_distribution"), 9, 3.5)

rd$Status <- factor(rd$fustat, levels = c(0, 1), labels = c("Alive", "Dead"))
p_stat <- ggplot(rd, aes(rank, futime, colour = risk, shape = Status)) + geom_point(size = 2, alpha = .85) +
  scale_colour_manual(values = COL, name = "Risk") + scale_shape_manual(values = c(Alive = 16, Dead = 17)) +
  geom_vline(xintercept = n_low, linetype = "dashed", colour = "grey50") +
  labs(title = "Survival status", x = "Patients (ranked)", y = "Survival time (years)") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_stat, file.path(ASSETS, "Survival_status"), 9, 3.5); save_fig(p_stat, file.path(args$outdir, "Survival_status"), 9, 3.5)

# ---- Step 3. 风险基因热图 ----
cat("Step 3/5: 基因热图...\n")
non_expr <- c("riskScore", "futime", "fustat", "rank", "id", "Status", "risk")
gcols <- names(rd)[sapply(rd, is.numeric) & !names(rd) %in% non_expr]
em <- t(scale(as.matrix(rd[, gcols, drop = FALSE])))
ha <- HeatmapAnnotation(Risk = rd$risk, col = list(Risk = COL))
ht <- Heatmap(em, name = "Z-score", top_annotation = ha, col = colorRamp2(c(-2, 0, 2), c("#0072B5", "white", "#BC3C29")),
              cluster_columns = FALSE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = "italic"))
for (dest in c(file.path(ASSETS, "Risk_heatmap"), file.path(args$outdir, "Risk_heatmap"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 9, height = 3); draw(ht); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 9, height = 3, units = "in", res = 300); draw(ht); dev.off()
}

# ---- Step 4. KM 曲线 ----
cat("Step 4/5: KM 生存曲线...\n")
fit <- survfit(Surv(futime, fustat) ~ risk, data = rd)
cox <- summary(coxph(Surv(futime, fustat) ~ risk, data = rd))
hr <- sprintf("HR = %.2f (%.2f-%.2f)\n%s", cox$conf.int[1], cox$conf.int[3], cox$conf.int[4],
              if (cox$coefficients[1, 5] < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", cox$coefficients[1, 5])))
km <- ggsurvplot(fit, data = rd, conf.int = TRUE, pval = hr, pval.size = 4.5, risk.table = TRUE,
                 legend.labs = c("Low risk", "High risk"), legend.title = "Risk", xlab = "Time (years)",
                 palette = unname(COL), risk.table.height = .26, ggtheme = theme_pub(base_size = 13, border = TRUE))
for (dest in c(file.path(ASSETS, "KM_curve"), file.path(args$outdir, "KM_curve"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 6.5, height = 6.5, onefile = FALSE); print(km); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 6.5, height = 6.5, units = "in", res = 300); print(km); dev.off()
}

# ---- Step 5. 时间依赖 ROC ----
cat("Step 5/5: 时间依赖 ROC...\n")
tmax <- max(rd$futime); times <- c(1, 3, 5)[c(1, 3, 5) < tmax]
tr <- timeROC(T = rd$futime, delta = rd$fustat, marker = rd$riskScore, cause = 1, weighting = "marginal", times = times, ROC = TRUE)
rocdf <- do.call(rbind, lapply(seq_along(times), function(i) data.frame(FPR = tr$FP[, i], TPR = tr$TP[, i],
  Time = sprintf("%d-year (AUC=%.3f)", times[i], tr$AUC[i]))))
p_roc <- ggplot(rocdf, aes(FPR, TPR, colour = Time)) + geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = pal_pub(length(times), "npg"), name = NULL) +
  labs(title = "Time-dependent ROC", x = "1 - Specificity", y = "Sensitivity") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.99, .02), legend.justification = c(1, 0))
save_fig(p_roc, file.path(ASSETS, "timeROC"), 5.5, 5.5); save_fig(p_roc, file.path(args$outdir, "timeROC"), 5.5, 5.5)
write.csv(data.frame(Time = times, AUC = tr$AUC[!is.na(tr$AUC)]), file.path(args$outdir, "roc_auc.csv"), row.names = FALSE)
cat("完成。预后五件套见", normalizePath(ASSETS), "\n")
