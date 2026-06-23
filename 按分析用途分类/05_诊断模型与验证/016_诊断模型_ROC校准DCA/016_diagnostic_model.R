# =============================================================================
# 编号       : 016
# 脚本名     : 诊断模型 ROC / 校准 / DCA / 列线图 (turnkey + 顶刊图)
# 分类       : 05_诊断模型与验证
# 用途       : 基于特征基因构建 logistic 诊断模型,输出列线图、校准曲线、决策曲线
#              (DCA)、ROC(联合+单基因)、OR 森林图、基因差异箱线图。
# 方法/包    : rms(lrm/nomogram/calibrate)+ rmda(decision_curve)+ pROC;绘图 theme_pub.R
# 结果图     : Nomogram;Calibration;DCA;ROC_combined;ROC_genes;OR_forest;Gene_boxplot
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 016_diagnostic_model.R
# 运行(自己): Rscript 016_diagnostic_model.R --input data/expr.csv --genes data/genes.csv
# 可选参数 : --ctrl _con --case _dis --seed 123
# 输入规格 : 表达矩阵 CSV(首列基因,样本名后缀分组)+ 诊断基因列表 CSV(首列基因名)。
# 整理日期 : 2026-06-23(turnkey 重构;rms/rmda 分析逻辑保持,base 图升级 ggplot)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(rms); library(rmda); library(pROC); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input = file.path(SCRIPT_DIR, "example_data", "Sample_Type_Matrix.csv"),
  genes = file.path(SCRIPT_DIR, "example_data", "diagnostic_genes.csv"),
  outdir = file.path(SCRIPT_DIR, "results"), ctrl = "_con", case = "_dis", seed = "123"))
set.seed(as.integer(args$seed))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 数据 + 模型 ----
cat("Step 1/6: 读取数据 + 拟合 logistic 模型...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
gl <- unique(trimws(as.character(read_table_smart(args$genes)[[1]])))
feat <- intersect(gl, rownames(expr)); if (length(feat) < 2) stop("诊断基因过少。")
df <- as.data.frame(t(as.matrix(expr[feat, , drop = FALSE])))
grp <- ifelse(grepl(paste0(args$ctrl, "$"), rownames(df)), 0L,
       ifelse(grepl(paste0(args$case, "$"), rownames(df)), 1L, NA))
if (any(is.na(grp))) stop("样本名后缀需为 ", args$ctrl, " / ", args$case)
df$Group <- grp
dd <- datadist(df); options(datadist = "dd")
form <- as.formula(paste("Group ~", paste(feat, collapse = " + ")))
fit <- lrm(form, data = df, x = TRUE, y = TRUE)
pred <- predict(fit, type = "fitted")
cat("  ", length(feat), "基因 ·", nrow(df), "样本 · 模型 AUC=", sprintf("%.3f", as.numeric(auc(roc(grp, pred, quiet = TRUE)))), "\n")

# ---- Step 2. 列线图(rms base,直接落 pdf/png)----
cat("Step 2/6: 列线图...\n")
nomo <- nomogram(fit, fun = plogis, fun.at = c(.01, .1, .3, .5, .7, .9, .99), lp = FALSE, funlabel = "Disease risk")
for (ext in c("pdf", "png")) {
  if (ext == "pdf") grDevices::cairo_pdf(file.path(ASSETS, "Nomogram.pdf"), width = 10, height = 6)
  else grDevices::png(file.path(ASSETS, "Nomogram.png"), width = 10, height = 6, units = "in", res = 300)
  plot(nomo, cex.axis = 0.8); dev.off()
}

# ---- Step 3. 校准曲线(ggplot)----
cat("Step 3/6: 校准曲线...\n")
cal <- calibrate(fit, method = "boot", B = 200)
cdf <- data.frame(predy = cal[, "predy"], apparent = cal[, "calibrated.orig"], corrected = cal[, "calibrated.corrected"])
p_cal <- ggplot(cdf, aes(predy)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(aes(y = apparent, colour = "Apparent"), linewidth = 0.8) +
  geom_line(aes(y = corrected, colour = "Bias-corrected"), linewidth = 0.9) +
  scale_colour_manual(values = c(Apparent = "#3C5488", `Bias-corrected` = "#E64B35"), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Calibration curve", x = "Predicted probability", y = "Observed probability") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.02, .98), legend.justification = c(0, 1))
save_fig(p_cal, file.path(ASSETS, "Calibration"), 5.5, 5.5); save_fig(p_cal, file.path(args$outdir, "Calibration"), 5.5, 5.5)

# ---- Step 4. DCA 决策曲线(ggplot)----
cat("Step 4/6: 决策曲线 DCA...\n")
dcdf <- df; dcdf$Model <- pred
dca <- decision_curve(Group ~ Model, data = dcdf, thresholds = seq(0, 1, 0.01),
                      family = binomial("logit"), bootstraps = 50)
dd2 <- dca$derived.data
p_dca <- ggplot(dd2[dd2$model != "None", ], aes(thresholds, NB, colour = model)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(All = "grey55", `Group ~ Model` = "#E64B35"),
                      labels = c(All = "Treat all", `Group ~ Model` = "Diagnostic model"), name = NULL) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  coord_cartesian(ylim = c(-0.05, max(dd2$NB, na.rm = TRUE) * 1.05)) +
  labs(title = "Decision curve analysis", x = "Threshold probability", y = "Net benefit") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.98, .98), legend.justification = c(1, 1))
save_fig(p_dca, file.path(ASSETS, "DCA"), 6, 5); save_fig(p_dca, file.path(args$outdir, "DCA"), 6, 5)

# ---- Step 5. ROC(联合 + 单基因)----
cat("Step 5/6: ROC...\n")
roc_model <- roc(grp, pred, quiet = TRUE)
roc_genes <- lapply(feat, function(g) { r <- roc(grp, df[[g]], quiet = TRUE); if (as.numeric(r$auc) < .5) roc(grp, -df[[g]], quiet = TRUE) else r })
names(roc_genes) <- feat
all_roc <- c(list(`Combined model` = roc_model), roc_genes)
aucs <- sapply(all_roc, function(r) as.numeric(r$auc))
p_roc <- pROC::ggroc(all_roc, legacy.axes = TRUE, linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  scale_colour_manual(values = c("#E64B35", pal_pub(length(feat), "lancet")),
                      labels = paste0(names(all_roc), " (", sprintf("%.3f", aucs), ")"), name = "AUC") +
  labs(title = "Diagnostic ROC", x = "1 - Specificity", y = "Sensitivity") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.99, .02), legend.justification = c(1, 0), legend.text = element_text(size = 8))
save_fig(p_roc, file.path(ASSETS, "ROC"), 6, 5.5); save_fig(p_roc, file.path(args$outdir, "ROC"), 6, 5.5)
write.csv(data.frame(predictor = names(aucs), AUC = aucs), file.path(args$outdir, "AUC.csv"), row.names = FALSE)

# ---- Step 6. OR 森林图 + 基因箱线图 ----
cat("Step 6/6: OR 森林图 + 箱线图...\n")
glm_fit <- glm(form, data = df, family = binomial("logit"))
ct <- summary(glm_fit)$coefficients; ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
or <- data.frame(Gene = rownames(ct), OR = exp(ct[, 1]),
                 lo = exp(ct[, 1] - 1.96 * ct[, 2]), hi = exp(ct[, 1] + 1.96 * ct[, 2]), P = ct[, 4])
or$Dir <- ifelse(or$OR > 1, "Risk", "Protective"); or <- or[order(or$OR), ]; or$Gene <- factor(or$Gene, levels = or$Gene)
write.csv(or, file.path(args$outdir, "OR_table.csv"), row.names = FALSE)
if (all(is.finite(c(or$OR, or$lo, or$hi))) && all(c(or$OR, or$lo, or$hi) > 0)) {
  p_or <- ggplot(or, aes(OR, Gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = .25, colour = "grey40") +
    geom_point(aes(fill = Dir), shape = 21, size = 4, colour = "black") +
    scale_fill_manual(values = c(Risk = "#BC3C29", Protective = "#0072B5"), name = NULL) +
    scale_x_log10() + labs(title = "Logistic OR (95% CI)", x = "Odds ratio (log scale)", y = NULL) +
    theme_pub(base_size = 12, border = TRUE) + theme(axis.text.y = element_text(face = "italic"))
  save_fig(p_or, file.path(ASSETS, "OR_forest"), 6, 4.5); save_fig(p_or, file.path(args$outdir, "OR_forest"), 6, 4.5)
}
bx <- do.call(rbind, lapply(feat, function(g) data.frame(Gene = g, Expr = df[[g]], Group = factor(grp, labels = c("Control", "Disease")))))
p_box <- ggplot(bx, aes(Gene, Expr, fill = Group)) +
  geom_boxplot(width = .6, outlier.size = .6, alpha = .85) +
  scale_fill_manual(values = pal_pub(2, "npg")) +
  labs(title = "Diagnostic gene expression", x = NULL, y = "Expression") +
  theme_pub(base_size = 12, border = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
save_fig(p_box, file.path(ASSETS, "Gene_boxplot"), 6, 5); save_fig(p_box, file.path(args$outdir, "Gene_boxplot"), 6, 5)
cat("完成。诊断模型图/表见", normalizePath(args$outdir), "及 assets/\n")
