# =============================================================================
# 编号       : 063
# 脚本名     : GEO 诊断模型外部验证 (turnkey + 顶刊图)
# 分类       : 05_诊断模型与验证
# 用途       : 在训练队列拟合 logistic 诊断模型,在独立验证队列做外部验证,输出
#              训练/验证 ROC 对比、验证集校准曲线、各队列 AUC。
# 方法/包    : rms::lrm + pROC;绘图共享 theme_pub.R
# 结果图     : ROC_train_vs_valid;Calibration_valid
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 063_diagnostic_validation.R
# 运行(自己): Rscript 063_diagnostic_validation.R --train data/train.csv --valid data/valid.csv --genes data/genes.csv
# 输入规格 : --train / --valid 两份表达矩阵(首列基因,样本名后缀分组 *_con/*_dis,基因需可对齐);
#            --genes 诊断基因列表。验证集省略时退化为训练集自评。
# 整理日期 : 2026-06-23(turnkey 重构;聚焦外部验证,与 016 内部评价互补)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(rms); library(pROC); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  train = file.path(SCRIPT_DIR, "example_data", "train_matrix.csv"),
  valid = file.path(SCRIPT_DIR, "example_data", "validation_matrix.csv"),
  genes = file.path(SCRIPT_DIR, "example_data", "diagnostic_genes.csv"),
  outdir = file.path(SCRIPT_DIR, "results"), ctrl = "_con", case = "_dis"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

prep <- function(path, feat) {
  e <- read_table_smart(path, row_names = TRUE); f <- intersect(feat, rownames(e))
  d <- as.data.frame(t(as.matrix(e[f, , drop = FALSE])))
  d$Group <- ifelse(grepl(paste0(args$ctrl, "$"), rownames(d)), 0L,
             ifelse(grepl(paste0(args$case, "$"), rownames(d)), 1L, NA))
  if (any(is.na(d$Group))) stop("样本名后缀需为 ", args$ctrl, " / ", args$case, " (", basename(path), ")")
  d
}

cat("Step 1/3: 拟合训练模型...\n")
gl <- unique(trimws(as.character(read_table_smart(args$genes)[[1]])))
tr <- prep(args$train, gl); feat <- setdiff(colnames(tr), "Group")
has_valid <- !is.null(args$valid) && file.exists(args$valid)
va <- if (has_valid) prep(args$valid, feat) else tr
feat <- intersect(feat, setdiff(colnames(va), "Group")); if (length(feat) < 2) stop("可对齐基因过少。")
dd <- datadist(tr); options(datadist = "dd")
fit <- lrm(as.formula(paste("Group ~", paste(feat, collapse = " + "))), data = tr, x = TRUE, y = TRUE)

cat("Step 2/3: 训练 / 验证 ROC...\n")
p_tr <- predict(fit, tr, type = "fitted"); p_va <- predict(fit, va, type = "fitted")
roc_tr <- roc(tr$Group, p_tr, quiet = TRUE); roc_va <- roc(va$Group, p_va, quiet = TRUE)
auc_tr <- as.numeric(roc_tr$auc); auc_va <- as.numeric(roc_va$auc)
cat("  训练 AUC=", sprintf("%.3f", auc_tr), if (has_valid) paste0(" · 验证 AUC=", sprintf("%.3f", auc_va)) else "", "\n")
write.csv(data.frame(cohort = c("Train", "Validation"), AUC = c(auc_tr, auc_va)), file.path(args$outdir, "AUC.csv"), row.names = FALSE)
rl <- list(); rl[[sprintf("Train (AUC=%.3f)", auc_tr)]] <- roc_tr
if (has_valid) rl[[sprintf("Validation (AUC=%.3f)", auc_va)]] <- roc_va
p_roc <- pROC::ggroc(rl, legacy.axes = TRUE, linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  scale_colour_manual(values = pal_pub(length(rl), "npg"), name = NULL) +
  labs(title = "Diagnostic model — external validation", x = "1 - Specificity", y = "Sensitivity") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.99, .02), legend.justification = c(1, 0))
save_fig(p_roc, file.path(ASSETS, "ROC_train_vs_valid"), 6, 5.5); save_fig(p_roc, file.path(args$outdir, "ROC_train_vs_valid"), 6, 5.5)

cat("Step 3/3: 验证集校准曲线...\n")
cohort <- if (has_valid) va else tr; pc <- if (has_valid) p_va else p_tr
br <- cut(pc, breaks = quantile(pc, seq(0, 1, .25), na.rm = TRUE), include.lowest = TRUE)
cdf <- aggregate(data.frame(pred = pc, obs = cohort$Group), list(bin = br), mean)
p_cal <- ggplot(cdf, aes(pred, obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#E64B35", linewidth = 0.8) + geom_point(colour = "#E64B35", size = 3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = paste0("Calibration (", if (has_valid) "validation" else "train", ")"),
       x = "Predicted probability", y = "Observed frequency") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_cal, file.path(ASSETS, "Calibration_valid"), 5.5, 5.5); save_fig(p_cal, file.path(args$outdir, "Calibration_valid"), 5.5, 5.5)
cat("完成。验证图/表见", normalizePath(args$outdir), "\n")
