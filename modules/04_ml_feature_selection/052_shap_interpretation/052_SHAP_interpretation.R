# =============================================================================
# 编号       : 052
# 脚本名     : SHAP 机器学习解释分析 (turnkey + 顶刊图)
# 分类       : 04_ml_feature_selection
# 用途       : 训练多种 ML 分类器,选最优后用 kernelshap/shapviz 解释特征贡献,
#              输出 ROC 与 SHAP 重要性/蜂群/依赖/瀑布/力图等解释图。
# 方法/包    : caret + kernelshap + shapviz + pROC;绘图共享 theme_pub.R
# 结果图     : Model_ROC;SHAP_importance_bar;SHAP_beeswarm;SHAP_dependence;SHAP_waterfall;SHAP_force
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 052_SHAP_interpretation.R
# 运行(自己): Rscript 052_SHAP_interpretation.R --input data/geneexp.csv --outdir results/run1
# 可选参数 : --train 0.7 --cv 5 --ctrl _con --case _tra --seed 12345
# 输入规格 : 表达矩阵 CSV(首列基因,样本列名后缀分组;默认对照 *_con、实验 *_tra)。
# 整理日期 : 2026-06-23(turnkey 重构;配色由"可爱风"升级为期刊主题)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(caret); library(kernelshap); library(shapviz); library(pROC); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input = file.path(SCRIPT_DIR, "example_data", "geneexp.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  train = "0.7", cv = "5", ctrl = "_con", case = "_tra", seed = "12345"))
set.seed(as.integer(args$seed))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 数据 ----
cat("Step 1/5: 读取数据...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
dat <- as.data.frame(t(as.matrix(expr)))
grp <- ifelse(grepl(paste0(args$ctrl, "$"), rownames(dat)), "Control",
       ifelse(grepl(paste0(args$case, "$"), rownames(dat)), "Treatment", NA))
if (any(is.na(grp))) stop("样本名后缀需为 ", args$ctrl, " / ", args$case)
dat$Group <- factor(grp, levels = c("Control", "Treatment"))
idx <- createDataPartition(dat$Group, p = as.numeric(args$train), list = FALSE)
tr <- dat[idx, ]; te <- dat[-idx, ]
cat("  训练", nrow(tr), "/ 测试", nrow(te), ";特征", ncol(dat) - 1, "\n")
ctrl <- trainControl(method = "repeatedcv", number = as.integer(args$cv), classProbs = TRUE)

# ---- Step 2. 多模型 + 选最优 ----
cat("Step 2/5: 训练 RF/SVM/XGB,按 AUC 选最优...\n")
mlist <- c(RF = "rf", SVM = "svmRadial", XGB = "xgbTree")
rocs <- list(); aucs <- c(); models <- list()
for (nm in names(mlist)) {
  md <- tryCatch(train(Group ~ ., tr, method = mlist[[nm]], trControl = ctrl), error = function(e) NULL)
  if (is.null(md)) next
  pr <- predict(md, te[, -ncol(te)], type = "prob")[, "Treatment"]
  rc <- roc(ifelse(te$Group == "Control", 0, 1), as.numeric(pr), quiet = TRUE)
  rocs[[nm]] <- rc; aucs[nm] <- as.numeric(rc$auc); models[[nm]] <- md
  cat("   ", nm, "AUC=", sprintf("%.3f", aucs[nm]), "\n")
}
best <- names(which.max(aucs)); cat("  最优模型:", best, "\n")

# ROC 图(多模型)
cols <- pal_pub(length(rocs), "npg")
p_roc <- pROC::ggroc(rocs, legacy.axes = TRUE, linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  scale_colour_manual(values = cols, labels = paste0(names(rocs), " (AUC=", sprintf("%.3f", aucs[names(rocs)]), ")"), name = NULL) +
  labs(title = "Model ROC", x = "1 - Specificity", y = "Sensitivity") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(.99, .02), legend.justification = c(1, 0))
save_fig(p_roc, file.path(ASSETS, "Model_ROC"), 6, 5.5); save_fig(p_roc, file.path(args$outdir, "Model_ROC"), 6, 5.5)

# ---- Step 3. kernelshap ----
cat("Step 3/5: 计算 SHAP 值 (kernelshap)...\n")
Xbg <- tr[, -ncol(tr)]
ks <- kernelshap(models[[best]], X = Xbg, bg_X = Xbg,
                 pred_fun = function(o, newdata) predict(o, newdata, type = "prob")[, "Treatment"], verbose = FALSE)
sv <- shapviz(ks, X = Xbg)
imp <- sort(colMeans(abs(sv$S)), decreasing = TRUE); top <- names(imp)
write.csv(data.frame(Gene = names(imp), MeanAbsSHAP = as.numeric(imp)), file.path(args$outdir, "SHAP_importance.csv"), row.names = FALSE)

# ---- Step 4. SHAP 图(theme_pub)----
cat("Step 4/5: SHAP 解释图...\n")
thm <- theme_pub(base_size = 12, border = TRUE)
save_fig(sv_importance(sv, kind = "bar", fill = pal_pub(1, "npg")) + thm + labs(title = "SHAP importance"),
         file.path(ASSETS, "SHAP_importance_bar"), 6.5, 5)
save_fig(sv_importance(sv, kind = "beeswarm") + scale_colour_viridis_c(option = "D") + thm + labs(title = "SHAP beeswarm"),
         file.path(ASSETS, "SHAP_beeswarm"), 6.5, 5)
save_fig(sv_dependence(sv, v = top[1]) + thm + labs(title = paste("SHAP dependence —", top[1])),
         file.path(ASSETS, "SHAP_dependence"), 6, 5)
for (f in c("SHAP_importance_bar", "SHAP_beeswarm", "SHAP_dependence")) file.copy(file.path(ASSETS, paste0(f, ".png")), args$outdir, overwrite = TRUE)

# ---- Step 5. 单样本瀑布/力图 ----
cat("Step 5/5: 单样本 瀑布/力图...\n")
save_fig(sv_waterfall(sv, row_id = 1) + thm + labs(title = "SHAP waterfall (sample 1)"),
         file.path(ASSETS, "SHAP_waterfall"), 6.5, 5)
save_fig(sv_force(sv, row_id = 1) + thm + labs(title = "SHAP force (sample 1)"),
         file.path(ASSETS, "SHAP_force"), 7, 3.5)
cat("完成。SHAP 图/表见", normalizePath(args$outdir), "及 assets/\n")
