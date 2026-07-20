# =============================================================================
# 编号       : 034
# 脚本名     : 多种机器学习方法比较 + 特征筛选 (turnkey + 顶刊图)
# 分类       : 04_ml_feature_selection
# 用途       : 用 caret 训练多种 ML 分类器(Lasso/ElasticNet/RF/SVM/LDA/GBM/NNet/
#              PLS/kNN/LogitBoost),比较测试集 ROC/AUC,并对各法 top 特征取交集。
# 方法/包    : caret + 各算法包;重要性用 caret::varImp(免 DALEX);ROC 用 pROC。
# 结果图     : ROC_overlay(多模型 ROC);AUC_leaderboard(排行榜);Feature_UpSet(特征交集)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 034_multiML_feature_selection.R
# 运行(自己): Rscript 034_multiML_feature_selection.R --input data/expr.csv --topgene 10
# 可选参数 : --train 0.7 --cv 5 --topgene 10 --seed 123
# 输入规格 : 表达矩阵 CSV(首列基因,样本名后缀 *_con/*_tre);建议输入为小候选集
#            (如上游交集基因),避免高维下 LDA/PLS 奇异。缺失算法包的方法会自动跳过。
# 整理日期 : 2026-06-23(turnkey 重构;以 caret::varImp 替代缺失的 DALEX,逻辑等价)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(caret); library(pROC); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input  = file.path(SCRIPT_DIR, "example_data", "Sample_Type_Matrix.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  train = "0.7", cv = "5", topgene = "10", seed = "123"))
set.seed(as.integer(args$seed)); TOPG <- as.integer(args$topgene)
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 数据 ----
cat("Step 1/4: 读取数据...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
rownames(expr) <- make.names(rownames(expr))
dat <- as.data.frame(t(as.matrix(expr)))
grp <- ifelse(grepl("_con$", rownames(dat), ignore.case = TRUE), "Control", "Treatment")
dat$Type <- factor(grp, levels = c("Control", "Treatment"))
idx <- createDataPartition(dat$Type, p = as.numeric(args$train), list = FALSE)
tr <- dat[idx, ]; te <- dat[-idx, ]; yte <- ifelse(te$Type == "Control", 0, 1)
cat("  训练", nrow(tr), "/ 测试", nrow(te), "样本;特征", ncol(dat) - 1, "\n")
ctrl <- trainControl(method = "repeatedcv", number = as.integer(args$cv), classProbs = TRUE, savePredictions = TRUE)

# ---- Step 2. 多方法训练 ----
methods <- data.frame(
  Name   = c("Lasso", "ElasticNet", "RF", "SVM", "LDA", "GBM", "NeuralNet", "PLS", "kNN", "LogitBoost"),
  Method = c("glmnet", "glmnet", "rf", "svmRadial", "lda", "gbm", "nnet", "pls", "kknn", "LogitBoost"),
  stringsAsFactors = FALSE)
cat("Step 2/4: 训练", nrow(methods), "种 ML 方法...\n")
fit_one <- function(i) {
  nm <- methods$Name[i]; mt <- methods$Method[i]
  tryCatch({
    extra <- if (mt == "nnet") list(trace = FALSE) else list()      # 仅 nnet 接受 trace
    grid <- if (nm == "Lasso") expand.grid(alpha = 1, lambda = seq(.001, .1, length = 10))
            else if (nm == "ElasticNet") expand.grid(alpha = .5, lambda = seq(.001, .1, length = 10)) else NULL
    a <- c(list(Type ~ ., data = tr, method = mt, trControl = ctrl), extra)
    if (!is.null(grid)) a$tuneGrid <- grid
    md <- do.call(train, a)
    pr <- predict(md, te[, -ncol(te)], type = "prob")[, "Treatment"]
    rc <- roc(yte, as.numeric(pr), quiet = TRUE)
    vi <- tryCatch(caret::varImp(md)$importance, error = function(e) NULL)
    top <- if (!is.null(vi)) { v <- rowMeans(vi, na.rm = TRUE); names(sort(v, decreasing = TRUE))[seq_len(min(TOPG, length(v)))] } else NULL
    cat("   ", nm, " AUC=", sprintf("%.3f", as.numeric(rc$auc)), "\n")
    list(name = nm, roc = rc, auc = as.numeric(rc$auc), imp = top)
  }, error = function(e) { cat("   跳过", nm, ":", conditionMessage(e), "\n"); NULL })
}
res <- Filter(Negate(is.null), lapply(seq_len(nrow(methods)), fit_one))
if (length(res) == 0) stop("无任何模型成功,请检查数据/依赖。")
nms <- vapply(res, `[[`, "", "name")
rocs <- setNames(lapply(res, `[[`, "roc"), nms)
aucs <- setNames(vapply(res, `[[`, 0, "auc"), nms)
imps <- setNames(lapply(res, `[[`, "imp"), nms); imps <- imps[!vapply(imps, is.null, TRUE)]
write.csv(data.frame(Method = names(aucs), AUC = as.numeric(aucs)), file.path(args$outdir, "model_AUC.csv"), row.names = FALSE)

# ---- Step 3. ROC 叠加 + AUC 排行榜 ----
cat("Step 3/4: ROC + AUC 排行榜...\n")
cols <- pal_pub(length(rocs), "npg")
p_roc <- pROC::ggroc(rocs, legacy.axes = TRUE, linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  scale_colour_manual(values = cols,
    labels = paste0(names(rocs), "  (AUC=", sprintf("%.3f", aucs[names(rocs)]), ")"), name = NULL) +
  labs(title = "Model comparison — ROC", x = "1 - Specificity", y = "Sensitivity") +
  theme_pub(base_size = 12, border = TRUE) + theme(legend.position = c(0.99, 0.02),
    legend.justification = c(1, 0), legend.text = element_text(size = 8))
save_fig(p_roc, file.path(ASSETS, "ROC_overlay"), 6.5, 6); save_fig(p_roc, file.path(args$outdir, "ROC_overlay"), 6.5, 6)

ad <- data.frame(Method = names(aucs), AUC = as.numeric(aucs))
ad <- ad[order(ad$AUC), ]; ad$Method <- factor(ad$Method, levels = ad$Method)
p_auc <- ggplot(ad, aes(AUC, Method)) +                       # lollipop(顶刊优于条形)
  geom_segment(aes(x = 0.5, xend = AUC, yend = Method, colour = AUC), linewidth = 1.1) +
  geom_point(aes(colour = AUC), size = 4.5) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.25, size = 3.0, fontface = "bold") +
  scale_colour_viridis_c(option = "D", guide = "none") +
  coord_cartesian(xlim = c(0.5, 1.03)) +
  labs(title = "Model AUC leaderboard", x = "Test AUC", y = NULL) +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_auc, file.path(ASSETS, "AUC_leaderboard"), 6.5, 5); save_fig(p_auc, file.path(args$outdir, "AUC_leaderboard"), 6.5, 5)

# ---- Step 4. 特征交集 UpSet ----
cat("Step 4/4: 特征交集...\n")
if (length(imps) >= 2 && requireNamespace("UpSetR", quietly = TRUE)) {
  inter <- Reduce(intersect, imps)
  writeLines(inter, file.path(args$outdir, "intersect_genes.txt"))
  cat("  各法 top", TOPG, "特征交集:", length(inter), "个\n")
  allg <- unique(unlist(imps)); mm <- as.data.frame(sapply(imps, function(g) as.integer(allg %in% g)))
  rownames(mm) <- allg
  for (dest in c(file.path(ASSETS, "Feature_UpSet"), file.path(args$outdir, "Feature_UpSet"))) {
    grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 10, height = 6)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", point.size = 2.6, line.size = 1,
                        mainbar.y.label = "Shared features", sets.x.label = "Top features / method")); dev.off()
    grDevices::png(paste0(dest, ".png"), width = 10, height = 6, units = "in", res = 300)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", point.size = 2.6, line.size = 1,
                        mainbar.y.label = "Shared features", sets.x.label = "Top features / method")); dev.off()
  }
}
cat("完成。AUC 表/图见", normalizePath(args$outdir), "\n")
