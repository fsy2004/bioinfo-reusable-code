# =============================================================================
# 编号       : 013
# 脚本名     : SVM-RFE 特征基因筛选 (turnkey + 顶刊图)
# 分类       : 04_机器学习筛选特征基因
# 用途       : 用 SVM 递归特征消除(线性核权重)对基因排序,并以交叉验证准确率-
#              特征数曲线选最优特征子集。
# 方法/包    : e1071(线性 SVM)+ SVM-RFE 递归消除;绘图共享 theme_pub.R
# 结果图     : SVMRFE_CV_accuracy(准确率 vs 特征数);SVMRFE_top_rank(top 特征排名)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 013_SVM_RFE_feature_selection.R
# 运行(自己): Rscript 013_SVM_RFE_feature_selection.R --input data/expr.csv --genes data/candidate.csv
# 可选参数 : --maxk 30(评估的最大特征数) --folds 5 --seed 12345
# 输入规格 : 同 012(表达矩阵 + 可选候选基因;样本名后缀分组)。
# 整理日期 : 2026-06-23(turnkey 重构;SVM-RFE 算法保持线性核权重递归消除)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(e1071); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input  = file.path(SCRIPT_DIR, "example_data", "Sample_Type_Matrix.csv"),
  genes  = file.path(SCRIPT_DIR, "example_data", "candidate_genes.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  maxk = "30", folds = "5", seed = "12345"))
set.seed(as.integer(args$seed)); FOLDS <- as.integer(args$folds)
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 数据 ----
cat("Step 1/4: 读取数据...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
if (!is.null(args$genes) && file.exists(args$genes)) {
  gl <- unique(trimws(as.character(read_table_smart(args$genes)[[1]])))
  hit <- intersect(gl, rownames(expr)); if (length(hit) >= 2) expr <- expr[hit, , drop = FALSE]
}
X <- scale(t(as.matrix(expr)))
y <- factor(sub(".*_([^_]+)$", "\\1", rownames(X)))
if (nlevels(y) != 2) stop("需两组(样本名后缀)。")
MAXK <- min(as.integer(args$maxk), ncol(X))
cat("  ", nrow(X), "样本 x", ncol(X), "基因\n")

# ---- Step 2. SVM-RFE 递归排名 ----
cat("Step 2/4: SVM-RFE 递归消除排名...\n")
svm_weights <- function(Xs, ys) {
  m <- svm(Xs, ys, cost = 10, scale = FALSE, type = "C-classification", kernel = "linear")
  w <- t(m$coefs) %*% m$SV; as.numeric(w)^2
}
surv <- seq_len(ncol(X)); ranked <- integer(0)
while (length(surv) > 0) {
  if (length(surv) == 1) { ranked <- c(surv, ranked); break }
  sc <- svm_weights(X[, surv, drop = FALSE], y)
  drop <- which.min(sc); ranked <- c(surv[drop], ranked); surv <- surv[-drop]
}
ranked_genes <- colnames(X)[ranked]   # 重要性从高到低
write.csv(data.frame(rank = seq_along(ranked_genes), gene = ranked_genes),
          file.path(args$outdir, "SVMRFE_ranking.csv"), row.names = FALSE)

# ---- Step 3. CV 准确率 vs 特征数 ----
cat("Step 3/4: 交叉验证选最优特征数...\n")
folds <- sample(rep(seq_len(FOLDS), length.out = nrow(X)))
cv_acc <- sapply(seq_len(MAXK), function(k) {
  feats <- ranked_genes[seq_len(k)]
  acc <- sapply(seq_len(FOLDS), function(f) {
    tr <- folds != f; te <- !tr
    m <- svm(X[tr, feats, drop = FALSE], y[tr], kernel = "linear", cost = 10, scale = FALSE, type = "C-classification")
    mean(predict(m, X[te, feats, drop = FALSE]) == y[te])
  })
  mean(acc)
})
best_k <- which.max(cv_acc); best_genes <- ranked_genes[seq_len(best_k)]
writeLines(best_genes, file.path(args$outdir, "SVMRFE_selected_genes.txt"))
cat("  最优特征数:", best_k, "(CV 准确率", sprintf("%.3f", max(cv_acc)), ")\n")

dcv <- data.frame(k = seq_len(MAXK), acc = cv_acc)
p_cv <- ggplot(dcv, aes(k, acc)) +
  geom_line(colour = pal_pub(1, "npg"), linewidth = 0.8) + geom_point(colour = pal_pub(1, "npg"), size = 1.8) +
  geom_vline(xintercept = best_k, linetype = "dashed", colour = "#3C5488") +
  annotate("point", x = best_k, y = max(cv_acc), size = 3.5, shape = 21, fill = "#E64B35", colour = "black") +
  annotate("text", x = best_k, y = max(cv_acc), label = paste0("  n=", best_k), hjust = 0, fontface = "bold", size = 3.5) +
  labs(title = "SVM-RFE feature selection", x = "Number of features", y = paste0(FOLDS, "-fold CV accuracy")) +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_cv, file.path(ASSETS, "SVMRFE_CV_accuracy"), 6.5, 5); save_fig(p_cv, file.path(args$outdir, "SVMRFE_CV_accuracy"), 6.5, 5)

# ---- Step 4. top 特征排名图 ----
cat("Step 4/4: top 特征排名...\n")
nshow <- min(best_k + 5, 20, length(ranked_genes))
dr <- data.frame(gene = ranked_genes[seq_len(nshow)], rank = seq_len(nshow))
dr$gene <- factor(dr$gene, levels = rev(dr$gene)); dr$selected <- dr$rank <= best_k
p_rank <- ggplot(dr, aes(rank, gene)) +
  geom_segment(aes(x = 0, xend = rank, yend = gene, colour = selected), linewidth = 0.8) +
  geom_point(aes(colour = selected), size = 3) +
  scale_colour_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "grey65"),
                      labels = c("not selected", "selected"), name = NULL) +
  scale_x_reverse() +
  labs(title = "SVM-RFE feature ranking", x = "RFE rank (1 = most important)", y = NULL) +
  theme_pub(base_size = 12, border = TRUE) + theme(axis.text.y = element_text(face = "italic"))
save_fig(p_rank, file.path(ASSETS, "SVMRFE_top_rank"), 6.5, 6); save_fig(p_rank, file.path(args$outdir, "SVMRFE_top_rank"), 6.5, 6)
cat("完成。排名/最优子集/图见", normalizePath(args$outdir), "\n")
