# =============================================================================
# 编号       : 014
# 脚本名     : 随机森林特征基因筛选 (turnkey + 顶刊图)
# 分类       : 04_ml_feature_selection
# 用途       : 用随机森林评估基因重要性(MeanDecreaseGini)筛选特征基因,输出
#              OOB 错误率曲线与重要性棒棒糖图。
# 方法/包    : randomForest;绘图共享 theme_pub.R(viridis)
# 结果图     : RF_OOB_error(错误率 vs 树数);RF_importance_lollipop(top 基因重要性)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 014_RandomForest_feature_selection.R
# 运行(自己): Rscript 014_RandomForest_feature_selection.R --input data/expr.csv --genes data/candidate.csv
# 可选参数 : --ntree 500 --top 15 --threshold 1 --seed 2025
# 输入规格 : 同 012(表达矩阵 + 可选候选基因;样本名后缀分组)。
# 整理日期 : 2026-06-23(turnkey 重构;随机森林逻辑保持原状)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(randomForest); library(ggplot2); library(reshape2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input  = file.path(SCRIPT_DIR, "example_data", "Sample_Type_Matrix.csv"),
  genes  = file.path(SCRIPT_DIR, "example_data", "candidate_genes.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  ntree = "500", top = "15", threshold = "1", seed = "2025"))
set.seed(as.integer(args$seed)); NTREE <- as.integer(args$ntree); TOP <- as.integer(args$top); THR <- as.numeric(args$threshold)
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 读取 ----
cat("Step 1/4: 读取数据...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
if (!is.null(args$genes) && file.exists(args$genes)) {
  gl <- unique(trimws(as.character(read_table_smart(args$genes)[[1]])))
  hit <- intersect(gl, rownames(expr)); if (length(hit) >= 2) expr <- expr[hit, , drop = FALSE]
}
gmat <- as.data.frame(t(as.matrix(expr)))
y <- factor(sub(".*_([^_]+)$", "\\1", rownames(gmat)))
if (nlevels(y) != 2) stop("需两组(样本名后缀)。")
orig <- colnames(gmat); colnames(gmat) <- make.names(orig); gmap <- setNames(orig, colnames(gmat))
cat("  ", nrow(gmat), "样本 x", ncol(gmat), "基因;分组:", paste(levels(y), table(y), sep = "=", collapse = " "), "\n")

# ---- Step 2. 随机森林 ----
cat("Step 2/4: 随机森林 (", NTREE, "棵树)...\n")
rf <- randomForest(x = gmat, y = y, ntree = NTREE, importance = FALSE)
imp <- importance(rf)[, "MeanDecreaseGini"]; names(imp) <- gmap[names(imp)]
imp <- sort(imp, decreasing = TRUE)
sel <- names(imp)[imp > THR]; if (length(sel) == 0) sel <- names(imp)[1]
write.csv(data.frame(gene = names(imp), MeanDecreaseGini = as.numeric(imp)),
          file.path(args$outdir, "RF_gene_importance.csv"), row.names = FALSE)
writeLines(sel, file.path(args$outdir, "RF_selected_genes.txt"))
cat("  重要性 >", THR, "的基因:", length(sel), "\n")

# ---- Step 3. OOB 错误率曲线 ----
cat("Step 3/4: 绘图...\n")
er <- as.data.frame(rf$err.rate); er$ntree <- seq_len(nrow(er))
erl <- reshape2::melt(er, id.vars = "ntree", variable.name = "Type", value.name = "Error")
erl$Type <- factor(erl$Type, levels = c("OOB", levels(y)))
p_err <- ggplot(erl, aes(ntree, Error, colour = Type)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("black", pal_pub(nlevels(y), "npg")),
                      labels = c("OOB", paste0(levels(y), " class"))) +
  labs(title = "Random forest OOB error", x = "Number of trees", y = "Error rate", colour = NULL) +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_err, file.path(ASSETS, "RF_OOB_error"), 6.5, 5); save_fig(p_err, file.path(args$outdir, "RF_OOB_error"), 6.5, 5)

# ---- Step 4. 重要性棒棒糖图 ----
df <- data.frame(Gene = names(imp), Importance = as.numeric(imp))[seq_len(min(TOP, length(imp))), ]
df$Gene <- factor(df$Gene, levels = rev(df$Gene))
p_imp <- ggplot(df, aes(Importance, Gene)) +
  geom_segment(aes(x = 0, xend = Importance, yend = Gene), colour = "grey80", linewidth = 0.8) +
  geom_point(aes(colour = Importance, size = Importance)) +
  scale_colour_viridis_c(option = "D", direction = -1, name = "Gini") +
  scale_size_continuous(range = c(3, 7), guide = "none") +
  labs(title = "Random forest gene importance", x = "Mean decrease Gini", y = NULL) +
  theme_pub(base_size = 12, border = TRUE) + theme(axis.text.y = element_text(face = "italic"))
save_fig(p_imp, file.path(ASSETS, "RF_importance_lollipop"), 6.5, 6); save_fig(p_imp, file.path(args$outdir, "RF_importance_lollipop"), 6.5, 6)
cat("完成。结果见", normalizePath(args$outdir), "\n")
