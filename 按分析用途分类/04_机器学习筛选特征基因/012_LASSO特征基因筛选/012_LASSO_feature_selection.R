# =============================================================================
# 编号       : 012
# 脚本名     : LASSO 回归特征基因筛选 (turnkey + 顶刊图)
# 分类       : 04_机器学习筛选特征基因
# 用途       : 用 LASSO(L1 正则 logistic)从候选基因中筛选与分组相关的特征基因,
#              输出交叉验证曲线与系数收缩路径图,以及入选基因列表/表达矩阵。
# 方法/包    : glmnet(family=binomial, alpha=1)+ cv.glmnet;绘图共享 theme_pub.R
# 结果图     : LASSO_CV_curve(偏差 vs logλ);LASSO_coefficient_path(系数路径)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 012_LASSO_feature_selection.R
# 运行(自己): Rscript 012_LASSO_feature_selection.R --input data/expr.csv --genes data/candidate.csv
# 可选参数 : --lambda min|1se(默认 min) --nfolds 10 --seed 123456
# 输入规格 : --input 表达矩阵 CSV(首列基因,样本列名后缀分组,如 *_con/*_tre);
#            --genes 候选基因列表 CSV(首列基因名,可选;省略则用全部基因)。
# 整理日期 : 2026-06-23(turnkey 重构;LASSO 逻辑保持原状)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(glmnet); library(ggplot2); library(ggrepel); library(reshape2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input  = file.path(SCRIPT_DIR, "example_data", "Sample_Type_Matrix.csv"),
  genes  = file.path(SCRIPT_DIR, "example_data", "candidate_genes.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  lambda = "min", nfolds = "10", seed = "123456"))
set.seed(as.integer(args$seed))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 读表达 + 可选候选基因 ----
cat("Step 1/4: 读取数据...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
if (!is.null(args$genes) && file.exists(args$genes)) {
  gl <- read_table_smart(args$genes); gl <- unique(trimws(as.character(gl[[1]])))
  hit <- intersect(gl, rownames(expr))
  if (length(hit) >= 2) { expr <- expr[hit, , drop = FALSE]; cat("  限定候选基因:", length(hit), "/", length(gl), "\n") }
}
x <- t(as.matrix(expr)); storage.mode(x) <- "double"
y <- sub(".*_([^_]+)$", "\\1", rownames(x))
if (length(unique(y)) != 2) stop("需恰好两组(样本名后缀),当前: ", paste(unique(y), collapse = ", "))
cat("  ", nrow(x), "样本 x", ncol(x), "基因;分组:", paste(names(table(y)), table(y), sep = "=", collapse = " "), "\n")

# ---- Step 2. LASSO + 交叉验证 ----
cat("Step 2/4: LASSO + ", args$nfolds, "折交叉验证...\n")
fit <- glmnet(x, y, family = "binomial", alpha = 1)
cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = as.integer(args$nfolds))
lam <- if (args$lambda == "1se") cv$lambda.1se else cv$lambda.min
co <- as.matrix(coef(fit, s = lam))[-1, 1]
sel <- names(co)[co != 0]; sel <- sel[order(abs(co[sel]), decreasing = TRUE)]
write.csv(data.frame(gene = sel, coef = co[sel]), file.path(args$outdir, "LASSO_selected_genes.csv"), row.names = FALSE)
cat("  入选特征基因:", length(sel), "(λ.", args$lambda, ")\n")

# ---- Step 3. CV 曲线 ----
cat("Step 3/4: 绘图...\n")
cvd <- data.frame(loglam = log(cv$lambda), cvm = cv$cvm, lo = cv$cvlo, up = cv$cvup,
                  n = cv$nzero)
p_cv <- ggplot(cvd, aes(loglam, cvm)) +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.04, colour = "grey70") +
  geom_point(colour = pal_pub(1, "npg"), size = 2) +
  geom_vline(xintercept = log(cv$lambda.min), linetype = "dashed", colour = "#E64B35") +
  geom_vline(xintercept = log(cv$lambda.1se), linetype = "dashed", colour = "#3C5488") +
  annotate("text", x = log(cv$lambda.min), y = max(cvd$up), label = "λ.min", colour = "#E64B35", hjust = -0.1, size = 3.3) +
  annotate("text", x = log(cv$lambda.1se), y = max(cvd$up), label = "λ.1se", colour = "#3C5488", hjust = -0.1, size = 3.3) +
  labs(title = "LASSO cross-validation", x = expression(log(lambda)), y = "Binomial deviance") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_cv, file.path(ASSETS, "LASSO_CV_curve"), 6.5, 5.5); save_fig(p_cv, file.path(args$outdir, "LASSO_CV_curve"), 6.5, 5.5)

# ---- Step 4. 系数收缩路径 ----
beta <- as.matrix(fit$beta); ll <- log(fit$lambda)
long <- reshape2::melt(data.frame(loglam = ll, t(beta)), id.vars = "loglam", variable.name = "gene", value.name = "coef")
rep_pts <- do.call(rbind, lapply(split(long, long$gene), function(d) d[which.max(abs(d$coef)), ]))
rep_lab <- rep_pts[rep_pts$gene %in% sel, ]   # 仅标注入选基因
p_path <- ggplot(long, aes(loglam, coef, group = gene, colour = gene)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  geom_vline(xintercept = log(lam), linetype = "dashed", colour = "grey40") +
  ggrepel::geom_text_repel(data = rep_lab, aes(label = gene), size = 2.8, fontface = "italic",
                           max.overlaps = 30, segment.size = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = colorRampPalette(pal_pub(name = "npg"))(nrow(beta)), guide = "none") +
  labs(title = "LASSO coefficient path", x = expression(log(lambda)), y = "Coefficient") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_path, file.path(ASSETS, "LASSO_coefficient_path"), 7, 5.8); save_fig(p_path, file.path(args$outdir, "LASSO_coefficient_path"), 7, 5.8)

write.csv(expr[intersect(sel, rownames(expr)), , drop = FALSE], file.path(args$outdir, "LASSO_gene_expression.csv"))
cat("完成。入选基因 + 图见", normalizePath(args$outdir), "\n")
