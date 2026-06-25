# =============================================================================
# 编号       : 056
# 脚本名     : GEO 多队列合并 + 批次校正 (turnkey + 顶刊图)
# 分类       : 03_transcriptomics_deg
# 用途       : 自动合并目录下多个 GEO 表达矩阵(按 geneSymbol 取交集),用 limma
#              removeBatchEffect 去批次效应,输出校正前/后 PCA 与箱线图(独立图)。
# 方法/包    : limma::removeBatchEffect + prcomp;绘图共享 theme_pub.R
# 结果图     : PCA_before / PCA_after(批次着色+椭圆);Boxplot_before / Boxplot_after
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 056_GEO_merge_batch_correction.R
# 运行(自己): Rscript 056_GEO_merge_batch_correction.R --input data/cohorts_dir --outdir results/run1
# 输入规格 : --input 指向一个【目录】,内含多份 CSV,每份=一个队列的表达矩阵
#            (首列列名将被统一为 geneSymbol,其余列=样本)。批次名 = 文件名。
# 整理日期 : 2026-06-23(turnkey 重构;批次校正逻辑保持原状)
# =============================================================================

# ---- turnkey preamble ----
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(limma); library(ggplot2); library(reshape2) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input  = file.path(SCRIPT_DIR, "example_data"),
  outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 读取并合并多队列 ----
cat("Step 1/4: 合并多队列...\n")
files <- list.files(args$input, pattern = "\\.csv$", full.names = TRUE)
if (length(files) < 2) stop("--input 目录内需 ≥2 份队列 CSV。")
dl <- lapply(files, function(f) { d <- read_table_smart(f); names(d)[1] <- "geneSymbol"; d })
merged <- Reduce(function(x, y) merge(x, y, by = "geneSymbol"), dl)
rn <- merged$geneSymbol; expr <- as.matrix(merged[, -1, drop = FALSE]); rownames(expr) <- rn
storage.mode(expr) <- "double"
batch <- factor(rep(tools::file_path_sans_ext(basename(files)), vapply(dl, function(x) ncol(x) - 1, 0)))
cat("  合并后:", nrow(expr), "共有基因 x", ncol(expr), "样本;", nlevels(batch), "个批次\n")
write.csv(data.frame(geneSymbol = rn, expr, check.names = FALSE),
          file.path(args$outdir, "merged_before_correction.csv"), row.names = FALSE)

# ---- Step 2. 批次校正 ----
cat("Step 2/4: removeBatchEffect 批次校正...\n")
corrected <- removeBatchEffect(expr, batch = batch)
write.csv(data.frame(geneSymbol = rn, corrected, check.names = FALSE),
          file.path(args$outdir, "merged_after_correction.csv"), row.names = FALSE)

# ---- Step 3. PCA(前/后,各自独立图)----
cat("Step 3/4: PCA...\n")
pca_plot <- function(mat, title) {
  pc <- prcomp(t(mat), scale. = TRUE); vp <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], Batch = batch)
  ggplot(df, aes(PC1, PC2, colour = Batch, fill = Batch)) +
    stat_ellipse(geom = "polygon", level = 0.9, alpha = 0.15, linetype = 2, show.legend = FALSE) +
    geom_point(size = 3.2, alpha = 0.9) +
    scale_colour_manual(values = pal_pub(nlevels(batch), "npg")) +
    scale_fill_manual(values = pal_pub(nlevels(batch), "npg")) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", vp[1]), y = sprintf("PC2 (%.1f%%)", vp[2])) +
    theme_pub(base_size = 12, border = TRUE)
}
for (cfg in list(list(m = expr, t = "PCA before batch correction", f = "PCA_before"),
                 list(m = corrected, t = "PCA after batch correction", f = "PCA_after"))) {
  p <- pca_plot(cfg$m, cfg$t)
  save_fig(p, file.path(ASSETS, cfg$f), 6.5, 5.5); save_fig(p, file.path(args$outdir, cfg$f), 6.5, 5.5)
}

# ---- Step 4. 箱线图(前/后)----
cat("Step 4/4: 样本表达箱线图...\n")
box_plot <- function(mat, title) {
  df <- reshape2::melt(data.frame(mat, check.names = FALSE)); df$Batch <- batch[match(df$variable, colnames(mat))]
  ggplot(df, aes(variable, value, fill = Batch)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
    scale_fill_manual(values = pal_pub(nlevels(batch), "npg")) +
    labs(title = title, x = NULL, y = "Expression") +
    theme_pub(base_size = 11, border = TRUE) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
for (cfg in list(list(m = expr, t = "Before batch correction", f = "Boxplot_before"),
                 list(m = corrected, t = "After batch correction", f = "Boxplot_after"))) {
  p <- box_plot(cfg$m, cfg$t)
  save_fig(p, file.path(ASSETS, cfg$f), 9, 4.5); save_fig(p, file.path(args$outdir, cfg$f), 9, 4.5)
}
write.csv(data.frame(Sample = colnames(expr), Batch = batch), file.path(args$outdir, "sample_batch_info.csv"), row.names = FALSE)
cat("完成。结果见", normalizePath(args$outdir), ";展示图见", normalizePath(ASSETS), "\n")
