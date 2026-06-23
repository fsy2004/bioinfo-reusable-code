# =============================================================================
# 编号       : 021
# 脚本名     : 免疫浸润可视化 (turnkey + 顶刊图)
# 分类       : 06_免疫浸润与免疫可视化
# 用途       : 对免疫细胞比例矩阵(CIBERSORT 等去卷积结果)做分组差异箱线图、
#              样本堆叠组成图、免疫细胞相关性热图。
# 方法/包    : ggplot2 + ggpubr(差异检验)+ ComplexHeatmap/corrplot;主题 theme_pub.R
# 结果图     : Immune_boxplot(分组+显著性);Immune_stackbar(组成);Immune_correlation(相关热图)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 021_immune_visualization.R
# 运行(自己): Rscript 021_immune_visualization.R --input data/CIBERSORT_results.csv
# 可选参数 : --ctrl _con --case _tre
# 输入规格 : CSV,首列 Sample(样本名,后缀分组),其余列=各免疫细胞比例(行和≈1)。
# 整理日期 : 2026-06-23(turnkey 重构;配色升级期刊风,base/pheatmap→ggplot/ComplexHeatmap)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(ggplot2); library(ggpubr); library(reshape2); library(ComplexHeatmap); library(circlize) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "CIBERSORT_results.csv"),
                      outdir = file.path(SCRIPT_DIR, "results"), ctrl = "_con", case = "_tre"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/4: 读取免疫比例矩阵...\n")
rt <- read_table_smart(args$input, row_names = TRUE)
grp <- ifelse(grepl(paste0(args$ctrl, "$"), rownames(rt)), "Control",
       ifelse(grepl(paste0(args$case, "$"), rownames(rt)), "Treat", NA))
if (any(is.na(grp))) stop("样本名后缀需为 ", args$ctrl, " / ", args$case)
cells <- colnames(rt); frac <- as.matrix(rt)
cat("  ", nrow(rt), "样本 x", length(cells), "免疫细胞\n")

# ---- Step 2. 分组差异箱线图 ----
cat("Step 2/4: 分组差异箱线图...\n")
dl <- reshape2::melt(data.frame(Sample = rownames(rt), Group = grp, rt, check.names = FALSE),
                     id.vars = c("Sample", "Group"), variable.name = "Cell", value.name = "Fraction")
p_box <- ggboxplot(dl, x = "Cell", y = "Fraction", fill = "Group", width = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = pal_pub(2, "npg")) +
  stat_compare_means(aes(group = Group), label = "p.signif", hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "*", "ns"))) +
  labs(title = "Immune cell infiltration by group", x = NULL, y = "Fraction") +
  theme_pub(base_size = 11, border = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_fig(p_box, file.path(ASSETS, "Immune_boxplot"), 8.5, 5.5); save_fig(p_box, file.path(args$outdir, "Immune_boxplot"), 8.5, 5.5)

# ---- Step 3. 堆叠组成图 ----
cat("Step 3/4: 样本堆叠组成图...\n")
dl$Sample <- factor(dl$Sample, levels = rownames(rt)[order(grp)])
p_stack <- ggplot(dl, aes(Sample, Fraction, fill = Cell)) +
  geom_col(width = 1, colour = "white", linewidth = 0.05) +
  scale_fill_manual(values = pal_pub(length(cells), "npg")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(title = "Immune cell composition", x = "Sample", y = "Proportion", fill = "Cell type") +
  theme_pub(base_size = 11, border = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.key.size = unit(.4, "cm"))
save_fig(p_stack, file.path(ASSETS, "Immune_stackbar"), 9, 5.5); save_fig(p_stack, file.path(args$outdir, "Immune_stackbar"), 9, 5.5)

# ---- Step 4. 相关性热图 ----
cat("Step 4/4: 免疫细胞相关性热图...\n")
cmat <- cor(frac, method = "spearman")
col_fun <- colorRamp2(c(-1, 0, 1), c("#3C5488", "white", "#E64B35"))
ht <- Heatmap(cmat, name = "Spearman r", col = col_fun, rect_gp = gpar(col = "white", lwd = 1),
              column_title = "Immune cell correlation", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
              cell_fun = function(j, i, x, y, w, h, fill) if (abs(cmat[i, j]) > 0.4 && i != j)
                grid.text(sprintf("%.1f", cmat[i, j]), x, y, gp = gpar(fontsize = 6)))
for (dest in c(file.path(ASSETS, "Immune_correlation"), file.path(args$outdir, "Immune_correlation"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 7, height = 6.5); draw(ht); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 7, height = 6.5, units = "in", res = 300); draw(ht); dev.off()
}
write.csv(as.data.frame(cmat), file.path(args$outdir, "immune_correlation.csv"))
cat("完成。免疫图/表见", normalizePath(args$outdir), "\n")
