# =============================================================================
# 编号       : 022
# 脚本名     : 分子对接结合能可视化 (turnkey + 顶刊图)
# 分类       : 07_分子对接与结合能可视化
# 用途       : 把化合物×靶点的对接结合能矩阵绘成热图(蓝=强结合)与排序气泡图。
# 方法/包    : ComplexHeatmap + ggplot2;主题 theme_pub.R
# 结果图     : Binding_heatmap(化合物×靶点);Binding_bubble(最强结合排序)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 022_docking_binding_energy.R
# 运行(自己): Rscript 022_docking_binding_energy.R --input data/binding_energy.csv
# 输入规格 : CSV,首列 Target(靶点名),其余列=各化合物的对接结合能(kcal/mol,越负越强)。
#            单化合物时也可只有一列。
# 整理日期 : 2026-06-23(turnkey 重构;扩展为矩阵热图,配色升级期刊风)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(ggplot2); library(ComplexHeatmap); library(circlize) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "binding_energy.csv"),
                      outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/3: 读取结合能矩阵...\n")
be <- read_table_smart(args$input, row_names = TRUE)
M <- as.matrix(be); storage.mode(M) <- "double"
cat("  ", nrow(M), "靶点 x", ncol(M), "化合物\n")

# ---- Step 2. 结合能热图(蓝=强结合)----
cat("Step 2/3: 结合能热图...\n")
rng <- range(M, na.rm = TRUE)
col_fun <- colorRamp2(c(rng[1], mean(rng), max(rng[2], -3)), c("#08519C", "#6BAED6", "#F7F7F7"))
ht <- Heatmap(M, name = "kcal/mol", col = col_fun, rect_gp = gpar(col = "white", lwd = 1.2),
              cluster_rows = TRUE, cluster_columns = ncol(M) > 2,
              column_title = "Docking binding energy", column_title_gp = gpar(fontsize = 13, fontface = "bold"),
              row_names_gp = gpar(fontsize = 10, fontface = "italic"), column_names_gp = gpar(fontsize = 10),
              cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.1f", M[i, j]), x, y, gp = gpar(fontsize = 8,
                col = ifelse(M[i, j] < mean(rng), "white", "black"))))
for (dest in c(file.path(ASSETS, "Binding_heatmap"), file.path(args$outdir, "Binding_heatmap"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = max(5, ncol(M) * 1.1 + 2), height = max(4, nrow(M) * 0.5 + 1.5)); draw(ht); dev.off()
  grDevices::png(paste0(dest, ".png"), width = max(5, ncol(M) * 1.1 + 2), height = max(4, nrow(M) * 0.5 + 1.5), units = "in", res = 300); draw(ht); dev.off()
}

# ---- Step 3. 最强结合排序气泡图 ----
cat("Step 3/3: 最强结合气泡图...\n")
best <- data.frame(Target = rownames(M), Compound = colnames(M)[apply(M, 1, which.min)],
                   BindingEnergy = apply(M, 1, min, na.rm = TRUE))
best <- best[order(best$BindingEnergy), ]; best$Target <- factor(best$Target, levels = rev(best$Target))
write.csv(best, file.path(args$outdir, "best_binding.csv"), row.names = FALSE)
p <- ggplot(best, aes(BindingEnergy, Target)) +
  geom_segment(aes(x = 0, xend = BindingEnergy, yend = Target), colour = "grey80", linewidth = 0.6) +
  geom_point(aes(fill = BindingEnergy), shape = 21, size = 8, colour = "black", stroke = 0.6) +
  geom_text(aes(label = sprintf("%.1f", BindingEnergy)), size = 2.8, fontface = "bold", colour = "white") +
  scale_fill_gradient(low = "#08519C", high = "#9ECAE1", name = "kcal/mol") +
  labs(title = "Strongest binding per target", x = "Binding energy (kcal/mol)", y = NULL) +
  theme_pub(base_size = 12, border = TRUE) + theme(axis.text.y = element_text(face = "italic"))
save_fig(p, file.path(ASSETS, "Binding_bubble"), 6.5, 5); save_fig(p, file.path(args$outdir, "Binding_bubble"), 6.5, 5)
cat("完成。对接图/表见", normalizePath(args$outdir), "\n")
