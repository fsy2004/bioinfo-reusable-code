# =============================================================================
# 编号       : 010
# 脚本名     : GEO 差异表达分析 — 火山图 / 热图 / PCA (turnkey + 顶刊图)
# 分类       : 03_GEO转录组整理与差异分析
# 用途       : 对表达矩阵做 limma 差异表达分析,输出 DEG 表与顶刊级独立图:
#              渐变火山图(标注 top 基因)、PCA 散点(95% 椭圆)、DEG 聚类热图。
# 方法/包    : limma(lmFit/eBayes) + prcomp + ComplexHeatmap;绘图共享 theme_pub.R
# 结果图     : DEG_volcano(渐变,top标注);DEG_PCA;DEG_heatmap
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 010_GEO_DEG_volcano_heatmap_PCA.R
# 运行(自己): Rscript 010_GEO_DEG_volcano_heatmap_PCA.R --input data/expr.csv --outdir results/run1
# 可选参数 : --logfc 0.5  --padj 0.05  --topn 20  --ctrl _con  --case _tre
# 输入规格 : CSV 表达矩阵,首列=基因名,其余每列=一个样本;样本名后缀区分分组
#            (默认对照 *_con、实验 *_tre)。表达值建议已 log2 归一化。
# 整理日期 : 2026-06-23(turnkey 重构;limma 分析逻辑保持原状)
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
suppressWarnings(suppressMessages({
  library(limma); library(ggplot2); library(ggrepel); library(ComplexHeatmap); library(circlize)
}))
set.seed(12345)

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input  = file.path(SCRIPT_DIR, "example_data", "expr_matrix.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  logfc = "0.5", padj = "0.05", topn = "20", ctrl = "_con", case = "_tre"))
LOGFC <- as.numeric(args$logfc); PADJ <- as.numeric(args$padj); TOPN <- as.integer(args$topn)
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 读表达矩阵 ----
cat("Step 1/6: 读取表达矩阵...\n")
expr_raw <- read_table_smart(args$input, row_names = TRUE)
expr_mat <- as.matrix(expr_raw); storage.mode(expr_mat) <- "double"
cat("  维度:", nrow(expr_mat), "基因 x", ncol(expr_mat), "样本\n")

# ---- Step 2. 自动分组(按样本名后缀)----
cat("Step 2/6: 判断分组...\n")
sn <- colnames(expr_mat)
grp <- ifelse(grepl(paste0(args$ctrl, "$"), sn, ignore.case = TRUE), "Control",
       ifelse(grepl(paste0(args$case, "$"), sn, ignore.case = TRUE), "Disease", NA))
if (any(is.na(grp))) stop(sprintf("有样本名不含分组后缀(%s / %s),请检查或用 --ctrl/--case 指定。", args$ctrl, args$case))
nC <- sum(grp == "Control"); nD <- sum(grp == "Disease")
cat("  Control =", nC, " Disease =", nD, "\n")

# ---- Step 3. limma 差异分析 ----
cat("Step 3/6: limma 差异表达分析...\n")
gl <- factor(grp, levels = c("Control", "Disease"))
design <- model.matrix(~0 + gl); colnames(design) <- c("Control", "Disease")
fit <- eBayes(contrasts.fit(lmFit(expr_mat, design), makeContrasts(Disease - Control, levels = design)))
res <- topTable(fit, adjust.method = "fdr", number = Inf)
res$Gene <- rownames(res)
write.csv(res[, c("Gene", setdiff(names(res), "Gene"))], file.path(args$outdir, "DE_results.csv"), row.names = FALSE)
sig <- subset(res, abs(logFC) > LOGFC & adj.P.Val < PADJ)
sig$Regulation <- ifelse(sig$logFC > 0, "Up", "Down")
write.csv(sig[, c("Gene", "logFC", "Regulation", "AveExpr", "t", "P.Value", "adj.P.Val")],
          file.path(args$outdir, "DE_significant_genes.csv"), row.names = FALSE)
nUp <- sum(sig$Regulation == "Up"); nDn <- sum(sig$Regulation == "Down")
cat("  显著 DEG:", nrow(sig), " (上调", nUp, "/ 下调", nDn, ")\n")

# ---- Step 4. 渐变火山图(标注 top 基因)----
cat("Step 4/6: 火山图...\n")
v <- res; v$Group <- "Not sig"
v$Group[v$logFC >  LOGFC & v$adj.P.Val < PADJ] <- "Up"
v$Group[v$logFC < -LOGFC & v$adj.P.Val < PADJ] <- "Down"
xl <- max(abs(v$logFC), na.rm = TRUE) * 1.1; ym <- max(-log10(v$adj.P.Val), na.rm = TRUE) * 1.05
lab <- rbind(head(v[v$Group == "Up", ][order(v[v$Group == "Up", ]$adj.P.Val), ], TOPN),
             head(v[v$Group == "Down", ][order(v[v$Group == "Down", ]$adj.P.Val), ], TOPN))
p_volcano <- ggplot(v, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(colour = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  geom_point(data = lab, shape = 21, fill = NA, colour = "black", stroke = 0.9,
             aes(size = -log10(adj.P.Val))) +
  scale_colour_gradientn(colours = c("#2B83BA", "#5AAE61", "#A6D96A", "#FFFFBF",
                                     "#FEE08B", "#FDAE61", "#D7191C"),
                         limits = c(-xl, xl), name = expression(log[2]~FC)) +
  scale_size_continuous(range = c(0.8, 6), name = expression(-log[10]~P[adj])) +
  geom_vline(xintercept = c(-LOGFC, LOGFC), linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  ggrepel::geom_text_repel(data = lab, aes(label = Gene), size = 2.8, fontface = "italic",
                           max.overlaps = Inf, box.padding = 0.4, segment.colour = "grey60",
                           segment.size = 0.25, colour = "black") +
  annotate("text", x = -xl * 0.7, y = ym, label = paste0("Down (", nDn, ")"),
           colour = "#2B83BA", fontface = "bold", size = 4) +
  annotate("text", x = xl * 0.7, y = ym, label = paste0("Up (", nUp, ")"),
           colour = "#D7191C", fontface = "bold", size = 4) +
  scale_x_continuous(limits = c(-xl, xl)) +
  labs(title = "Differential expression volcano", x = expression(log[2]~fold~change),
       y = expression(-log[10]~italic(P)[adj])) +
  theme_pub(base_size = 12, border = TRUE) +
  guides(colour = guide_colorbar(barwidth = 0.9, barheight = 7, frame.colour = "black", order = 1))
save_fig(p_volcano, file.path(ASSETS, "DEG_volcano"), width = 8.5, height = 7)
save_fig(p_volcano, file.path(args$outdir, "DEG_volcano"), width = 8.5, height = 7)

# ---- Step 5. PCA ----
cat("Step 5/6: PCA...\n")
pca <- prcomp(t(expr_mat), scale. = TRUE)
vp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pdf_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = gl)
p_pca <- ggplot(pdf_df, aes(PC1, PC2, colour = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.18, linewidth = 0.7) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = pal_pub(2, "npg")) +
  scale_fill_manual(values = pal_pub(2, "npg")) +
  labs(title = "PCA of samples", x = sprintf("PC1 (%.1f%%)", vp[1]), y = sprintf("PC2 (%.1f%%)", vp[2])) +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_pca, file.path(ASSETS, "DEG_PCA"), width = 6.5, height = 5.5)
save_fig(p_pca, file.path(args$outdir, "DEG_PCA"), width = 6.5, height = 5.5)

# ---- Step 6. DEG 热图(ComplexHeatmap)----
cat("Step 6/6: 热图...\n")
if (nrow(sig) > 0) {
  ord <- sig[order(sig$logFC), "Gene"]; n2 <- min(50, floor(length(ord) / 2))
  sel <- if (length(ord) > 2 * n2) c(head(ord, n2), tail(ord, n2)) else ord
  hm <- t(scale(t(expr_mat[sel, , drop = FALSE])))
  ha <- HeatmapAnnotation(Group = gl, col = list(Group = setNames(pal_pub(2, "npg"), c("Control", "Disease"))),
                          annotation_name_side = "left", annotation_legend_param = list(Group = list(title = "Group")))
  col_fun <- colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
  ht <- Heatmap(hm, name = "Z-score", col = col_fun, top_annotation = ha,
                cluster_columns = TRUE, show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                column_title = sprintf("Top DEGs (Control n=%d | Disease n=%d)", nC, nD),
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                heatmap_legend_param = list(title = "Z-score"))
  for (dest in c(file.path(ASSETS, "DEG_heatmap"), file.path(args$outdir, "DEG_heatmap"))) {
    grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 9, height = 9); draw(ht); dev.off()
    grDevices::png(paste0(dest, ".png"), width = 9, height = 9, units = "in", res = 300); draw(ht); dev.off()
  }
  cat("  热图基因数:", length(sel), "\n")
} else cat("  无显著 DEG,跳过热图。\n")

cat("完成。结果见", normalizePath(args$outdir), ";展示图见", normalizePath(ASSETS), "\n")
