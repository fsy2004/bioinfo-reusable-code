# =============================================================================
# 编号       : R010
# 脚本名     : 差异分析.R
# 分类       : 03_GEO转录组整理与差异分析
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 进行转录组差异表达分析，输出差异基因表、火山图、热图和 PCA 图。
# 结果图     : 火山图；热图；PCA图；气泡图/点图
# 主要 R 包  : circlize; ComplexHeatmap; ggplot2; ggrepel; grid; limma; pheatmap; RColorBrewer
# 整理时间   : 2026-05-10
# =============================================================================
# ======================== 差异表达分析 ========================
# 目的：读取表达矩阵，进行差异分析并生成可视化图表
# 输入：Sample Type Matrix.csv
# 输出：差异基因列表、火山图、PCA图、热图等
# ==============================================================

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(ComplexHeatmap)
library(circlize)

set.seed(12345)

# 设置工作目录
work_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/14.差异分析"
setwd(work_dir)

cat("================== 差异表达分析开始 ==================\n\n")

# 设置阈值参数
threshold_logFC <- 0.5       # log2折叠变化阈值
threshold_adjP  <- 0.05      # 调整后P值阈值
max_display_genes <- 50      # 热图中每侧（上调和下调）展示的基因最大数量

# =============== 步骤1：读取表达矩阵 ===============
cat("步骤1: 读取表达矩阵...\n")

file_expr <- "Sample Type Matrix.csv"

# 自动检测分隔符
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")

# 读取表达文件
expr_raw <- read.table(file_expr, header = TRUE, sep = sep,
                       check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_raw) <- expr_raw[, 1]
expr_mat <- expr_raw[, -1, drop = FALSE]

# 强制转为数值型矩阵
expr_mat <- as.matrix(expr_mat)
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- rownames(expr_raw)
colnames(expr_mat) <- colnames(expr_raw)[-1]

if (!is.numeric(expr_mat)) stop("表达矩阵存在非数值列！")

cat("表达矩阵维度:", nrow(expr_mat), "基因 x", ncol(expr_mat), "样本\n\n")

# =============== 步骤2：自动判断分组信息 ===============
cat("步骤2: 自动判断分组信息...\n")

sample_names <- colnames(expr_mat)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Disease", "Unknown"))

if (any(group_info == "Unknown")) stop("分组未知的样本存在，请检查样本名后缀！")

num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Disease")

cat("分组情况：Control =", num_ctrl, ", Disease =", num_treat, "\n\n")

# =============== 步骤3：差异表达分析（LIMMA） ===============
cat("步骤3: 进行差异表达分析...\n")

group_labels <- factor(group_info, levels = c("Control", "Disease"))
design_mat <- model.matrix(~0 + group_labels)
colnames(design_mat) <- c("Control", "Disease")

fit <- lmFit(expr_mat, design_mat)
contrast_mat <- makeContrasts(Disease - Control, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

all_diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)

cat("差异分析完成！FDR<0.05基因数：", sum(all_diff_results$adj.P.Val < 0.05), "\n\n")

# 保存全部差异分析结果
write.csv(cbind(Gene = rownames(all_diff_results), all_diff_results),
          file = "DE_results.csv", row.names = FALSE)

# =============== 步骤4：筛选显著差异基因 ===============
cat("步骤4: 筛选显著差异基因...\n")

significant_DEGs <- all_diff_results[with(all_diff_results,
                                          (abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP)), ]

output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)

# 计算标准误差
SE <- ifelse(as.numeric(output_DEGs[, "t"]) != 0,
             abs(as.numeric(output_DEGs[, "logFC"]) / as.numeric(output_DEGs[, "t"])),
             NA)
output_DEGs <- cbind(output_DEGs, SE = SE)
output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors = FALSE)

# 添加上下调标记
output_DEGs$Regulation <- ifelse(as.numeric(output_DEGs$logFC) > 0, "Up", "Down")

# 调整列顺序
desired_order <- c("Gene", "logFC", "Regulation", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]

write.csv(output_DEGs, file = "DE_significant_genes.csv", row.names = FALSE)

cat("显著差异基因数量：", nrow(significant_DEGs), "\n")
cat("  上调基因：", sum(output_DEGs$Regulation == "Up"), "\n")
cat("  下调基因：", sum(output_DEGs$Regulation == "Down"), "\n\n")

# =============== 步骤5：绘制火山图 ===============
cat("步骤5: 绘制火山图...\n")

volcano_data <- all_diff_results
volcano_data$Gene <- rownames(volcano_data)

# 分组标记
volcano_data$Group <- "Not Significant"
volcano_data$Group[volcano_data$logFC > threshold_logFC & volcano_data$adj.P.Val < threshold_adjP] <- "Up regulated"
volcano_data$Group[volcano_data$logFC < -threshold_logFC & volcano_data$adj.P.Val < threshold_adjP] <- "Down regulated"

# 统计数量
up_count <- sum(volcano_data$Group == "Up regulated")
down_count <- sum(volcano_data$Group == "Down regulated")
not_count <- sum(volcano_data$Group == "Not Significant")

# 绘制火山图
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up regulated" = "#FF4500",
               "Down regulated" = "#1E90FF",
               "Not Significant" = "#808080"),
    labels = c(paste0("Up regulated (", up_count, ")"),
               paste0("Down regulated (", down_count, ")"),
               paste0("Not Significant (", not_count, ")"))
  ) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(threshold_adjP),
             linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "#2F4F4F"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

ggsave("DE_volcano.pdf", volcano_plot, width = 8, height = 8, dpi = 300)
cat("已保存：DE_volcano.pdf\n\n")

# =============== 步骤5.5：绘制渐变色火山图（标注Top20基因） ===============
cat("步骤5.5: 绘制渐变色火山图（标注Top20上调和下调基因）...\n")

# 准备数据
volcano_data2 <- all_diff_results
volcano_data2$Gene <- rownames(volcano_data2)

# 筛选显著差异基因
sig_genes <- volcano_data2[volcano_data2$adj.P.Val < threshold_adjP &
                           abs(volcano_data2$logFC) > threshold_logFC, ]

# 上调基因：按-log10(adj.P.Val)排序，取前20
up_genes <- sig_genes[sig_genes$logFC > 0, ]
up_genes <- up_genes[order(-(-log10(up_genes$adj.P.Val))), ]
up_top20 <- head(up_genes$Gene, 20)

# 下调基因：按-log10(adj.P.Val)排序，取前20
down_genes <- sig_genes[sig_genes$logFC < 0, ]
down_genes <- down_genes[order(-(-log10(down_genes$adj.P.Val))), ]
down_top20 <- head(down_genes$Gene, 20)

# 合并Top40基因
label_genes <- c(up_top20, down_top20)
volcano_data2$label_gene <- ifelse(volcano_data2$Gene %in% label_genes,
                                   volcano_data2$Gene, "")

cat("标注基因数量：上调", length(up_top20), "个，下调", length(down_top20), "个\n")

# 计算坐标轴范围
x_limit <- max(abs(volcano_data2$logFC), na.rm = TRUE) * 1.1
y_max <- max(-log10(volcano_data2$adj.P.Val), na.rm = TRUE) * 1.05

# 计算气泡大小图例的breaks
y_values <- -log10(volcano_data2$adj.P.Val)
y_values <- y_values[is.finite(y_values)]
y_range <- range(y_values, na.rm = TRUE)
size_breaks <- pretty(y_range, n = 5)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max(y_values)]

# 绘制渐变色火山图
gradient_volcano <- ggplot(volcano_data2, aes(x = logFC, y = -log10(adj.P.Val))) +
  # 绘制点，颜色按logFC渐变，大小按-log10(p值)
  geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  # 为标注基因添加黑色外框
  geom_point(data = subset(volcano_data2, label_gene != ""),
             aes(size = -log10(adj.P.Val)),
             shape = 21, fill = NA, color = "black", stroke = 1.2) +
  # 颜色渐变：深蓝-青绿-黄绿-黄-橙-红
  scale_color_gradientn(
    colors = c("#2B83BA", "#5AAE61", "#A6D96A", "#FFFFBF",
               "#FEE08B", "#FDAE61", "#D7191C"),
    values = scales::rescale(c(-6, -4, -2, 0, 2, 4, 6)),
    limits = c(-x_limit, x_limit),
    name = "log2FC"
  ) +
  # 点大小映射
  scale_size_continuous(
    range = c(0.8, 7),
    name = "-log10(p_val)",
    breaks = size_breaks,
    labels = as.character(round(size_breaks, 0))
  ) +
  # 添加阈值线
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(threshold_adjP),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  # 添加基因标签
  geom_text_repel(
    data = subset(volcano_data2, label_gene != ""),
    aes(label = label_gene),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.3,
    fontface = "italic",
    color = "black",
    show.legend = FALSE
  ) +
  # 添加"Down"和"Up"标注箭头
  annotate("segment", x = -x_limit * 0.55, xend = -x_limit * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#2B83BA", linewidth = 1.2) +
  annotate("text", x = -x_limit * 0.72, y = y_max * 0.97,
           label = paste0("Down (", down_count, ")"), color = "#2B83BA",
           fontface = "bold", size = 4.5, vjust = -0.8) +
  annotate("segment", x = x_limit * 0.55, xend = x_limit * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#D7191C", linewidth = 1.2) +
  annotate("text", x = x_limit * 0.72, y = y_max * 0.97,
           label = paste0("Up (", up_count, ")"), color = "#D7191C",
           fontface = "bold", size = 4.5, vjust = -0.8) +
  # 添加p值阈值标注
  annotate("text", x = x_limit * 0.98, y = -log10(threshold_adjP),
           label = paste0("p = ", threshold_adjP),
           hjust = 1, vjust = -0.5, size = 3.5, color = "grey30") +
  # 设置坐标轴
  scale_x_continuous(limits = c(-x_limit, x_limit), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0.02, 0)) +
  # 标题和标签
  labs(
    title = "Volcano Plot (Top 20 Genes Labeled)",
    x = "avg_log2FC",
    y = "-log10(p_val)"
  ) +
  # 主题设置
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin = margin(15, 15, 10, 10)
  ) +
  # 图例设置
  guides(
    color = guide_colorbar(
      title = "avg_log2FC",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1.2,
      barheight = 10,
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black",
      ticks.linewidth = 0.5,
      order = 1
    ),
    size = guide_legend(
      title = "-log10(p_val)",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, color = "grey30"),
      order = 2
    )
  )

ggsave("DE_volcano_gradient.pdf", gradient_volcano, width = 9, height = 7.5, dpi = 300)
cat("已保存：DE_volcano_gradient.pdf\n\n")

# =============== 步骤6：绘制PCA图 ===============
cat("步骤6: 绘制PCA图...\n")

pca_result <- prcomp(t(expr_mat), scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(expr_mat),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = factor(group_info, levels = c("Control", "Disease"))
)

pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.2, linewidth = 1) +
  geom_point(size = 3, shape = 15) +
  scale_color_manual(values = c("Control" = "#E64B35", "Disease" = "#4DBBD5")) +
  scale_fill_manual(values = c("Control" = "#E64B35", "Disease" = "#4DBBD5")) +
  labs(x = paste0("PCA1(", pca_var_perc[1], "%)"),
       y = paste0("PCA2(", pca_var_perc[2], "%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("DE_PCA.pdf", pca_plot, width = 7, height = 6, dpi = 300)
cat("已保存：DE_PCA.pdf\n\n")

# =============== 步骤7：绘制热图 ===============
cat("步骤7: 绘制差异基因热图...\n")

if (nrow(significant_DEGs) > 0) {
  ordered_DEGs <- significant_DEGs[order(as.numeric(significant_DEGs$logFC)), ]
  ordered_gene_names <- rownames(ordered_DEGs)
  total_DEG_count <- length(ordered_gene_names)

  if (total_DEG_count > (max_display_genes * 2)) {
    selected_genes <- ordered_gene_names[c(1:max_display_genes,
                                           (total_DEG_count - max_display_genes + 1):total_DEG_count)]
  } else {
    selected_genes <- ordered_gene_names
  }

  heatmap_expr <- expr_mat[selected_genes, ]
  heatmap_scaled <- t(scale(t(heatmap_expr)))

  sample_annotation <- data.frame(
    Group = factor(group_info, levels = c("Control", "Disease"))
  )
  rownames(sample_annotation) <- colnames(expr_mat)

  ha_top <- HeatmapAnnotation(
    Group = sample_annotation$Group,
    col = list(Group = c("Control" = "#66C2A5", "Disease" = "#FC8D62")),
    annotation_name_side = "left"
  )

  color_palette <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

  gene_logFC <- significant_DEGs[selected_genes, "logFC"]
  gene_labels <- paste0(selected_genes, ifelse(gene_logFC > 0, " (Up)", " (Down)"))
  names(gene_labels) <- selected_genes

  ht <- Heatmap(
    heatmap_scaled,
    name = "Z-score",
    col = color_palette,
    top_annotation = ha_top,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    row_labels = gene_labels[rownames(heatmap_scaled)],
    column_title = paste("Differential Expression Heatmap\nControl:", num_ctrl, "| Disease:", num_treat),
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )

  pdf("DE_heatmap.pdf", width = 12, height = 10)
  draw(ht)
  dev.off()

  cat("已保存：DE_heatmap.pdf\n\n")
} else {
  cat("警告：未发现显著差异基因，跳过热图绘制\n\n")
}

# =============== 步骤8：输出差异基因列表 ===============
cat("步骤8: 输出差异基因列表...\n")

diff_gene_list <- data.frame(gene = rownames(significant_DEGs))
write.table(diff_gene_list,
            file = "DEG_geneList.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

cat("已保存：DEG_geneList.txt\n\n")

# =============== 汇总报告 ===============
cat("================== 分析汇总 ==================\n")
cat("输入文件：", file_expr, "\n")
cat("样本总数：", ncol(expr_mat), "\n")
cat("  Control：", num_ctrl, "\n")
cat("  Disease：", num_treat, "\n")
cat("基因总数：", nrow(expr_mat), "\n")
cat("显著差异基因数：", nrow(significant_DEGs), "\n")
cat("  上调：", sum(output_DEGs$Regulation == "Up"), "\n")
cat("  下调：", sum(output_DEGs$Regulation == "Down"), "\n")
cat("==========================================\n\n")

cat("差异分析完成！\n")
