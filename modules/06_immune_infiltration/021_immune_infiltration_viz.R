# =============================================================================
# 编号       : R021
# 脚本名     : 免疫浸润分析可视化.R
# 分类       : 06_immune_infiltration
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 对免疫浸润评分结果进行差异、相关性和分组可视化。
# 结果图     : 热图；条形图/柱状图；气泡图/点图；箱线图；相关性图
# 主要 R 包  : broom; dplyr; ggplot2; ggpubr; ggsci; limma; linkET; pheatmap; RColorBrewer; reshape2; scales; tidyverse
# 整理时间   : 2026-05-10
# =============================================================================
# ==== 免疫浸润分析完整流程（可视化 + 基因-免疫相关性） ====

# ============================ 加载必要包 =============================
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)
library(ggsci)
library(limma)
library(pheatmap)
library(scales)
library(tidyverse)
library(linkET)

# ============================ 参数设置 =============================
work_dir       <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/23.免疫浸润可视化"
input_file     <- "Sample Type Matrix.csv"
cibersort_file <- "CIBERSORT-Results.csv"
gene_list_file <- "IntersectionGenes.csv"

setwd(work_dir)

# ============================================================
# 第一部分：免疫浸润可视化（箱线图 + 堆叠条形图）
# ============================================================
cat("==== 第一部分：免疫浸润可视化 ====\n")

# 读取CIBERSORT结果
rt <- read.table(cibersort_file, header=TRUE, sep=",", check.names=FALSE, row.names=1)

# 根据样本后缀判断分组
rt$Group <- ifelse(grepl("_con$", rownames(rt), ignore.case=TRUE), "Control",
                   ifelse(grepl("_tre$", rownames(rt), ignore.case=TRUE), "Treat", NA))
rt$Sample <- rownames(rt)

control_samples <- rownames(rt)[rt$Group == "Control"]
treat_samples   <- rownames(rt)[rt$Group == "Treat"]
all_samples_ordered <- c(control_samples, "gap", treat_samples)

# 转为长格式
data_long <- melt(rt, id.vars=c("Sample","Group"), variable.name="Immune", value.name="Fraction")
data_long$Sample <- factor(data_long$Sample, levels=all_samples_ordered)

countControl <- length(unique(data_long$Sample[data_long$Group == "Control"]))
countTreat   <- length(unique(data_long$Sample[data_long$Group == "Treat"]))

# ---- 箱线图：两组免疫细胞比较 ----
boxplot_immune <- ggboxplot(
  data_long, x="Immune", y="Fraction", fill="Group",
  palette=c("Control"="#FFC0CB", "Treat"="#87CEFA"),
  xlab="", ylab="Fraction", legend.title="Group", notch=FALSE, width=0.8
) +
  stat_compare_means(
    aes(group=Group), label="p.signif",
    symnum.args=list(cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","ns")),
    label.y.npc="top", vjust=-0.5
  ) +
  theme_classic(base_size=14) +
  theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1),
        axis.line=element_line(color="black", size=1),
        axis.ticks=element_line(color="black", size=1),
        plot.margin=margin(20,10,10,10)) +
  labs(title="Immune Cell Comparison",
       subtitle=paste0(" (Control n=", countControl, ", Treat n=", countTreat, ")")) +
  coord_cartesian(clip="off")

ggsave("immune_diff-points_n.pdf", boxplot_immune, width=8, height=6)
cat(">> 箱线图已保存为：immune_diff-points_n.pdf\n")

# 导出统计表
summary_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune, Group) %>%
  summarise(MeanFraction=mean(Fraction, na.rm=TRUE), MedianFraction=median(Fraction, na.rm=TRUE),
            SD=sd(Fraction, na.rm=TRUE), Count=n(), .groups="drop")
write.csv(summary_table, file="boxplot_summary_table.csv", row.names=FALSE)

pvalue_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune) %>%
  do(tidy(wilcox.test(Fraction ~ Group, data=.))) %>%
  dplyr::select(Immune, p.value)
write.csv(pvalue_table, file="immune_pvalues.csv", row.names=FALSE)
cat(">> 统计汇总表和p值表已保存\n")

# ---- 堆叠条形图 ----
immune_types <- unique(data_long$Immune)
nColors <- length(immune_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

barplot_immune <- ggplot(data_long, aes(x=Sample, y=Fraction, fill=Immune)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=myColors) +
  scale_x_discrete(drop=FALSE) +
  theme_minimal(base_size=18) +
  theme(text=element_text(face="bold"), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x=element_blank()) +
  labs(x=NULL, y="Relative Percent", fill="Immune\nCell Type",
       title="Immune Cell Distribution", subtitle="Control vs. Treat") +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.05)))

control_count <- length(control_samples)
treat_count   <- length(treat_samples)

barplot_immune <- barplot_immune +
  annotate("segment", x=0.5, xend=control_count+0.5, y=-0.04, yend=-0.04, color="#D65DB1", size=5) +
  annotate("text", x=(control_count)/2+0.5, y=-0.08, label="Control", color="#D65DB1", size=7, fontface="bold") +
  annotate("segment", x=control_count+1.5, xend=control_count+treat_count+1.5, y=-0.04, yend=-0.04, color="#0089BA", size=5) +
  annotate("text", x=control_count+(treat_count)/2+1.5, y=-0.08, label="Treat", color="#0089BA", size=7, fontface="bold")

ggsave("barplot_two_lines_set3_gap.pdf", barplot_immune, width=12, height=7.5)
cat(">> 堆叠条形图已保存为：barplot_two_lines_set3_gap.pdf\n")

# ============================================================
# 第二部分：关键基因与免疫浸润细胞的相关性分析
# ============================================================
cat("\n==== 第二部分：基因-免疫浸润细胞相关性分析 ====\n")

# 读取基因表达数据
exprRaw <- read.table(input_file, header=TRUE, sep=",", check.names=FALSE)
exprMat <- as.matrix(exprRaw)
rownames(exprMat) <- exprMat[, 1]
exprVals <- exprMat[, -1, drop=FALSE]
rowNamesTmp <- rownames(exprVals)
colNamesTmp <- colnames(exprVals)
numExpr <- matrix(as.numeric(as.matrix(exprVals)), nrow=nrow(exprVals),
                  dimnames=list(rowNamesTmp, colNamesTmp))
numExpr <- avereps(numExpr)

# 读取目标基因列表（CSV格式，有表头）
gene_list <- read.csv(gene_list_file, header=TRUE, stringsAsFactors=FALSE)
targetGene <- as.character(gene_list[, 1])
targetGene <- targetGene[targetGene != ""]

# 筛选目标基因
matched_genes <- targetGene[targetGene %in% rownames(numExpr)]
cat("目标基因数：", length(targetGene), "，匹配到：", length(matched_genes), "\n")
exprSelected <- t(numExpr[matched_genes, , drop=FALSE])

# 从CIBERSORT结果提取免疫浸润矩阵（宽格式）
immune_cols <- setdiff(colnames(rt), c("Sample", "Group"))
immuneMat <- as.matrix(rt[, immune_cols])
immuneMat <- apply(immuneMat, 2, as.numeric)
rownames(immuneMat) <- rownames(rt)

# 取共同样本
commonSamples <- intersect(rownames(exprSelected), rownames(immuneMat))
cat("共同样本数：", length(commonSamples), "\n")
exprSelected <- exprSelected[commonSamples, , drop=FALSE]
immuneMat    <- immuneMat[commonSamples, , drop=FALSE]

# 移除标准差为0的免疫细胞列
validImmune <- immuneMat[, apply(immuneMat, 2, sd) > 0, drop=FALSE]
cat("有效免疫细胞类型数：", ncol(validImmune), "\n")

# 计算Spearman相关性
cor_mat     <- matrix(NA, nrow=ncol(exprSelected), ncol=ncol(validImmune),
                      dimnames=list(colnames(exprSelected), colnames(validImmune)))
p_val_matrix <- matrix(NA, nrow=ncol(exprSelected), ncol=ncol(validImmune),
                       dimnames=list(colnames(exprSelected), colnames(validImmune)))

for (i in seq_len(ncol(exprSelected))) {
  for (j in seq_len(ncol(validImmune))) {
    test_result <- cor.test(exprSelected[,i], validImmune[,j], method="spearman", use="pairwise.complete.obs")
    cor_mat[i,j]     <- test_result$estimate
    p_val_matrix[i,j] <- test_result$p.value
  }
}

# 创建标签矩阵（相关系数 + 显著性星号）
labels_mat <- matrix("", nrow=nrow(cor_mat), ncol=ncol(cor_mat),
                     dimnames=dimnames(cor_mat))
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    stars <- ""
    if (!is.na(p_val_matrix[i,j])) {
      if (p_val_matrix[i,j] < 0.001) stars <- "***"
      else if (p_val_matrix[i,j] < 0.01) stars <- "**"
      else if (p_val_matrix[i,j] < 0.05) stars <- "*"
    }
    cor_val <- ifelse(is.na(cor_mat[i,j]), "", round(cor_mat[i,j], 2))
    labels_mat[i,j] <- paste0(cor_val, "\n", stars)
  }
}

# 调试输出
cat("cor_mat 维度：", nrow(cor_mat), "x", ncol(cor_mat), "\n")
cat("cor_mat 行名：", paste(rownames(cor_mat), collapse=", "), "\n")
cat("labels_mat 维度：", nrow(labels_mat), "x", ncol(labels_mat), "\n")

# 绘制相关性热图
gene_num   <- nrow(cor_mat)
immune_num <- ncol(cor_mat)
pdf_width  <- max(8, 4 + 0.5 * immune_num)
pdf_height <- max(5, 3 + 0.6 * gene_num)
my_col_fun <- colorRampPalette(c("#6b8e23", "white", "#ba6262"))(100)

# 确保cor_mat转为data.frame再转回matrix，避免pheatmap的维度问题
cor_mat_plot <- as.matrix(as.data.frame(cor_mat))

pdf("Gene_Immune_Correlation_Heatmap.pdf", width=pdf_width, height=pdf_height)
pheatmap(
  cor_mat_plot, color=my_col_fun, display_numbers=labels_mat,
  number_color="black", cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize_number=7, fontsize_row=12, fontsize_col=9, border_color="grey90",
  main="Hub Genes - Immune Cell Correlation Heatmap",
  angle_col=45
)
dev.off()
cat(">> 相关性热图已保存为：Gene_Immune_Correlation_Heatmap.pdf\n")

# 导出相关性结果表
corrDF <- data.frame()
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    corrDF <- rbind(corrDF, data.frame(
      Gene=rownames(cor_mat)[i], ImmuneCell=colnames(cor_mat)[j],
      Correlation=cor_mat[i,j], P.Value=p_val_matrix[i,j],
      Significance=ifelse(p_val_matrix[i,j] < 0.05, "Yes", "No")
    ))
  }
}
write.csv(corrDF, file="Gene_Immune_Correlation_Results.csv", row.names=FALSE)
cat(">> 相关性结果表已保存为：Gene_Immune_Correlation_Results.csv\n")

# ============================================================
# 第三部分：方块气泡图（Square Bubble Chart）基因-免疫细胞相关性
# ============================================================
cat("\n==== 第三部分：方块气泡图 ====\n")

# 准备数据（直接用corrDF，免疫细胞名保持原样）
bubble_df <- corrDF[!is.na(corrDF$P.Value) & !is.na(corrDF$Correlation), ]
bubble_df$AbsCorr <- abs(bubble_df$Correlation)

# 显著性标记
bubble_df$Sig <- ifelse(bubble_df$P.Value < 0.001, "***",
                 ifelse(bubble_df$P.Value < 0.01, "**",
                 ifelse(bubble_df$P.Value < 0.05, "*", "")))

# 调试：打印数据确认
cat("bubble_df 行数：", nrow(bubble_df), "\n")
cat("基因：", paste(unique(bubble_df$Gene), collapse=", "), "\n")
cat("免疫细胞：", paste(unique(bubble_df$ImmuneCell), collapse=", "), "\n")

# 免疫细胞按名称排序
bubble_df$ImmuneCell <- factor(bubble_df$ImmuneCell,
                               levels = rev(sort(unique(bubble_df$ImmuneCell))))
bubble_df$Gene <- factor(bubble_df$Gene, levels = sort(unique(as.character(bubble_df$Gene))))

# 图尺寸
n_genes  <- length(unique(bubble_df$Gene))
n_immune <- length(unique(bubble_df$ImmuneCell))
fig_width  <- max(6, 3 + n_genes * 1.5)
fig_height <- max(8, 2 + n_immune * 0.42)

pdf("Gene_Immune_SquareBubble.pdf", width = fig_width-2, height = fig_height-4)
p_bubble <- ggplot(bubble_df, aes(x = Gene, y = ImmuneCell)) +
  geom_point(aes(size = AbsCorr, fill = Correlation),
             shape = 22, color = "grey30", stroke = 0.3) +
  scale_fill_viridis_c(option = "viridis",
                       limits = c(-max(bubble_df$AbsCorr), max(bubble_df$AbsCorr)),
                       oob = squish, name = "Correlation",
                       breaks = pretty(range(bubble_df$Correlation), n = 5)) +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  geom_text(data = bubble_df[bubble_df$Sig != "", ],
            aes(label = Sig), color = "white", size = 4, fontface = "bold", vjust = 0.5) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 1.2),
    panel.grid.major = element_line(color = "black", size = 0.6),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "right",
    legend.key.height = unit(1.8, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )
print(p_bubble)
dev.off()
cat(">> 方块气泡图已保存为：Gene_Immune_SquareBubble.pdf\n")

# ============================================================
# 第四部分：linkET相关性图（cor2.pdf 和 cor_newStyle.pdf）
# ============================================================
cat("\n==== 第四部分：linkET相关性图 ====\n")

# 准备数据：计算基因与免疫细胞的相关性（用于linkET）
corrDF_linkET <- data.frame()

for (gene in colnames(exprSelected)) {
  for (cellType in colnames(validImmune)) {
    if (sd(validImmune[, cellType]) == 0) {
      next
    }
    immuneVec <- as.numeric(validImmune[, cellType])
    geneExpr  <- as.numeric(exprSelected[, gene])
    testResult <- cor.test(immuneVec, geneExpr, method = "spearman")

    tempDF <- data.frame(
      spec = gene,
      env  = cellType,
      r    = as.numeric(testResult$estimate),
      p    = as.numeric(testResult$p.value)
    )
    corrDF_linkET <- rbind(corrDF_linkET, tempDF)
  }
}

# 根据P值判断相关性方向
corrDF_linkET$pd <- ifelse(corrDF_linkET$p < 0.05,
                           ifelse(corrDF_linkET$r > 0, "positive correlation", "negative correlation"),
                           "not significant")

# 将相关系数取绝对值
corrDF_linkET$r <- abs(corrDF_linkET$r)

# 为相关系数添加分档信息
corrDF_linkET <- corrDF_linkET %>%
  mutate(rd = cut(r,
                  breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))

# ---------- 图形风格A：Spectral配色风格（cor_newStyle.pdf）---------- #
cat("绘制图形风格A（cor_newStyle.pdf）...\n")

plotA <- qcorrplot(correlate(validImmune, method = "spearman"), type = "lower", diag = FALSE) +
  geom_tile() +
  geom_couple(aes(colour = pd, size = rd),
              data = corrDF_linkET, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) +
  scale_size_manual(values = c(0.5, 1, 2, 3)) +
  scale_colour_manual(values = c("positive correlation" = "#FF1493",
                                 "negative correlation" = "#00CED1",
                                 "not significant"  = "#999999")) +
  guides(size = guide_legend(title = "abs(Cor)", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "P-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3)) +
  labs(x = "immune infiltrating cells", y = "immune infiltrating cells",
       title = "Gene-immune infiltrating cells Correlation") +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),
        axis.text.x  = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5))

pdf(file = "cor_newStyle.pdf", width = 12, height = 7)
print(plotA)
dev.off()
cat(">> cor_newStyle.pdf 已保存\n")

# ---------- 图形风格B：自定义主题 + Pastel2配色（cor2.pdf）---------- #
cat("绘制图形风格B（cor2.pdf）...\n")

# 定义自定义主题函数
theme_cute <- function() {
  theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "#f0f0f0", size = 0.5),
      panel.grid.minor = element_line(colour = "#f0f0f0", size = 0.25),
      axis.text        = element_text(size = 10, colour = "#555555"),
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title       = element_text(size = 12, face = "bold", colour = "#555555"),
      plot.title       = element_text(size = 18, face = "bold", colour = "#d35400", hjust = 0.5)
    )
}

plotB <- qcorrplot(correlate(validImmune, method = "spearman"), type = "lower", diag = FALSE) +
  geom_tile() +
  geom_couple(aes(colour = pd, size = rd), data = corrDF_linkET, curvature = nice_curvature()) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(9, "Pastel2")),
    name = "Cell-Cell Correlation"
  ) +
  scale_size_manual(
    name   = "abs(Cor)",
    values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3),
    labels = c("< 0.2" = "< 0.2", "0.2 - 0.4" = "0.2 - 0.4", "0.4 - 0.6" = "0.4 - 0.6", ">= 0.6" = "≥ 0.6")
  ) +
  scale_colour_manual(
    name   = "p-value",
    values = c("positive correlation" = "#F28C8C",
               "negative correlation" = "#8AB8FF",
               "not significant"      = "#B2BABB"),
    labels = c("positive correlation" = "Positive",
               "negative correlation" = "Negative",
               "not significant"      = "Not significant")
  ) +
  guides(
    size   = guide_legend(title = "abs(Cor)", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "p-value", override.aes = list(size = 3), order = 1),
    fill   = guide_colorbar(title = "Cell-Cell Correlation", order = 3)
  ) +
  ggtitle("Single Gene–immune infiltrating cells Correlation") +
  theme_cute() +
  labs(x = NULL, y = NULL)

pdf(file = "cor2.pdf", width = 9, height = 7)
print(plotB)
dev.off()
cat(">> cor2.pdf 已保存\n")

# ============================================================
# 第五部分：基因-免疫细胞相关性棒棒糖图（每个基因一个图）
# ============================================================
cat("\n==== 第五部分：基因-免疫细胞相关性棒棒糖图 ====\n")

# 为每个基因生成棒棒糖图
for (gene in colnames(exprSelected)) {
  cat(sprintf("正在生成基因 %s 的棒棒糖图...\n", gene))

  # 准备该基因的相关性数据
  gene_data <- data.frame()

  for (cellType in colnames(validImmune)) {
    if (sd(validImmune[, cellType]) == 0) {
      next
    }
    immuneVec <- as.numeric(validImmune[, cellType])
    geneExpr  <- as.numeric(exprSelected[, gene])
    testResult <- cor.test(immuneVec, geneExpr, method = "spearman")

    tempDF <- data.frame(
      env = cellType,
      r   = as.numeric(testResult$estimate),
      p   = as.numeric(testResult$p.value)
    )
    gene_data <- rbind(gene_data, tempDF)
  }

  # 如果没有有效数据，跳过
  if (nrow(gene_data) == 0) {
    cat(sprintf("  警告：基因 %s 没有有效的相关性数据，跳过\n", gene))
    next
  }

  # 1. 渐变色参数
  nColors <- 500
  myGradient <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(nColors)

  # 2. 主图点颜色基于p值归一化
  p_min <- min(gene_data$p, na.rm = TRUE)
  p_max <- max(gene_data$p, na.rm = TRUE)

  if (p_max > p_min) {
    norm_p <- (gene_data$p - p_min) / (p_max - p_min)
  } else {
    norm_p <- rep(0.5, nrow(gene_data))
  }

  col_index <- round(norm_p * (nColors - 1)) + 1
  col_index[col_index > nColors] <- nColors
  col_index[col_index < 1] <- 1
  gene_data$points.color <- myGradient[col_index]

  # 3. 点大小分档
  p.cex <- seq(2.5, 5.5, length = 5)
  fcex <- function(x) {
    x <- abs(x)
    cex <- ifelse(x < 0.1, p.cex[1],
                  ifelse(x < 0.2, p.cex[2],
                         ifelse(x < 0.3, p.cex[3],
                                ifelse(x < 0.4, p.cex[4], p.cex[5]))))
    return(cex)
  }
  gene_data$points.cex <- fcex(gene_data$r)

  # 4. 排序
  gene_data <- gene_data[order(gene_data$r), ]

  # 5. 画图
  pdf_filename <- sprintf("%s_Lollipop.pdf", gene)
  pdf(file = pdf_filename, width = 9.5, height = 7)

  xlim <- ceiling(max(abs(gene_data$r)) * 10) / 10
  layout(mat = matrix(c(1, 1, 1, 1, 1, 0, 2, 0, 3, 0), nc = 2),
         width = c(8, 2.2), heights = c(1, 2, 1, 2, 1))
  par(bg = "white", las = 1, mar = c(5, 18, 2, 4), cex.axis = 1.5, cex.lab = 2)

  # 主体lollipop plot
  plot(1, type = "n", xlim = c(-xlim, xlim), ylim = c(0.5, nrow(gene_data) + 0.5),
       xlab = "Correlation Coefficient", ylab = "", yaxt = "n", yaxs = "i", axes = FALSE)
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4],
       col = "#F5F5F5", border = NA)
  grid(ny = nrow(gene_data), col = "white", lty = 1, lwd = 2)

  segments(x0 = 0, y0 = 1:nrow(gene_data), x1 = gene_data$r, y1 = 1:nrow(gene_data), lwd = 4)
  points(x = gene_data$r, y = 1:nrow(gene_data), col = gene_data$points.color,
         pch = 16, cex = gene_data$points.cex)

  text(par('usr')[1], 1:nrow(gene_data), gene_data$env, adj = 1, xpd = TRUE, cex = 1.5)
  pvalue.text <- ifelse(gene_data$p < 0.001, "<0.001", sprintf("%.03f", gene_data$p))
  redcutoff_cor <- 0
  redcutoff_pvalue <- 0.05
  text(par('usr')[2], 1:nrow(gene_data), pvalue.text, adj = 0, xpd = TRUE,
       col = ifelse(abs(gene_data$r) > redcutoff_cor & gene_data$p < redcutoff_pvalue, "red", "black"),
       cex = 1.5)
  axis(1, tick = FALSE)

  # 点大小图例
  par(mar = c(0, 4, 3, 4))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("left", legend = c(0.1, 0.2, 0.3, 0.4, 0.5), col = "black", pt.cex = p.cex,
         pch = 16, bty = "n", cex = 2, title = "abs(r)")

  # P值渐变色图注
  par(mar = c(0, 2, 4, 6), cex.axis = 1.5, cex.main = 2)
  image(
    x = 1,
    y = seq(0, 1, length = nColors),
    z = matrix(seq(0, 1, length = nColors), nrow = 1),
    col = myGradient,
    axes = FALSE,
    xlab = "", ylab = "", main = "P value"
  )
  axis(4, at = seq(0, 1, length = 5), labels = round(seq(0, 1, length = 5), 2),
       las = 2, tick = FALSE)

  dev.off()
  cat(sprintf(">> %s 已保存\n", pdf_filename))
}

cat("\n==== 全流程执行完毕！====\n")
