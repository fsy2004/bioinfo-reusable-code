# =============================================================================
# 编号       : R007
# 脚本名     : KEGG和GO.R
# 分类       : 02_GO_KEGG富集分析
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 对候选基因进行 GO/KEGG 富集分析，输出通路表和条形图、气泡图等。
# 结果图     : 条形图/柱状图；气泡图/点图
# 主要 R 包  : circlize; clusterProfiler; ComplexHeatmap; dplyr; enrichplot; ggplot2; ggpubr; KEGGREST; org.Hs.eg.db; RColorBrewer
# 整理时间   : 2026-05-10
# =============================================================================
# =================== 参数设置部分 ===================
# 1. p值过滤阈值
filter_pvalue <- 0.05         # p值过滤条件
filter_padjust <- 0.05        # 矫正p值过滤条件

# 2. 展示通路数量
max_show_num <- 30            # 最多展示的通路数量

# 3. 颜色方式
color_by <- "p.adjust"        # 默认颜色映射
if(filter_padjust > 0.05){    # 若矫正p值阈值较宽松，则用原始p值上色
  color_by <- "pvalue"
}

# 4. 文件与工作目录
data_file <- "IntersectionGenes.csv"   # 输入文件名
output_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/10.富集分析KEGG和GO" # 工作目录
output_kegg <- "KEGG.csv"              # 富集分析输出文件
setwd(output_dir)                      # 设置工作目录

# =================== 加载所需包 ===================
library(clusterProfiler)    # 富集分析包
library(org.Hs.eg.db)       # 人类基因注释包
library(enrichplot)         # 富集分析可视化
library(ggplot2)            # 画图
library(circlize)           # 环形图可视化
library(RColorBrewer)       # 颜色方案
library(dplyr)              # 数据处理
library(ggpubr)             # ggplot美化
library(KEGGREST)           # KEGG数据库接口
library(ComplexHeatmap)     # 热图绘制

# =================== 步骤1：数据读取与预处理 ===================
cat("步骤1/5：开始读取数据...\n")   # 进度条
if(!file.exists(data_file)){          # 判断输入文件是否存在
  stop("输入文件不存在，请检查路径和文件名！")
}
raw_data <- read.csv(data_file, header=TRUE, stringsAsFactors=FALSE)  # 读取CSV文件

# =================== 步骤2：基因名转换与检查 ===================
cat("步骤2/5：进行基因名到ENTREZ ID的转换...\n")
gene_symbols <- unique(as.character(raw_data$Gene))    # 提取Gene列的基因名

# 判断基因数是否为空
if(length(gene_symbols)==0){
  stop("未检测到基因名，请检查输入文件内容！")
}

# 转换基因名到ENTREZ ID
entrez_result <- mget(gene_symbols, org.Hs.egSYMBOL2EG, ifnotfound=NA)   # 用org.Hs.eg.db查ID
entrez_id_vec <- as.character(entrez_result)      # 转成字符向量
map_table <- data.frame(Symbol=gene_symbols, EntrezID=entrez_id_vec)   # 合成表格

# 去除无效ID
valid_gene_ids <- entrez_id_vec[!is.na(entrez_id_vec) & entrez_id_vec!="NA"] # 只保留有效ID

# 冗余的变量和判断（安全性代码）
if(length(valid_gene_ids) == 0){
  stop("所有基因名转换失败，无有效ENTREZ ID。")
}
redundant_tmp <- sum(is.na(entrez_id_vec)) # 冗余：统计NA数
if(redundant_tmp > 0) {cat(sprintf("有%d个基因未能转换为ID。\n", redundant_tmp))}

# =================== 步骤3：KEGG富集分析 ===================
cat("步骤3/5：正在进行KEGG富集分析...\n")

# 仅在有有效基因时分析
if(length(valid_gene_ids)>0){
  kegg_result <- enrichKEGG(gene=valid_gene_ids, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)  # 运行KEGG分析
  kegg_df <- as.data.frame(kegg_result)    # 结果转为数据框
  # 基因ID回转成基因名
  kegg_df$GeneSymbol <- sapply(kegg_df$geneID, function(x) {
    gene_list <- unlist(strsplit(x, "/"))
    paste(map_table$Symbol[match(gene_list, map_table$EntrezID)], collapse="/")
  })
  # 按阈值过滤
  filtered_kegg <- kegg_df %>% filter(pvalue < filter_pvalue & p.adjust < filter_padjust)
  
  # 如果没有显著通路，给提示
  if(nrow(filtered_kegg)==0){
    cat("未检出显著富集通路，输出原始结果。\n")
  }
  
  # 添加category和subcategory列（从kegg_result中提取）
  if(nrow(filtered_kegg) > 0) {
    # 从KEGG数据库获取pathway的分类信息
    pathway_class <- KEGGREST::keggList("pathway", "hsa")

    # 为每个通路添加分类信息
    filtered_kegg$category <- NA
    filtered_kegg$subcategory <- NA

    for(i in 1:nrow(filtered_kegg)) {
      pathway_id <- filtered_kegg$ID[i]
      tryCatch({
        pathway_info <- KEGGREST::keggGet(pathway_id)[[1]]
        if(!is.null(pathway_info$CLASS)) {
          class_parts <- strsplit(pathway_info$CLASS, "; ")[[1]]
          if(length(class_parts) >= 1) filtered_kegg$category[i] <- class_parts[1]
          if(length(class_parts) >= 2) filtered_kegg$subcategory[i] <- class_parts[2]
        }
      }, error = function(e) {
        # 如果获取失败，保持NA
      })
    }

    # 重新排列列顺序，将category和subcategory放在前面
    filtered_kegg <- filtered_kegg[, c("category", "subcategory", setdiff(names(filtered_kegg), c("category", "subcategory")))]
  }

  # 保存结果到文件
  write.table(filtered_kegg, file=output_kegg, sep=",", quote=FALSE, row.names=FALSE)
} else {
  filtered_kegg <- data.frame()
}

# =================== 步骤4：可视化分析 ===================
cat("步骤4/5：生成KEGG可视化图形...\n")

# 判断是否有显著结果用于展示
show_count <- min(nrow(filtered_kegg), max_show_num)
if(show_count == 0){
  show_count <- 1   # 至少展示1个以避免错误
}

# ----------- 1. 柱状图 -----------
pdf(file="barplot.pdf", width=8, height=7)     # 创建PDF
tryCatch({
  barplot(kegg_result, drop=TRUE, showCategory=show_count, label_format=130, color=color_by)
}, error=function(e){cat("barplot绘图出错。\n")})
dev.off()  # 关闭PDF

# ----------- 2. 气泡图 -----------
pdf(file="bubble.pdf", width=8, height=7)
tryCatch({
  dotplot(kegg_result, showCategory=show_count, orderBy="GeneRatio", label_format=130, color=color_by)
}, error=function(e){cat("dotplot绘图出错。\n")})
dev.off()

# ----------- 3. KEGG分类分析图 (B) -----------
pdf(file="KEGG_classification_analysis.pdf", width=14, height=10)
tryCatch({
  if(nrow(filtered_kegg) > 0 && "category" %in% colnames(filtered_kegg)) {
    # 直接使用CSV文件中的category列
    kegg_class <- data.frame(
      ID = filtered_kegg$ID,
      Description = filtered_kegg$Description,
      Count = filtered_kegg$Count,
      p.adjust = filtered_kegg$p.adjust,
      Category = filtered_kegg$category
    )

    # 取每个类别中最显著的通路（最多10个）
    top_pathways_by_category <- kegg_class %>%
      group_by(Category) %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      ungroup() %>%
      arrange(Category, p.adjust)

    # 按类别分组排序，确保同一颜色的条形聚集在一起
    top_pathways_by_category$Description <- factor(top_pathways_by_category$Description,
                                                   levels = top_pathways_by_category$Description)

    # 设置颜色方案
    category_colors <- c(
      "Environmental Information Processing" = "#00BA38",
      "Cellular Processes" = "#619CFF",
      "Organismal Systems" = "#F8766D",
      "Human Diseases" = "#00BFC4",
      "Metabolism" = "#C77CFF",
      "Genetic Information Processing" = "#E76BF3",
      "Others" = "#FF9999"
    )

    # 创建分类分析条形图（按类别分组显示具体通路）
    p_class <- ggplot(top_pathways_by_category, aes(x = Count, y = Description, fill = Category)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = category_colors, name = "KEGG Classification") +
      labs(x = "Count", y = "", title = "KEGG Classification Analysis") +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 11),
        # 移除网格线但保留坐标轴线
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.5)
      ) +
      geom_text(aes(label = Count), hjust = -0.1, size = 3)

    print(p_class)
  } else {
    plot.new()
    text(0.5,0.5,"没有显著通路可进行分类分析",cex=1.2)
  }
}, error=function(e){cat("KEGG分类分析图绘制出错。\n")})
dev.off()

# =================== 步骤5：GO富集分析 ===================
cat("步骤5/8：正在进行GO富集分析...\n")

# 使用相同的valid_gene_ids进行GO富集分析
if(length(valid_gene_ids) > 0) {
  go_result <- enrichGO(
    gene = valid_gene_ids,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    ont = "all",
    readable = TRUE
  )

  go_df <- as.data.frame(go_result)

  # 按阈值过滤
  filtered_go <- go_df %>% filter(pvalue < filter_pvalue & p.adjust < filter_padjust)

  # 如果没有显著GO term，给提示
  if(nrow(filtered_go) == 0) {
    cat("未检出显著富集GO term，输出原始结果。\n")
  }

  # 保存GO结果到文件
  write.table(filtered_go, file="GO_results.csv", sep=",", quote=FALSE, row.names=FALSE)

} else {
  filtered_go <- data.frame()
}

# =================== 步骤6：GO可视化分析 ===================
cat("步骤6/8：生成GO可视化图形...\n")

# 判断是否有显著结果用于展示
if(nrow(filtered_go) > 0) {

  # ----------- 1. GO柱状图（分面显示BP/CC/MF）-----------
  pdf("GO_barplot.pdf", width=8, height=10)
  tryCatch({
    bar_plot <- barplot(
      go_result,
      drop = TRUE,
      showCategory = 10,
      label_format = 50,
      split = "ONTOLOGY",
      color = color_by
    ) +
      facet_grid(ONTOLOGY ~ ., scale = 'free') +
      scale_fill_gradientn(colors = c("#FF6666", "#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
    print(bar_plot)
  }, error=function(e){cat("GO柱状图绘图出错。\n")})
  dev.off()

  # ----------- 2. GO气泡图（分面显示BP/CC/MF）-----------
  pdf("GO_bubble.pdf", width=8, height=10)
  tryCatch({
    bubble_plot <- dotplot(
      go_result,
      showCategory = 10,
      orderBy = "GeneRatio",
      label_format = 50,
      split = "ONTOLOGY",
      color = color_by
    ) +
      facet_grid(ONTOLOGY ~ ., scale = 'free') +
      scale_color_gradientn(colors = c("#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
    print(bubble_plot)
  }, error=function(e){cat("GO气泡图绘图出错。\n")})
  dev.off()

  # ----------- 3. GO分组条形图 -----------
  pdf("GO_grouped_barplot.pdf", width=11, height=8)
  tryCatch({
    top_go <- filtered_go %>% group_by(ONTOLOGY) %>% slice_head(n = 10)

    group_bar_plot <- ggbarplot(
      top_go,
      x = "Description", y = "Count", fill = "ONTOLOGY", color = "white",
      xlab = "", palette = "aaas",
      legend = "right", sort.val = "desc", sort.by.groups = TRUE,
      position = position_dodge(0.9)
    ) +
      rotate_x_text(75) +
      theme(panel.background = element_blank(),
            axis.text.x = element_text(size = 10, color = "black")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      scale_x_discrete(expand = c(0, 0)) +
      geom_text(
        aes(label = Count),
        position = position_dodge(0.9), vjust = -0.3, size = 3
      )
    print(group_bar_plot)
  }, error=function(e){cat("GO分组条形图绘图出错。\n")})
  dev.off()

  # ----------- 4. GO自定义条形图和气泡图（分面）-----------
  tryCatch({
    # 准备数据
    filtered_go_sub <- filtered_go %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

    top_go_custom <- filtered_go_sub %>%
      group_by(ONTOLOGY) %>%
      slice_min(order_by = p.adjust, n = 10) %>%
      ungroup() %>%
      arrange(ONTOLOGY, p.adjust)

    top_go_custom <- top_go_custom %>%
      group_by(ONTOLOGY) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
      ungroup()

    # GeneRatio转数值
    top_go_custom$GeneRatio_num <- sapply(top_go_custom$GeneRatio, function(x) {
      sp <- unlist(strsplit(as.character(x), "/"))
      as.numeric(sp[1]) / as.numeric(sp[2])
    })

    # 自定义条形图
    bar_colors <- brewer.pal(7, "YlOrRd")
    p_bar <- ggplot(top_go_custom, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
      geom_bar(stat = "identity", width = 0.8) +
      coord_flip() +
      facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
      scale_fill_gradientn(colors = bar_colors) +
      labs(x = "GO Term", y = "-log10(adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        strip.text = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, face = "bold")
      )
    ggsave("GO_barplot_custom.pdf", p_bar, width = 12, height = 12)

    # 自定义气泡图
    dot_colors <- brewer.pal(7, "Spectral")
    p_dot <- ggplot(top_go_custom, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
      geom_point(alpha = 0.8) +
      facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
      scale_color_gradientn(colors = dot_colors) +
      labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
        strip.text = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, face = "bold")
      )
    ggsave("GO_dotplot_custom.pdf", p_dot, width = 12, height = 12)

  }, error=function(e){cat("GO自定义图形绘图出错。\n")})

} else {
  cat("没有显著的GO富集结果，跳过GO可视化。\n")
}

# =================== 步骤7：流程结束 ===================
cat("步骤7/7：全部KEGG和GO分析与绘图完成！请查看结果文件和PDF图片。\n")

# 检查下内存
gc()
