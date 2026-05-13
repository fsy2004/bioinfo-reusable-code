# =============================================================================
# 编号       : R011
# 脚本名     : venn交集.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 计算多个基因集合交集/两两交集，并生成 Venn/集合图与交集基因表。
# 结果图     : Venn图；UpSet图
# 主要 R 包  : ComplexUpset; dplyr; ggplot2; ggvenn; RColorBrewer; viridis
# 整理时间   : 2026-05-10
# =============================================================================
library(ggvenn)
library(RColorBrewer)
library(viridis)
library(ComplexUpset)
library(ggplot2)
library(dplyr)

# 1. 设置工作目录
setwd("C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/15.差异基因和药物分子基因取交集")

# 2. 创建输出文件夹
output_folder <- "output_folder"
if (!dir.exists(output_folder)) dir.create(output_folder)

# 3. 读取差异基因文件（DE_significant_genes.csv）
de_file <- "DE_significant_genes.csv"
de_data <- read.csv(de_file, header=TRUE, stringsAsFactors=FALSE)
cat("差异基因总数：", nrow(de_data), "\n")

# 4. 读取药物靶基因文件（CTD_Target_Genes.csv）
ctd_file <- "CTD_Target_Genes.csv"
ctd_data <- read.csv(ctd_file, header=TRUE, stringsAsFactors=FALSE)
cat("药物靶基因总数：", nrow(ctd_data), "\n")

# 5. 处理CTD数据，提取基因和作用方式
ctd_data$Gene.Symbol <- trimws(ctd_data$Gene.Symbol)
ctd_data$Interaction.Actions <- trimws(ctd_data$Interaction.Actions)

# 判断药物作用方向：increases^expression为上调，decreases^expression为下调
ctd_data$Drug_Direction <- ifelse(grepl("increases\\^expression", ctd_data$Interaction.Actions), "Up",
                                   ifelse(grepl("decreases\\^expression", ctd_data$Interaction.Actions), "Down", "Other"))

# 6. 取DE基因和CTD靶基因的交集
de_genes <- trimws(de_data$Gene)
ctd_genes <- unique(trimws(ctd_data$Gene.Symbol))
common_genes <- intersect(de_genes, ctd_genes)
cat("差异基因与药物靶基因的交集数：", length(common_genes), "\n")

# 7. 筛选表达方向一致的基因
# 提取差异基因的调控方向
de_regulation <- de_data %>% select(Gene, Regulation)
de_regulation$Gene <- trimws(de_regulation$Gene)

# 提取CTD基因的药物作用方向（去重，优先保留有明确方向的）
ctd_regulation <- ctd_data %>%
  filter(Drug_Direction %in% c("Up", "Down")) %>%
  select(Gene.Symbol, Drug_Direction, Interaction.Actions) %>%
  distinct(Gene.Symbol, .keep_all = TRUE)
colnames(ctd_regulation) <- c("Gene", "CTD_Direction", "CTD_Interaction_Actions")

# 合并两个数据框
merged_data <- merge(de_regulation, ctd_regulation, by = "Gene")

# 筛选方向一致的基因（差异基因上调且药物上调，或差异基因下调且药物下调）
consistent_genes <- merged_data %>%
  filter((Regulation == "Up" & CTD_Direction == "Up") |
         (Regulation == "Down" & CTD_Direction == "Down"))

cat("表达方向一致的基因数：", nrow(consistent_genes), "\n")
cat("  - 上调一致：", sum(consistent_genes$Regulation == "Up"), "\n")
cat("  - 下调一致：", sum(consistent_genes$Regulation == "Down"), "\n")

# 分离上调和下调的基因
consistent_up <- consistent_genes %>% filter(Regulation == "Up")
consistent_down <- consistent_genes %>% filter(Regulation == "Down")

# 保存表达方向一致的基因列表（合并表格）
write.csv(consistent_genes,
          file = file.path(output_folder, "Consistent_Direction_Genes.csv"),
          row.names = FALSE,
          quote = FALSE)
cat("表达方向一致的基因已保存至：", file.path(output_folder, "Consistent_Direction_Genes.csv"), "\n")

# 8. 准备绘制Venn图的数据（区分上调和下调）
# 上调基因
gene_list_up <- list()
gene_list_up[["DE Up"]] <- de_genes[de_data$Regulation == "Up"]
gene_list_up[["CTD Up"]] <- ctd_genes[ctd_genes %in% ctd_regulation$Gene[ctd_regulation$CTD_Direction == "Up"]]

# 下调基因
gene_list_down <- list()
gene_list_down[["DE Down"]] <- de_genes[de_data$Regulation == "Down"]
gene_list_down[["CTD Down"]] <- ctd_genes[ctd_genes %in% ctd_regulation$Gene[ctd_regulation$CTD_Direction == "Down"]]

# 9. 利用 viridis 调色板生成颜色
myColors_up <- viridis(2)
myColors_down <- viridis(2)

# 10. 绘制两个Venn图并保存为PDF（一个文件包含两个图）
pdf(file = file.path(output_folder, "venn.pdf"), width = 16, height = 8)

par(mfrow = c(1, 2))  # 设置1行2列的布局

# 绘制上调基因的Venn图
print(ggvenn(gene_list_up, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 1.5,
       fill_color = myColors_up,
       set_name_color = "black",
       set_name_size = 6,
       text_size = 4,
       text_color = "black") +
  ggtitle("Up-regulated Genes"))

# 绘制下调基因的Venn图
print(ggvenn(gene_list_down, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 1.5,
       fill_color = myColors_down,
       set_name_color = "black",
       set_name_size = 6,
       text_size = 4,
       text_color = "black") +
  ggtitle("Down-regulated Genes"))

dev.off()

cat("Venn图保存至：", file.path(output_folder, "venn.pdf"), "\n")

# 11. 读取化合物和疾病交集基因文件（IntersectionGenes.csv）
intersection_file <- "IntersectionGenes.csv"
if (file.exists(intersection_file)) {
  intersection_data <- read.csv(intersection_file, header=TRUE, stringsAsFactors=FALSE)
  intersection_genes <- trimws(intersection_data$Gene)
  cat("\n化合物和疾病交集基因数：", length(intersection_genes), "\n")
} else {
  stop("错误：未找到IntersectionGenes.csv文件")
}

# 12. 将表达方向一致的基因与化合物-疾病交集基因取交集
final_genes <- intersect(consistent_genes$Gene, intersection_genes)
cat("最终三方交集基因数：", length(final_genes), "\n")

# 13. 保存最终交集基因列表
final_genes_df <- consistent_genes %>%
  filter(Gene %in% final_genes)

write.csv(final_genes_df,
          file = file.path(output_folder, "Final_Intersection_Genes.csv"),
          row.names = FALSE,
          quote = FALSE)
cat("最终交集基因已保存至：", file.path(output_folder, "Final_Intersection_Genes.csv"), "\n")

# 14. 绘制第二个Venn图（两个圈：方向一致基因 vs 化合物-疾病基因）
gene_list2 <- list()
gene_list2[["DE(CTD Consistent)"]] <- consistent_genes$Gene
gene_list2[["Compound-Disease"]] <- intersection_genes

# 利用 viridis 调色板生成颜色
nColors2 <- length(gene_list2)
myColors2 <- viridis(nColors2)

# 绘制第二个Venn图并保存为PDF
pdf(file = file.path(output_folder, "venn_final.pdf"), width = 10, height = 10)

ggvenn(gene_list2, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 1.5,
       fill_color = myColors2,
       set_name_color = "black",
       set_name_size = 6,
       text_size = 4,
       text_color = "black")

dev.off()

cat("第二个Venn图保存至：", file.path(output_folder, "venn_final.pdf"), "\n")

# 15. 计算方向一致基因与化合物-疾病基因的交集（这就是最终三方交集）
# final_genes 已经在步骤12中计算过了，这里不需要重复计算

# 16. 保存最终交集基因列表（已在步骤13完成）
# Final_Intersection_Genes.csv 就是 Consistent_Direction_Genes 与 IntersectionGenes 的交集

# 17. 输出统计摘要
cat("\n========== 统计摘要 ==========\n")
cat("1. 差异基因总数：", length(de_genes), "\n")
cat("   - 上调：", sum(de_data$Regulation == "Up"), "\n")
cat("   - 下调：", sum(de_data$Regulation == "Down"), "\n")
cat("2. 药物靶基因总数：", length(ctd_genes), "\n")
cat("3. 差异基因与药物靶基因交集：", length(common_genes), "\n")
cat("4. 表达方向一致的基因数：", nrow(consistent_genes), "\n")
cat("   - 上调一致：", nrow(consistent_up), "\n")
cat("   - 下调一致：", nrow(consistent_down), "\n")
cat("5. 化合物-疾病交集基因数：", length(intersection_genes), "\n")
cat("6. 最终三方交集基因数（方向一致 ∩ 化合物-疾病）：", length(final_genes), "\n")
cat("   - 上调：", sum(final_genes_df$Regulation == "Up"), "\n")
cat("   - 下调：", sum(final_genes_df$Regulation == "Down"), "\n")
cat("==============================\n")
