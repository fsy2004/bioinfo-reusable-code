# =============================================================================
# 编号       : R004
# 脚本名     : 提取基因.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 从 GeneCards 检索结果中提取疾病相关基因。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : 
# 整理时间   : 2026-05-10
# =============================================================================
# =====================================
# 从GeneCards结果中提取基因
# 基于Relevance Score阈值筛选
# =====================================

# 清空工作环境
rm(list = ls())

# 设置工作路径
setwd("C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/07.genecard寻找疾病的相关靶点")

# =====================================
# 用户配置参数
# =====================================
# 只有Relevance Score >= 此值的基因才会被提取
RELEVANCE_THRESHOLD <- 6  # 可以修改此数值

# 设置要提取的基因类型（Category）
GENE_CATEGORIES <- c("Protein Coding", "RNA Gene")  # 只提取这两种类型

# 输入和输出文件
input_file <- "GeneCards-SearchResults.csv"
output_file <- "Extracted_Genes.csv"

# =====================================
# 进度显示函数
# =====================================
show_progress <- function(current, total, step_name = "") {
  percent <- round(current / total * 100)
  cat("\r", step_name, " ", percent, "%", sep = "")
  if (current == total) {
    cat(" - 完成\n")
  }
  flush.console()
}

# =====================================
# 主分析流程
# =====================================
cat("\n==================== GeneCards基因提取分析 ====================\n\n")

cat("配置参数:\n")
cat("  Relevance Score阈值: >=", RELEVANCE_THRESHOLD, "\n")
cat("  基因类型过滤:", paste(GENE_CATEGORIES, collapse = ", "), "\n")
cat("  输入文件:", input_file, "\n")
cat("  输出文件:", output_file, "\n\n")

# =====================================
# 步骤1: 加载数据
# =====================================
cat("步骤1/4: 加载GeneCards数据\n")

if (!file.exists(input_file)) {
  stop("错误: 未找到输入文件! 请检查文件路径。")
}

# 显示加载进度
for (i in 1:100) {
  show_progress(i, 100, "加载数据:")
  if (i == 50) {
    data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
  }
  Sys.sleep(0.008)
}

cat("已加载记录数:", nrow(data), "\n")
cat("列名:", paste(colnames(data), collapse = ", "), "\n\n")

# =====================================
# 步骤2: 按Relevance Score过滤
# =====================================
cat("步骤2/4: 按Relevance Score过滤数据\n")

# 清理列名（移除空格）
colnames(data) <- gsub(" ", ".", colnames(data))

# 检查Relevance.score列是否存在
if (!"Relevance.score" %in% colnames(data)) {
  stop("错误: 数据中未找到'Relevance score'列!")
}

# 将Relevance.score转换为数值型（以防读取为字符型）
data$Relevance.score <- as.numeric(data$Relevance.score)

# 显示过滤进度
for (i in 1:100) {
  show_progress(i, 100, "过滤数据:")
  Sys.sleep(0.008)
}

# 根据阈值过滤数据
filtered_data <- data[data$Relevance.score >= RELEVANCE_THRESHOLD, ]

cat("按Relevance Score过滤:\n")
cat("  过滤前记录数:", nrow(data), "\n")
cat("  过滤后记录数:", nrow(filtered_data), "\n")
cat("  已过滤掉:", nrow(data) - nrow(filtered_data), "条记录\n\n")

# =====================================
# 步骤2.5: 按基因类型过滤
# =====================================
cat("步骤2.5: 按基因类型过滤\n")

# 检查Category列是否存在
if (!"Category" %in% colnames(filtered_data)) {
  stop("错误: 数据中未找到'Category'列!")
}

# 显示过滤前的类型分布
cat("过滤前的基因类型分布:\n")
print(table(filtered_data$Category))
cat("\n")

# 按基因类型过滤
before_category_filter <- nrow(filtered_data)
filtered_data <- filtered_data[filtered_data$Category %in% GENE_CATEGORIES, ]

# 显示类型过滤进度
for (i in 1:100) {
  show_progress(i, 100, "类型过滤:")
  Sys.sleep(0.008)
}

cat("按基因类型过滤:\n")
cat("  过滤前记录数:", before_category_filter, "\n")
cat("  过滤后记录数:", nrow(filtered_data), "\n")
cat("  已过滤掉:", before_category_filter - nrow(filtered_data), "条记录\n\n")

# =====================================
# 步骤3: 提取基因信息
# =====================================
cat("步骤3/4: 提取基因信息\n")

# 显示提取进度
for (i in 1:100) {
  show_progress(i, 100, "提取基因:")
  Sys.sleep(0.008)
}

# 选择相关列
gene_data <- filtered_data[, c("Gene.Symbol", "Description", "Category",
                                "Uniprot.ID", "Relevance.score")]

# 按Relevance Score降序排列
gene_data <- gene_data[order(-gene_data$Relevance.score), ]

# 重置行名
rownames(gene_data) <- NULL

cat("已提取基因数:", nrow(gene_data), "\n")
cat("Relevance Score范围:",
    round(min(gene_data$Relevance.score), 2), "-",
    round(max(gene_data$Relevance.score), 2), "\n\n")

# =====================================
# 步骤4: 保存结果
# =====================================
cat("步骤4/4: 保存结果\n")

# 显示保存进度
for (i in 1:100) {
  show_progress(i, 100, "保存文件:")
  if (i == 50) {
    write.csv(gene_data, file = output_file, row.names = FALSE)
  }
  Sys.sleep(0.008)
}

cat("结果已保存到:", output_file, "\n")
cat("文件大小:", round(file.size(output_file)/1024, 2), "KB\n\n")

# =====================================
# 统计摘要
# =====================================
cat("==================== 统计摘要 ====================\n")
cat("配置参数:\n")
cat("  Relevance Score阈值:             >=", RELEVANCE_THRESHOLD, "\n")
cat("  基因类型过滤:                    ", paste(GENE_CATEGORIES, collapse = ", "), "\n")
cat("\n结果:\n")
cat("  输入总记录数:                    ", nrow(data), "\n")
cat("  按Score过滤后:                   ", before_category_filter, "\n")
cat("  按类型过滤后:                    ", nrow(gene_data), "\n")
cat("  最终提取率:                      ", round(nrow(gene_data)/nrow(data)*100, 2), "%\n")
cat("\nRelevance Score统计:\n")
cat("  最小值:                          ", round(min(gene_data$Relevance.score), 2), "\n")
cat("  最大值:                          ", round(max(gene_data$Relevance.score), 2), "\n")
cat("  平均值:                          ", round(mean(gene_data$Relevance.score), 2), "\n")
cat("  中位数:                          ", round(median(gene_data$Relevance.score), 2), "\n")
cat("==================================================\n\n")

# =====================================
# 显示样本结果
# =====================================
cat("Relevance Score前10位的基因:\n")
cat("------------------------------------------------------------\n")
print(head(gene_data[, c("Gene.Symbol", "Relevance.score", "Description")], 10))
cat("\n")

# =====================================
# 基因类别分布
# =====================================
if ("Category" %in% colnames(gene_data)) {
  cat("基因类别分布:\n")
  cat("------------------------------------------------------------\n")
  category_table <- table(gene_data$Category)
  category_df <- data.frame(
    Category = names(category_table),
    Count = as.numeric(category_table),
    Percentage = round(as.numeric(category_table) / nrow(gene_data) * 100, 2)
  )
  print(category_df)
  cat("\n")
}

cat("==================== 分析完成! ====================\n")
cat("输出文件:", output_file, "\n")
cat("基因类型: 仅包含", paste(GENE_CATEGORIES, collapse = "和"), "\n")
cat("如需修改Relevance Score阈值或基因类型，请修改脚本顶部的\n")
cat("'RELEVANCE_THRESHOLD'和'GENE_CATEGORIES'变量。\n")
cat("===================================================\n\n")

# =====================================
# 创建仅包含基因符号的文件
# =====================================
gene_symbols_only <- data.frame(Gene.Symbol = gene_data$Gene.Symbol)
gene_symbol_file <- "GeneCards_Gene.csv"
write.csv(gene_symbols_only, file = gene_symbol_file, row.names = FALSE)
