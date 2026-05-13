# =============================================================================
# 编号       : R001
# 脚本名     : CTD中查询与某个化合物相关的靶点基因.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 从 CTD 数据库结果中提取化合物相关靶点基因并导出靶点表。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : dplyr; stringr
# 整理时间   : 2026-05-10
# =============================================================================
# ======================== CTD化合物靶点基因提取 ========================
# 目的：从CTD数据库导出的CSV文件中提取化合物相关的靶点基因
# ======================================================================

library(dplyr)
library(stringr)

set.seed(12345)

# 设置工作路径
work_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/03.CTD数据库"
setwd(work_dir)

cat("================== CTD化合物靶点基因提取开始 ==================\n\n")

# =============== 步骤1：自动识别CTD开头的CSV文件 ===============
cat("步骤1: 自动识别CTD开头的CSV文件...\n")

# 列出所有CTD开头的CSV文件
ctd_files <- list.files(path = work_dir, pattern = "^CTD.*\\.csv$", full.names = FALSE)

if (length(ctd_files) == 0) {
  stop("错误：未找到CTD开头的CSV文件！")
} else if (length(ctd_files) > 1) {
  cat("找到多个CTD文件：\n")
  print(ctd_files)
  cat("\n使用第一个文件：", ctd_files[1], "\n")
  ctd_file <- ctd_files[1]
} else {
  cat("找到CTD文件：", ctd_files[1], "\n")
  ctd_file <- ctd_files[1]
}

cat("\n")

# =============== 步骤2：读取CTD数据 ===============
cat("步骤2: 读取CTD数据...\n")

# 读取CSV文件
ctd_data <- read.csv(ctd_file, header = TRUE, stringsAsFactors = FALSE,
                     fileEncoding = "UTF-8")

cat("CTD数据维度:", nrow(ctd_data), "条记录 x", ncol(ctd_data), "列\n")

# 查看列名
cat("CTD数据列名：\n")
print(colnames(ctd_data))
cat("\n")

# 查看前几行数据
cat("CTD数据预览（前5行）：\n")
print(head(ctd_data, 5))
cat("\n")

# =============== 步骤3：提取Gene Symbol和Interaction Actions ===============
cat("步骤3: 提取Gene Symbol和Interaction Actions列...\n")

# 检查必需的列是否存在
if (!("Gene.Symbol" %in% colnames(ctd_data))) {
  stop("错误：未找到'Gene.Symbol'列！")
}

if (!("Interaction.Actions" %in% colnames(ctd_data))) {
  stop("错误：未找到'Interaction.Actions'列！")
}

# 提取两列数据
result_data <- ctd_data %>%
  select(`Gene Symbol` = Gene.Symbol,
         `Interaction Actions` = Interaction.Actions) %>%
  # 去除Gene Symbol为空的行
  filter(!is.na(`Gene Symbol`) & `Gene Symbol` != "")

cat("提取后的数据维度:", nrow(result_data), "条记录\n")
cat("独特基因数量:", length(unique(result_data$`Gene Symbol`)), "\n\n")

# 查看提取结果
cat("提取结果预览（前10行）：\n")
print(head(result_data, 10))
cat("\n")

# =============== 步骤4：统计信息 ===============
cat("步骤4: 统计Interaction Actions分布...\n")

# 统计不同Interaction Actions的数量
action_summary <- result_data %>%
  group_by(`Interaction Actions`) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

cat("Interaction Actions类型数量:", nrow(action_summary), "\n")
cat("Top 10 Interaction Actions：\n")
print(head(action_summary, 10))
cat("\n")

# 统计每个基因的记录数
gene_summary <- result_data %>%
  group_by(`Gene Symbol`) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

cat("Top 10 基因（按记录数）：\n")
print(head(gene_summary, 10))
cat("\n")

# =============== 步骤5：保存结果 ===============
cat("步骤5: 保存结果文件...\n")

# 生成输出文件名（基于输入文件名）
output_file <- "CTD_Target_Genes.csv"

# 保存结果
write.csv(result_data,
          file = output_file,
          row.names = FALSE,
          fileEncoding = "UTF-8")

cat("已保存文件：", output_file, "\n\n")

# =============== 步骤6：生成汇总报告 ===============
cat("================== 分析汇总 ==================\n")
cat("输入文件:", ctd_file, "\n")
cat("输出文件:", output_file, "\n")
cat("总记录数:", nrow(result_data), "\n")
cat("独特基因数:", length(unique(result_data$`Gene Symbol`)), "\n")
cat("Interaction Actions类型数:", nrow(action_summary), "\n")
cat("==========================================\n\n")

cat("分析完成！\n")
