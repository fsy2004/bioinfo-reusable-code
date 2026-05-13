# =============================================================================
# 编号       : R002
# 脚本名     : Swiss Target Prediction.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 整理 SwissTargetPrediction 预测结果，提取药物/化合物潜在靶点。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : dplyr; tidyr
# 整理时间   : 2026-05-10
# =============================================================================
# ======================== Swiss Target Prediction 靶基因提取 ========================
# 目的：从Swiss Target Prediction数据库导出的CSV文件中提取化合物的靶基因
# 输入：SwissTarget开头的CSV文件（自动识别）
# 输出基因列表（Common name列）
# ===================================================================================

library(dplyr)
library(tidyr)

set.seed(12345)

# 设置工作路径
work_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/04.Swiss Target Prediction"
setwd(work_dir)

cat("================== Swiss Target Prediction 靶基因提取开始 ==================\n\n")

# =============== 步骤1：自动识别SwissTarget开头的CSV文件 ===============
cat("步骤1: 自动识别SwissTarget开头的CSV文件...\n")

swiss_files <- list.files(path = work_dir, pattern = "^SwissTarget.*\\.csv$",
                          full.names = FALSE, ignore.case = TRUE)

if (length(swiss_files) == 0) {
  stop("错误：未找到SwissTarget开头的CSV文件！")
} else if (length(swiss_files) > 1) {
  cat("找到多个SwissTarget文件：\n")
  print(swiss_files)
  cat("使用第一个文件：", swiss_files[1], "\n")
  swiss_file <- swiss_files[1]
} else {
  cat("找到SwissTarget文件：", swiss_files[1], "\n")
  swiss_file <- swiss_files[1]
}
cat("\n")

# =============== 步骤2：读取数据 ===============
cat("步骤2: 读取Swiss Target Prediction数据...\n")

swiss_data <- read.csv(swiss_file, header = TRUE, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8")

cat("数据维度:", nrow(swiss_data), "条记录 x", ncol(swiss_data), "列\n")
cat("列名：", paste(colnames(swiss_data), collapse = ", "), "\n\n")

# =============== 步骤3：提取Common name列并处理多基因拆分 ===============
cat("步骤3: 提取Common name列，拆分多基因并去重...\n")

if (!("Common.name" %in% colnames(swiss_data))) {
  stop("错误：未找到'Common.name'列！")
}

# 有些行包含多个基因（空格分隔，如 "CCND1 CDK4"），需要拆分
gene_list <- swiss_data %>%
  select(Gene = Common.name) %>%
  filter(!is.na(Gene) & Gene != "") %>%
  mutate(Gene = strsplit(Gene, "\\s+")) %>%
  unnest(Gene) %>%
  filter(Gene != "") %>%
  distinct(Gene) %>%
  arrange(Gene)

cat("原始记录数:", nrow(swiss_data), "\n")
cat("去重后独特基因数:", nrow(gene_list), "\n\n")

cat("基因列表预览（前20个）：\n")
print(head(gene_list, 20))
cat("\n")

# =============== 步骤4：保存结果 ===============
cat("步骤4: 保存结果文件...\n")

output_file <- "SwissTarget_Genes.csv"

write.csv(gene_list, file = output_file, row.names = FALSE, fileEncoding = "UTF-8")

cat("已保存文件：", output_file, "\n\n")

# =============== 汇总 ===============
cat("================== 分析汇总 ==================\n")
cat("输入文件:", swiss_file, "\n")
cat("输出文件:", output_file, "\n")
cat("原始记录数:", nrow(swiss_data), "\n")
cat("去重后基因数:", nrow(gene_list), "\n")
cat("==========================================\n\n")

cat("分析完成！\n")
