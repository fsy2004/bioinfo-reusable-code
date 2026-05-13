# =============================================================================
# 编号       : R006
# 脚本名     : venn交集.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 计算多个基因集合交集/两两交集，并生成 Venn/集合图与交集基因表。
# 结果图     : Venn图
# 主要 R 包  : ggplot2; ggvenn; viridis
# 整理时间   : 2026-05-10
# =============================================================================
# ======================== 疾病与化合物靶点取交集 ========================
# 目的：读取CSV格式的基因列表，计算交集并绘制韦恩图
# 输入：路径下的CSV文件（包含表头）
# 输出：交集基因列表、韦恩图
# ========================================================================

library(ggvenn)
library(viridis)
library(ggplot2)

set.seed(12345)

# 设置工作目录
work_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/09.疾病与化合物靶点取交集"
setwd(work_dir)

cat("================== 疾病与化合物靶点取交集开始 ==================\n\n")

# 创建输出文件夹
output_folder <- "output_folder"
if (!dir.exists(output_folder)) dir.create(output_folder)

# =============== 步骤1：读取所有CSV文件 ===============
cat("步骤1: 读取所有CSV文件...\n")

csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("错误：未找到CSV文件！")
}

gene_list <- list()

for (file in csv_files) {
  # 获取文件名（不带后缀）
  file_name <- tools::file_path_sans_ext(basename(file))

  # 读取CSV文件（包含表头）
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

  # 取第一列作为基因
  genes <- as.character(data[, 1])
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)

  # 将基因数据添加到gene_list中
  gene_list[[file_name]] <- genes

  # 打印日志
  cat("  -", file_name, "：", length(genes), "个基因\n")
}

cat("\n共读取", length(gene_list), "个基因集合\n\n")

# =============== 步骤2：计算全局交集 ===============
cat("步骤2: 计算全局交集...\n")

intersect_genes <- Reduce(intersect, gene_list)
cat("全局交集基因数：", length(intersect_genes), "\n\n")

# =============== 步骤3：保存全局交集基因 ===============
cat("步骤3: 保存全局交集基因...\n")

global_intersect_df <- data.frame(Gene = intersect_genes, stringsAsFactors = FALSE)
write.csv(global_intersect_df,
          file = file.path(output_folder, "IntersectionGenes.csv"),
          row.names = FALSE)

cat("已保存：", file.path(output_folder, "IntersectionGenes.csv"), "\n\n")

# =============== 步骤4：计算两两交集 ===============
cat("步骤4: 计算每两个文件之间的交集...\n")

pairwise_intersections <- data.frame(
  File1 = character(),
  File2 = character(),
  Intersection_Count = integer(),
  Intersection_Genes = character(),
  stringsAsFactors = FALSE
)

files <- names(gene_list)
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    file1 <- files[i]
    file2 <- files[j]

    common_genes <- intersect(gene_list[[file1]], gene_list[[file2]])

    pairwise_intersections <- rbind(pairwise_intersections,
                                    data.frame(
                                      File1 = file1,
                                      File2 = file2,
                                      Intersection_Count = length(common_genes),
                                      Intersection_Genes = paste(common_genes, collapse = ";"),
                                      stringsAsFactors = FALSE
                                    ))
  }
}

write.csv(pairwise_intersections,
          file = file.path(output_folder, "PairwiseIntersectionGenes.csv"),
          row.names = FALSE)

cat("已保存：", file.path(output_folder, "PairwiseIntersectionGenes.csv"), "\n\n")

# =============== 步骤5：绘制韦恩图 ===============
cat("步骤5: 绘制韦恩图...\n")

# 使用viridis调色板生成颜色
nColors <- length(gene_list)
myColors <- viridis(nColors)

pdf(file = file.path(output_folder, "venn.pdf"), width = 10, height = 10)

ggvenn(gene_list,
       show_percentage = TRUE,
       stroke_color = "white",
       stroke_size = 1.5,
       fill_color = myColors,
       set_name_color = "black",
       set_name_size = 8,
       text_size = 6,
       text_color = "black")

dev.off()

cat("已保存：", file.path(output_folder, "venn.pdf"), "\n\n")

# =============== 步骤6：生成汇总报告 ===============
cat("================== 分析汇总 ==================\n")
cat("输入文件数量：", length(csv_files), "\n")
for (set_name in names(gene_list)) {
  cat("  -", set_name, "：", length(gene_list[[set_name]]), "个基因\n")
}
cat("全局交集基因数：", length(intersect_genes), "\n")
cat("输出文件夹：", output_folder, "\n")
cat("==========================================\n\n")

cat("处理完成！\n")
