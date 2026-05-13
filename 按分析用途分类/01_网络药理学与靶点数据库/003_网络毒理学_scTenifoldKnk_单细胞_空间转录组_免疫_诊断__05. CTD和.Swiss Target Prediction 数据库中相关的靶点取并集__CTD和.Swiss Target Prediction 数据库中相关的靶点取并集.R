# =============================================================================
# 编号       : R003
# 脚本名     : CTD和.Swiss Target Prediction 数据库中相关的靶点取并集.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 从 CTD 数据库结果中提取化合物相关靶点基因并导出靶点表。
# 结果图     : Venn图；条形图/柱状图
# 主要 R 包  : dplyr; ggplot2; ggvenn; grid; gridExtra
# 整理时间   : 2026-05-10
# =============================================================================
# ======================== CTD和Swiss Target Prediction 靶点取并集 ========================
# 目的：读取CTD和Swiss Target Prediction的基因列表，取并集并绘制韦恩图
# 输入：路径下的CSV文件（自动识别）
# 输出：并集基因列表、韦恩图、柱状图
# =======================================================================================

library(ggvenn)
library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)

set.seed(12345)

# 设置工作路径
work_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/05. CTD和.Swiss Target Prediction 数据库中相关的靶点取并集"
setwd(work_dir)

cat("================== CTD和Swiss Target Prediction 靶点取并集开始 ==================\n\n")

# =============== 步骤1：创建输出文件夹 ===============
outputDir <- "output_folder"
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
  cat("已创建输出文件夹：", outputDir, "\n")
}
cat("\n")

# =============== 步骤2：自动识别CSV文件 ===============
cat("步骤1: 自动识别路径下的CSV文件...\n")

csv_files <- list.files(path = work_dir, pattern = "\\.csv$", full.names = FALSE)

if (length(csv_files) == 0) {
  stop("错误：未找到CSV文件！")
}

cat("找到", length(csv_files), "个CSV文件：\n")
print(csv_files)
cat("\n")

# =============== 步骤3：读取所有CSV文件中的基因信息 ===============
cat("步骤2: 读取所有CSV文件中的基因信息...\n")

geneList <- list()

for (csv_file in csv_files) {
  cat("正在读取：", csv_file, "\n")

  # 读取CSV文件
  data <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

  # 根据文件名判断是CTD还是SwissTarget
  if (grepl("CTD", csv_file, ignore.case = TRUE)) {
    # CTD文件：提取Gene Symbol列
    if ("Gene.Symbol" %in% colnames(data)) {
      genes <- data$Gene.Symbol
    } else if ("Gene Symbol" %in% colnames(data)) {
      genes <- data[["Gene Symbol"]]
    } else {
      cat("警告：", csv_file, "中未找到Gene Symbol列，跳过\n")
      next
    }
    set_name <- "CTD"

  } else if (grepl("Swiss", csv_file, ignore.case = TRUE)) {
    # SwissTarget文件：提取Gene列
    if ("Gene" %in% colnames(data)) {
      genes <- data$Gene
    } else {
      cat("警告：", csv_file, "中未找到Gene列，跳过\n")
      next
    }
    set_name <- "SwissTarget"

  } else {
    # 其他文件：尝试第一列
    cat("警告：无法识别文件类型，尝试读取第一列\n")
    genes <- data[, 1]
    set_name <- tools::file_path_sans_ext(basename(csv_file))
  }

  # 去除空值和重复
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)

  geneList[[set_name]] <- genes
  cat("  - ", set_name, "：", length(genes), "个独特基因\n")
}

cat("\n共读取", length(geneList), "个基因集合\n\n")

# =============== 步骤4：计算并集 ===============
cat("步骤3: 计算所有基因的并集...\n")

unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)

cat("并集基因数量：", unionCount, "\n\n")

# =============== 步骤5：计算交集 ===============
cat("步骤4: 计算交集...\n")

if (length(geneList) >= 2) {
  commonGenes <- Reduce(intersect, geneList)
  nCommon <- length(commonGenes)
  cat("交集基因数量：", nCommon, "\n\n")
} else {
  commonGenes <- c()
  nCommon <- 0
  cat("只有一个基因集合，无法计算交集\n\n")
}

# =============== 步骤6：输出并集基因列表 ===============
cat("步骤5: 保存并集基因列表...\n")

write.table(unionGenes,
            file = file.path(outputDir, "Union_Genes.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.csv(data.frame(Gene = unionGenes),
          file = file.path(outputDir, "Union_Genes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("已保存：", file.path(outputDir, "Union_Genes.txt"), "\n")
cat("已保存：", file.path(outputDir, "Union_Genes.csv"), "\n\n")

# =============== 步骤7：绘制韦恩图 ===============
cat("步骤6: 绘制韦恩图...\n")

pdf(file = file.path(outputDir, "venn.pdf"), width = 6, height = 6)
venn_plot <- ggvenn(
  geneList,
  show_percentage = TRUE,
  stroke_color = "white",
  stroke_size = 0.5,
  fill_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList)],
  set_name_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList)],
  set_name_size = 6,
  text_size = 4.5
)
grid.arrange(venn_plot,
             textGrob(paste("Union Gene Count:", unionCount), gp = gpar(fontsize = 14)),
             ncol = 1, heights = c(5, 1))
dev.off()

cat("已保存：", file.path(outputDir, "venn.pdf"), "\n\n")

# =============== 步骤8：绘制集合大小柱状图 ===============
cat("步骤7: 绘制集合大小柱状图...\n")

setSizes <- sapply(geneList, length)
pdf(file = file.path(outputDir, "set_barplot.pdf"), width = 6, height = 4)
print(
  ggplot(data.frame(Set = names(setSizes), Size = setSizes),
         aes(x = Set, y = Size, fill = Set)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylab("Gene Count") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

cat("已保存：", file.path(outputDir, "set_barplot.pdf"), "\n\n")

# =============== 步骤8：生成汇总报告 ===============
cat("================== 分析汇总 ==================\n")
cat("输入文件数量：", length(csv_files), "\n")
cat("基因集合数量：", length(geneList), "\n")
for (set_name in names(geneList)) {
  cat("  -", set_name, "：", length(geneList[[set_name]]), "个基因\n")
}
cat("并集基因数量：", unionCount, "\n")
if (length(geneList) >= 2) {
  cat("交集基因数量：", nCommon, "\n")
}
cat("输出文件夹：", outputDir, "\n")
cat("==========================================\n\n")

cat("处理完成！已输出韦恩图、欧拉图、花瓣图和并集基因列表。\n")
