# =============================================================================
# 编号       : R012
# 脚本名     : 机器学习LASSO回归筛选基因.R
# 分类       : 04_机器学习筛选特征基因
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 使用 LASSO 回归筛选特征基因并输出模型相关可视化。
# 结果图     : 气泡图/点图
# 主要 R 包  : cowplot; ggplot2; ggrepel; glmnet; RColorBrewer; reshape2
# 整理时间   : 2026-05-10
# =============================================================================

# 加载必要的包
library(glmnet)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
###########################################
#1. 设置工作目录并读取数据
###########################################

working_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/16.lasso回归筛选特征基因"  # 工作目录
inputFile <- file.path(working_dir, "Sample Type Matrix.csv")  # 基因表达数据文件
geneListFile <- file.path(working_dir, "Final_Intersection_Genes.csv")  # 要分析的基因列表文件

###########################################
# 2. 初始化进度条及步骤计数
###########################################
total_steps <- 15                           # 定义总步骤数
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)  # 初始化文本进度条
current_step <- 0                           # 初始化当前步骤计数

###########################################
# 3. 设置随机数种子，保证结果可重复
###########################################
current_step <- current_step + 1            # 步骤 1：更新步骤计数
setTxtProgressBar(pb, current_step)          # 更新进度条
set.seed(123456)                               # 设置随机种子

###########################################
# 4. 读取数据
###########################################
current_step <- current_step + 1            # 步骤 2：更新步骤计数
setTxtProgressBar(pb, current_step)          # 更新进度条

if (dir.exists(working_dir)) {               # 判断工作目录是否存在
  setwd(working_dir)                         # 存在则设置工作目录
} else {
  stop("工作目录不存在，请检查路径！")         # 不存在则终止程序
}

if (!file.exists(inputFile)) {                # 判断数据文件是否存在
  stop("数据文件不存在，请检查路径：", inputFile)  # 不存在则终止程序
}

if (!file.exists(geneListFile)) {            # 判断基因列表文件是否存在
  stop("基因列表文件不存在，请检查路径：", geneListFile)  # 不存在则终止程序
}

message("读取基因表达数据...")
data_raw <- read.csv(inputFile, header = TRUE, check.names = FALSE, row.names = 1)  # 读取数据文件
message(sprintf("原始数据：%d 个基因，%d 个样本", nrow(data_raw), ncol(data_raw)))

message("读取基因列表...")
geneListData <- read.csv(geneListFile, header = TRUE, stringsAsFactors = FALSE)
geneList <- trimws(geneListData[, 1])  # 提取第一列作为基因名
geneList <- geneList[geneList != ""]
geneList <- unique(geneList)
message(sprintf("共读取 %d 个基因", length(geneList)))

# 筛选指定的基因
message("筛选指定的基因...")
availableGenes <- intersect(geneList, rownames(data_raw))
message(sprintf("在表达数据中找到 %d 个基因", length(availableGenes)))

if (length(availableGenes) == 0) {
  stop("指定的基因列表中没有在表达数据中找到的基因！请检查基因名称是否匹配。")
}

# 只保留指定的基因
data_raw <- data_raw[availableGenes, , drop = FALSE]
message(sprintf("筛选后数据：%d 个基因，%d 个样本", nrow(data_raw), ncol(data_raw)))

data_transposed <- t(data_raw)                # 转置数据：行代表样品，列代表基因

 ###########################################
 # 5. 构建 LASSO 回归模型
 ###########################################
 current_step <- current_step + 1            # 步骤 3：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 x_matrix <- as.matrix(data_transposed)        # 将转置数据转换为矩阵（基因表达矩阵）
 sample_groups <- gsub("(.*)_(.*)", "\\2", rownames(data_transposed))  # 从样品名称中提取分组信息
 lasso_model <- glmnet(x_matrix, sample_groups, family = "binomial", alpha = 1)  # 构建 LASSO 模型

 ###########################################
 # 6. 使用交叉验证选择最佳 lambda 值
 ###########################################
 current_step <- current_step + 1            # 步骤 4：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 cv_model <- cv.glmnet(x_matrix, sample_groups, family = "binomial", alpha = 1,
                       type.measure = "deviance", nfolds = 10)  # 交叉验证

 ###########################################
 # 7. 提取非零系数的特征基因，并输出为 CSV 文件
 ###########################################
 current_step <- current_step + 1            # 步骤 5：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 coef_matrix <- coef(lasso_model, s = cv_model$lambda.min)  # 获取lambda.min下的系数（更多基因）
 nonzero_idx <- which(coef_matrix != 0)        # 查找非零系数的索引
 gene_names <- rownames(coef_matrix)[nonzero_idx]  # 提取对应基因名称

 # 直接使用lambda.min对应的基因：与CV_curve_plot_A_gene_count.pdf中左边虚线对应
 cat("使用lambda.min对应的基因（与CV曲线左边虚线对应）...\n")

 # 获取lambda.min下的系数（完整提取）
 coef_lambda_min_full <- coef(lasso_model, s = cv_model$lambda.min)
 # 转为向量格式并去除截距
 coef_lambda_min <- as.numeric(coef_lambda_min_full[-1])
 gene_names_list <- rownames(coef_lambda_min_full)[-1]
 names(coef_lambda_min) <- gene_names_list

 # 找到非零系数的基因
 nonzero_mask <- coef_lambda_min != 0
 nonzero_genes_min <- gene_names_list[nonzero_mask]

 # 按照系数绝对值大小排序
 coef_nonzero <- coef_lambda_min[nonzero_mask]
 all_selected_genes_sorted <- nonzero_genes_min[order(abs(coef_nonzero), decreasing = TRUE)]

 genes_df <- data.frame(gene = all_selected_genes_sorted)
 write.csv(genes_df, file = "LASSO.gene.csv", quote = FALSE, row.names = FALSE)  # 保存为 CSV 文件
 cat("已保存", nrow(genes_df), "个基因到 LASSO.gene.csv\n")
 cat("这些基因对应CV_curve_plot_A_gene_count.pdf中左边虚线（lambda.min）的基因\n")

 # 提取LASSO基因的表达矩阵
 cat("\n正在提取LASSO基因的表达矩阵...\n")

 # 从原始数据中提取这些基因的表达数据
 lasso_gene_names <- genes_df$gene

 # 检查基因是否在data_raw中存在
 lasso_genes_in_data <- lasso_gene_names[lasso_gene_names %in% rownames(data_raw)]

 if (length(lasso_genes_in_data) > 0) {
   # 提取这些基因的表达矩阵
   lasso_gene_exp <- data_raw[lasso_genes_in_data, ]

   # 保存为CSV文件
   write.csv(lasso_gene_exp, file = "LASSO.gene.exp.csv", quote = FALSE)
   cat("已保存", nrow(lasso_gene_exp), "个基因的表达矩阵到 LASSO.gene.exp.csv\n")
   cat("矩阵维度: ", nrow(lasso_gene_exp), "个基因 x ", ncol(lasso_gene_exp), "个样本\n")
 } else {
   warning("警告: 无法在表达矩阵中找到选中的基因!")
 }


 ###########################################
 # 8. 整理交叉验证数据，准备绘图
 ###########################################
 current_step <- current_step + 1            # 步骤 6：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 cv_data <- data.frame(lambda = cv_model$lambda,
                       cvm = cv_model$cvm,
                       cvup = cv_model$cvup,
                       cvlo = cv_model$cvlo)  # 构造 CV 数据框

 # 计算每个lambda对应的基因数量
 gene_counts <- apply(lasso_model$beta, 2, function(x) sum(x != 0))  # 计算每个lambda下非零系数的基因数量
 cv_data$gene_count <- gene_counts[1:length(cv_model$lambda)]  # 添加基因数量列

 cv_data$log_lambda <- log(cv_data$lambda)     # 计算 lambda 的对数值

 ###########################################
 # 9. 绘制交叉验证（CV）曲线，采用高级图形风格
 ###########################################
 current_step <- current_step + 1            # 步骤 7：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 cv_plot <- ggplot(cv_data, aes(x = gene_count, y = cvm)) +  # 使用基因数量作为x轴
   geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.5, color = "#7F7F7F") +  # 添加误差条
   geom_point(size = 3, color = "#0072B2") +  # 添加数据点
   geom_vline(xintercept = length(coef(lasso_model, s = cv_model$lambda.min)[coef(lasso_model, s = cv_model$lambda.min) != 0]) - 1, linetype = "dashed", color = "red", size = 1) +  # 标记 lambda.min对应的基因数（红色虚线）
   geom_vline(xintercept = length(coef(lasso_model, s = cv_model$lambda.1se)[coef(lasso_model, s = cv_model$lambda.1se) != 0]) - 1, linetype = "dashed", color = "red", size = 1) +  # 标记 lambda.1se对应的基因数（红色虚线）
   labs(x = "Number of Selected Genes", y = "Mean Cross-Validated Deviance",
        title = "CV Curve for LASSO Model (Gene Count)") +  # 修改坐标轴标签和标题
   theme_bw(base_size = 16) +                 # 使用黑白主题
   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
         panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
         panel.grid.minor = element_line(color = "grey90", linetype = "dotted"),
         plot.title = element_text(face = "bold", hjust = 0.5, color = "#333333"))

 ###########################################
 # 10. 绘制 LASSO 系数路径图，采用新思路及高级美学风格
 ###########################################
 current_step <- current_step + 1            # 步骤 8：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 # 将 LASSO 模型的系数矩阵转换为长格式数据，便于 ggplot 绘图
 lasso_long <- melt(data.frame(log_lambda = log(lasso_model$lambda),
                               t(as.matrix(lasso_model$beta))), id.vars = "log_lambda")
 # 对每个变量，选取绝对值最大的点作为代表
 representative_points <- lapply(unique(lasso_long$variable), function(v) {
   sub_data <- subset(lasso_long, variable == v)  # 提取当前变量数据
   idx <- which.max(abs(sub_data$value))          # 找到绝对值最大的索引
   data.frame(log_lambda = sub_data$log_lambda[idx], variable = v, value = sub_data$value[idx])
 })
 representative_df <- do.call(rbind, representative_points)  # 合并代表性数据
 n_vars <- length(unique(lasso_long$variable))  # 计算变量总数
 color_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(n_vars)  # 生成配色方案
 # 绘制 LASSO 系数路径图
 lasso_plot <- ggplot(lasso_long, aes(x = log_lambda, y = value, group = variable, color = variable)) +
   geom_line(size = 1.2, alpha = 0.9) +         # 绘制路径线
   geom_point(data = representative_df, aes(x = log_lambda, y = value), size = 3) +  # 添加代表性数据点
   geom_text_repel(data = representative_df, aes(x = log_lambda, y = value, label = variable),
                   size = 4, fontface = "bold", box.padding = 0.35, point.padding = 0.5,
                   segment.color = "#555555", segment.size = 0.8) +  # 添加不重叠标签
   geom_vline(xintercept = log(cv_model$lambda.min), linetype = "dashed", size = 1, color = "#D55E00") +  # 标记 lambda.min
   geom_vline(xintercept = log(cv_model$lambda.1se), linetype = "dashed", size = 1, color = "#D55E00") +  # 标记 lambda.1se
   labs(x = "Log(lambda)", y = "Coefficient Value", title = "LASSO Coefficient Path") +  # 添加标签和标题
   scale_color_manual(values = color_palette) +  # 应用自定义颜色
   theme_bw(base_size = 16) +                    # 使用黑白主题
   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
         panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
         panel.grid.minor = element_line(color = "grey90", linetype = "dotted"),
         legend.position = "none",
         plot.title = element_text(face = "bold", hjust = 0.5, color = "#333333"))


 ###########################################
 # 11.1 单独输出图 A（CV 曲线）为 PDF
 ###########################################
 current_step <- current_step + 1            # 步骤 9.1：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 pdf("CV_curve_plot_A_gene_count.pdf", width = 6, height = 5.5)  # 打开 PDF 设备（基因数量版本）
 print(cv_plot)                               # 输出 CV 曲线图到 PDF 文件
 dev.off()                                    # 关闭 PDF 设备

 # 创建原版CV曲线图（使用Log(lambda)作为x轴）并保存
 cv_plot_original <- ggplot(cv_data, aes(x = log_lambda, y = cvm)) +  # 使用log_lambda作为x轴
   geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05, color = "#7F7F7F") +  # 添加误差条
   geom_point(size = 3, color = "#0072B2") +  # 添加数据点
   geom_vline(xintercept = log(cv_model$lambda.min), linetype = "dashed", color = "red", size = 1) +  # 标记 lambda.min（红色虚线）
   geom_vline(xintercept = log(cv_model$lambda.1se), linetype = "dashed", color = "red", size = 1) +  # 标记 lambda.1se（红色虚线）
   labs(x = "Log(lambda)", y = "Mean Cross-Validated Deviance",
        title = "CV Curve for LASSO Model") +  # 原版标题
   theme_bw(base_size = 16) +                 # 使用黑白主题
   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
         panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
         panel.grid.minor = element_line(color = "grey90", linetype = "dotted"),
         plot.title = element_text(face = "bold", hjust = 0.5, color = "#333333"))

 pdf("CV_curve_plot_A.pdf", width = 6, height = 5.5)  # 打开 PDF 设备（原版Log(lambda)）
 print(cv_plot_original)                      # 输出原版CV曲线图到 PDF 文件
 dev.off()                                    # 关闭 PDF 设备

 ###########################################
 # 11.2 单独输出图 B（LASSO 系数路径）为 PDF
 ###########################################
 current_step <- current_step + 1            # 步骤 9.2：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 pdf("LASSO_coefficient_path_plot_B.pdf", width = 7, height = 6)  # 打开 PDF 设备
 print(lasso_plot)                            # 输出 LASSO 系数路径图到 PDF 文件
 dev.off()                                    # 关闭 PDF 设备

 ###########################################
 # 12. 输出表格：交叉验证结果、非零系数详细表及完整系数矩阵（均为 CSV 格式）
 ###########################################
 current_step <- current_step + 1            # 步骤 12：更新步骤计数
 setTxtProgressBar(pb, current_step)          # 更新进度条
 # 输出交叉验证结果表为 CSV 文件
 write.csv(cv_data, file = "cv_results_table_enhanced.csv", quote = FALSE, row.names = FALSE)


 # 输出完整系数矩阵表（不含截距）为 CSV 文件
 all_coef_matrix <- as.matrix(lasso_model$beta)            # 提取所有系数矩阵
 all_coef_df <- as.data.frame(all_coef_matrix)             # 转换为数据框
 all_coef_df$gene <- rownames(all_coef_df)                  # 添加 gene 列
 all_coef_df <- all_coef_df[, c("gene", setdiff(names(all_coef_df), "gene"))]  # 调整列顺序
 write.csv(all_coef_df, file = "all_coefficients_table_enhanced.csv", quote = FALSE, row.names = FALSE)

 ###########################################
 # 15. 完成所有步骤，关闭进度条并输出完成信息
 ###########################################
 close(pb)                                  # 关闭进度条
 cat("\n所有步骤已完成！\n")                # 控制台输出完成提示
