# =============================================================================
# 编号       : R031
# 脚本名     : 去除弱工具变量.R
# 分类       : 09_孟德尔随机化_GWAS处理
# 项目来源   : MR_GWAS暴露数据处理
# 用途       : 计算 F 统计量并剔除弱工具变量。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : pkg
# 整理时间   : 2026-05-10
# =============================================================================
# ------------------------ 第0步：安装和加载必要的包 ------------------------
# 这部分代码确保我们拥有所有必要的库，如果没有安装，将自动安装
# 如果你已经安装了这些库，可以跳过这一步

# 定义需要的包名
required_packages <- c("ieugwasr", "httr", "jsonlite")

# 检查每个包是否安装，未安装则安装并加载
for(pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)  # 安装缺少的包
    library(pkg, character.only = TRUE)  # 加载已安装的包
  } else {
    library(pkg, character.only = TRUE)  # 如果已安装，直接加载
  }
}

redundant_info <- "环境配置完毕" 

# ------------------------ 第1步：设置工作目录并读取数据 ------------------------
# 设置工作目录
working_directory <- "H:\\常用分析生信\\34，孟德尔随机化分析处理GWAS数据\\05.去除弱工具变量"
if (dir.exists(working_directory)) {
  setwd(working_directory)  # 设置工作目录
} else {
  stop("错误：工作目录不存在！")  # 如果目录不存在，停止程序
}

# 设置输入文件名
inputFile <- "filtered.eaf.csv"
if (!file.exists(inputFile)) {
  stop("错误：输入文件不存在，程序终止！")  # 如果文件不存在，停止程序
}

# 读取 CSV 文件
data_input <- read.csv(inputFile, header = TRUE, sep = ",", check.names = FALSE)
if (nrow(data_input) == 0) {
  stop("错误：输入数据为空，无法继续处理！")  # 如果数据为空，停止程序
}

print(head(data_input))  # 打印数据前几行供检查

# ------------------------ 第2步：计算 R² 和 F 值 ------------------------
# 初始化进度条：计算R2和F值，共有两个步骤
total_steps <- 2  # 设置总步骤数
progress_bar <- txtProgressBar(min = 0, max = total_steps, style = 3)  # 创建进度条

# 步骤1：计算 R²
data_input$R2 <- (2 * data_input$beta.exposure * data_input$beta.exposure * data_input$eaf.exposure * (1 - data_input$eaf.exposure)) / 
  (2 * data_input$beta.exposure * data_input$beta.exposure * data_input$eaf.exposure * (1 - data_input$eaf.exposure) + 
     2 * data_input$se.exposure * data_input$se.exposure * data_input$samplesize.exposure * data_input$eaf.exposure * (1 - data_input$eaf.exposure))
# R²（决定系数）是衡量工具变量解释因变量变异的比例。这个公式是在利用工具变量的β值、标准误、样本量和暴露等位点频率（EAF）来计算R²值。公式的分母是R²的计算公式中的方程，分子是与β、样本大小、EAF相关的部分。

setTxtProgressBar(progress_bar, 1)  # 更新进度条至步骤1完成
message("步骤1：R² 值计算完成。")  # 输出步骤1的完成信息

# 步骤2：计算 F 值
data_input$F <- data_input$R2 * (data_input$samplesize.exposure - 2) / (1 - data_input$R2)
# F值计算是基于R²和样本量（samplesize.exposure）。F值用于检验工具变量的有效性。如果F值大于10，通常表示工具变量足够强，有效用于孟德尔随机化分析。

setTxtProgressBar(progress_bar, 2)  # 更新进度条至步骤2完成
message("步骤2：F 值计算完成。")  # 输出步骤2的完成信息

# ------------------------ 第3步：筛选 F 值大于 10 的数据 ------------------------
# 步骤3：筛选 F 值大于10的数据，并输出到新文件
filtered_data <- data_input[as.numeric(data_input$F) > 10, ]  # 筛选 F 值大于 10 的数据
if (nrow(filtered_data) == 0) {
  stop("错误：没有 F 值大于 10 的数据！")  # 如果筛选后的数据为空，停止程序
}

# 输出结果到 CSV 文件
output_file <- "exposure.F.csv"
write.csv(filtered_data, file = output_file, row.names = FALSE)  # 将结果写入文件

# 判断文件是否成功保存
if (file.exists(output_file)) {
  message("步骤3：筛选结果已成功保存至文件：", output_file)  # 输出成功提示
} else {
  stop("错误：输出文件未成功生成！")  # 如果输出文件不存在，停止程序
}

# ------------------------ 第4步：关闭进度条并输出结束信息 ------------------------
close(progress_bar)  # 关闭进度条
message("所有步骤已成功完成！")  # 输出完成信息

dummy_var1 <- data_input$F  # 获取 F 值的副本
dummy_var2 <- mean(dummy_var1, na.rm = TRUE)  # 计算 F 值的平均值
message("操作完成：F 值的平均值为 ", dummy_var2)  

# ------------------------ 代码结束 ------------------------
