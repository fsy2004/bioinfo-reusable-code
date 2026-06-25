# =============================================================================
# 编号       : R029
# 脚本名     : 02去除连锁不平衡.R
# 分类       : 09_mendelian_randomization
# 项目来源   : MR_GWAS暴露数据处理
# 用途       : 对候选 SNP 去除连锁不平衡，保留独立工具变量。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : pkg
# 整理时间   : 2026-05-10
# =============================================================================
# 设置工作目录，确保后续文件操作在指定目录下进行  # 设置工作目录
setwd("H:\\常用分析生信\\34，孟德尔随机化分析处理GWAS数据\\03去除连锁不平衡性")  # 指定工作目录

# 定义所需使用的包名称，并检查是否已经安装  # 定义必需的包
required_packages <- c("VariantAnnotation", "gwasglue", "TwoSampleMR")  # 列出所需包

# 遍历每个包，若未安装则进行安装，再加载包  # 检查和加载包
for(pkg in required_packages) {  # 循环处理每个包
  if (!require(pkg, character.only = TRUE)) {  # 如果包未加载成功
    install.packages(pkg)  # 安装该包
    library(pkg, character.only = TRUE)  # 加载包
  } else {  # 如果包已存在
    library(pkg, character.only = TRUE)  # 直接加载包
  }
}

# 初始化进度条与步骤 
total_steps <- 4  
progress <- txtProgressBar(min = 0, max = total_steps, style = 3)  

# =================== 第一步：读取暴露数据 ===================  #
step <- 1  
setTxtProgressBar(progress, step)  # 更新进度条显示
input_file <- "filteredExposureData.csv"  # 定义输入文件名称

# 判断输入文件是否存在，若不存在则终止程序  # 检查输入文件是否存在
if (!file.exists(input_file)) {  # 如果文件不存在
  stop("输入文件不存在，程序终止！")  # 抛出错误信息并停止执行
}

# 调用 read_exposure_data 函数读取CSV文件中的暴露数据  # 读取暴露数据
exposure_data <- read_exposure_data(
  filename = input_file,               # 指定输入文件
  sep = ",",                           # 指定分隔符为逗号
  snp_col = "SNP",                     # 指定SNP所在列
  beta_col = "beta.exposure",          # 指定效应值所在列
  se_col = "se.exposure",              # 指定标准误所在列
  pval_col = "pval.exposure",          # 指定p值所在列
  effect_allele_col = "effect_allele.exposure",  # 指定效应等位基因列
  other_allele_col = "other_allele.exposure",      # 指定非效应等位基因列
  samplesize_col = "samplesize.exposure",          # 指定样本量所在列
  chr_col = "chr.exposure",            # 指定染色体列
  pos_col = "pos.exposure",            # 指定位点列
  clump = FALSE                        # 关闭默认的clumping操作
)  # 执行数据读取操作

# =================== 第二步：数据有效性检查 ===================  # 步骤2说明
step <- 2  # 更新当前步骤为2
setTxtProgressBar(progress, step)  # 更新进度条显示

# 判断读取的数据是否为空  # 检查数据是否为空
if (nrow(exposure_data) == 0) {  # 如果数据行数为0
  stop("读取的暴露数据为空，无法继续处理！")  # 抛出错误信息
}

# 冗余操作：对数据做摘要统计并打印（便于调试）  # 冗余检查步骤
data_summary <- summary(exposure_data)  # 获取数据摘要
print(data_summary)  # 打印数据摘要信息

# =================== 第三步：剔除连锁不平衡的SNP（Clumping） ===================  # 步骤3说明
step <- 3  # 当前步骤更新为3
setTxtProgressBar(progress, step)  # 更新进度条显示

# 使用 clump_data 函数对暴露数据进行 clumping 操作，移除连锁不平衡的 SNP  # 执行clumping操作
processed_data <- clump_data(
  exposure_data,    # 输入原始暴露数据
  clump_kb = 10000, # 设置clumping窗口大小为 10000 kb   最小500
  clump_r2 = 0.001  # 设置clumping相关系数阈值为 0.001  最大0.01
)  # 返回经过处理的暴露数据

# 添加冗余判断：检查处理后的数据是否发生变化（仅作示例）  # 冗余逻辑检查
if (identical(exposure_data, processed_data)) {  # 如果处理前后数据完全一致
  message("警告：clumping操作未能改变数据！")  # 打印警告信息
} else {  # 否则
  message("clumping操作已成功执行。")  # 打印确认信息
}

# =================== 第四步：保存处理结果 ===================  # 步骤4说明
step <- 4  # 当前步骤更新为4
setTxtProgressBar(progress, step)  # 更新进度条显示

output_file <- "filtered_clumped.csv"  # 定义输出文件名称

# 判断处理后的数据是否非空，若非空则写入CSV文件，否则给出警告  # 判断数据是否为空再保存
if (nrow(processed_data) > 0) {  # 如果处理数据的行数大于0
  write.csv(processed_data, file = output_file, row.names = FALSE)  # 保存处理数据到CSV文件，不写入行名
  message("输出文件已成功生成：", output_file)  # 打印生成成功信息
} else {  # 如果处理数据为空
  warning("处理后的数据为空，输出文件未生成！")  # 输出警告信息
}

# 关闭进度条  # 结束进度条
close(progress)  # 关闭进度条显示

# 附加判断：确认输出文件是否存在，并给出最终提示  # 最终输出文件检查
if (file.exists(output_file)) {  # 如果输出文件存在
  message("数据处理流程完成，输出文件生成成功！")  # 打印成功完成信息
} else {  # 如果输出文件不存在
  message("数据处理流程完成，但未检测到输出文件，请检查代码流程！")  # 打印错误提示信息
}
