# =============================================================================
# 编号       : R028
# 脚本名     : 01gwas数据库，VCF格式.筛选强关联的SNP.R
# 分类       : 09_孟德尔随机化_GWAS处理
# 项目来源   : MR_GWAS暴露数据处理
# 用途       : 处理 GWAS/VCF 暴露数据，筛选显著关联 SNP。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : gwasglue; pbapply; TwoSampleMR; VariantAnnotation
# 整理时间   : 2026-05-10
# =============================================================================
# 安装并加载进度条相关的包
# 如果尚未安装 pbapply 包，可以取消下面的注释并运行
# install.packages("pbapply")
library(VariantAnnotation)  # 用于处理VCF文件的包
library(gwasglue)           # 用于GWAS数据转换的包
library(TwoSampleMR)        # 用于双样本孟德尔随机化分析的包
library(pbapply)            # 用于添加进度条的包

# 设置工作目录到数据文件所在位置，确保路径正确
setwd("H:\\常用分析生信\\34，孟德尔随机化分析处理GWAS数据\\02gwas.mrcieu数据筛选，筛选强关联的SNP") 

# 定义暴露数据的文件名路径
inputFilePath <- "ukb-b-4522.vcf.gz"  # 曝露数据文件路径

# 判断文件是否存在，如果不存在则给出提示
if (!file.exists(inputFilePath)) {
  stop("文件不存在，请检查文件路径是否正确。")
} else {
  message("文件加载成功，开始读取数据...")
}

# 读取VCF文件格式的暴露数据，并进行数据转换
vcfData <- readVcf(inputFilePath)

# 使用gwasvcf_to_TwoSampleMR函数将VCF数据转化为TwoSampleMR包可以使用的数据格式
transformedExposureData <- gwasvcf_to_TwoSampleMR(vcf = vcfData, type = "exposure")

# 将转换后的数据保存为CSV文件，便于后续分析
write.csv(transformedExposureData, file = "transformedExposureData.csv")
message("暴露数据已保存为 'transformedExposureData.csv'")

# 使用进度条进行过滤操作
# 我们通过判断数据中每一行的p值来筛选数据
message("开始对暴露数据进行过滤...")

# 设置进度条的参数
totalRows <- nrow(transformedExposureData)  # 总行数

# 使用 pbapply 包提供的 pbsapply 函数，添加进度条
filteredExposureData <- pbsapply(1:totalRows, function(i) {
  if (transformedExposureData$pval.exposure[i] < 5e-06) {
    return(transformedExposureData[i, ])
  } else {
    return(NULL)
  }
}, simplify = FALSE)

# 将结果合并回数据框
filteredExposureData <- do.call(rbind, filteredExposureData)

# 如果过滤后的数据为空，提醒用户
if (nrow(filteredExposureData) == 0) {
  stop("没有数据符合筛选条件，检查过滤条件是否正确。")
} else {
  message("数据过滤完成，共保留了 ", nrow(filteredExposureData), " 条记录。")
}

# 保存过滤后的数据为新的CSV文件
write.csv(filteredExposureData, file = "filteredExposureData.csv", row.names = FALSE)
message("过滤后的数据已保存为 'filteredExposureData.csv'")
