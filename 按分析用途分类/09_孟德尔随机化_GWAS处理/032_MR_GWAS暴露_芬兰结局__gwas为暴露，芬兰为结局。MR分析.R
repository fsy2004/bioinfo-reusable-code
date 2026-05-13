# =============================================================================
# 编号       : R032
# 脚本名     : gwas为暴露，芬兰为结局。MR分析.R
# 分类       : 09_孟德尔随机化_GWAS处理
# 项目来源   : MR_GWAS暴露_芬兰结局
# 用途       : 将 GWAS 暴露与 FinnGen 结局数据进行孟德尔随机化分析并输出敏感性分析图表。
# 结果图     : 散点图；森林图；漏斗图；Leave-one-out敏感性图；曼哈顿图；QQ图；Radial MR图
# 主要 R 包  : gwasglue; qqman; RadialMR; TwoSampleMR; VariantAnnotation
# 整理时间   : 2026-05-10
# =============================================================================
# ------------------------------ #
# 1. 加载所需的R包并设置工作目录
# ------------------------------ #
library(VariantAnnotation)  # 加载用于变异注释的包
library(gwasglue)           # 加载用于GWAS数据整合的包
library(TwoSampleMR)        # 加载用于孟德尔随机化分析的包
library(qqman)              # 加载用于生成曼哈顿图和QQ图的包
library(RadialMR)           # 加载用于Radial MR分析的包

setwd("H:\\常用分析生信\\35.孟德尔随机化分析，暴露因素GWAS，结局芬兰")  # 设置工作目录

# ------------------------------ #
# 2. 获取所有 .gz 文件并初始化进度条
# ------------------------------ #
gz_files <- list.files(pattern = "\\.gz$")  # 获取当前工作目录下所有以 .gz 结尾的文件
num_files <- length(gz_files)                # 计算文件总数

if (num_files > 0) {  # 判断是否存在文件
  pb <- txtProgressBar(min = 0, max = num_files, style = 3)  # 初始化文本进度条
  
  # ------------------------------ #
  # 3. 循环处理每个 .gz 文件
  # ------------------------------ #
  for (file_idx in seq_along(gz_files)) {
    current_file <- gz_files[file_idx]  # 当前处理的文件名
    
    # 输出处理进度和当前步骤信息
    message(sprintf("步骤1：正在处理文件 %d/%d: %s", file_idx, num_files, current_file))
    
    # 生成输出文件名前缀：将 ".gz" 去除
    out_prefix <- sub("\\.gz$", "", current_file)
    
    # ------------------------------ #
    # 4. 读取并准备数据
    # ------------------------------ #
    # 4.1 读取暴露数据（exposure data）
    exposure_data <- read_exposure_data(
      filename = "exposure.F.csv",                # 暴露数据文件名
      sep = ",",                                  # CSV分隔符
      snp_col = "SNP",                            # SNP列名
      beta_col = "beta.exposure",                 # beta值列名
      se_col = "se.exposure",                     # 标准误列名
      pval_col = "pval.exposure",                 # p值列名
      effect_allele_col = "effect_allele.exposure",# 效应等位基因列名
      other_allele_col = "other_allele.exposure",  # 其他等位基因列名
      eaf_col = "eaf.exposure",                   # 效应等位基因频率列名
      phenotype_col = "exposure",                 # 表型列名
      samplesize_col = "samplesize.exposure",     # 样本量列名
      chr_col = "chr.exposure",                   # 染色体列名
      pos_col = "pos.exposure",                   # 位置信息列名
      clump = FALSE                               # 是否进行clumping（不进行）
    )
    # 判断暴露数据是否正确读取
    if (is.null(exposure_data) || nrow(exposure_data) == 0) {
      message(sprintf("警告：暴露数据为空，文件 %s 跳过。", current_file))
      next  # 跳过当前文件继续下一个循环
    }
    
    # 4.2 读取结局数据（outcome data），当前文件作为结局数据
    outcome_data <- read_outcome_data(
      snps = exposure_data$SNP,   # 使用暴露数据中的SNP列表
      filename = current_file,    # 当前结局数据文件
      sep = "\t",                 # 指定分隔符为制表符
      snp_col = "rsids",          # 结局数据中SNP列名
      beta_col = "beta",          # beta值列名
      se_col = "sebeta",          # 标准误列名
      effect_allele_col = "alt",  # 效应等位基因列名
      other_allele_col = "ref",   # 其他等位基因列名
      pval_col = "pval",          # p值列名
      eaf_col = "af_alt"          # 效应等位基因频率列名
    )
    # 判断结局数据是否正确读取
    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      message(sprintf("警告：结局数据为空，文件 %s 跳过。", current_file))
      next  # 跳过当前文件继续下一个循环
    }
    
    # 4.3 为结局数据添加结果名称
    outcome_data$outcome <- "BREAST cancer"  # 设置结局名称
    
    # 4.4 合并暴露和结局数据
    combined_data <- harmonise_data(exposure_data, outcome_data)
    
    # ------------------------------ #
    # 5. 孟德尔随机化（MR）分析及结果保存
    # ------------------------------ #
    # 5.1 输出用于MR分析的工具变量表
    mr_instruments <- combined_data[combined_data$mr_keep == "TRUE", ]
    write.csv(mr_instruments, file = paste0(out_prefix, "_SNP_table.csv"), row.names = FALSE)
    
    # 5.2 进行MR分析
    mr_results <- mr(combined_data)
    
    # 5.3 计算并保存OR值
    or_results <- generate_odds_ratios(mr_results)
    write.csv(or_results, file = paste0(out_prefix, "_MR_results.csv"), row.names = FALSE)
    
    # 5.4 进行异质性分析
    heterogeneity_results <- mr_heterogeneity(combined_data)
    write.csv(heterogeneity_results, file = paste0(out_prefix, "_heterogeneity.csv"), row.names = FALSE)
    
    # 5.5 进行多效性检验
    pleiotropy_results <- mr_pleiotropy_test(combined_data)
    write.csv(pleiotropy_results, file = paste0(out_prefix, "_pleiotropy.csv"), row.names = FALSE)
    
    # ------------------------------ #
    # 6. 生成各类图形并保存为PDF文件
    # ------------------------------ #
    # 6.1 绘制散点图（Scatter Plot）
    pdf(file = paste0(out_prefix, "_scatter_plot.pdf"), width = 7, height = 6.5)
    scatter_plot <- mr_scatter_plot(mr_results, combined_data)  # 生成散点图对象
    if (!is.null(scatter_plot)) {
      print(scatter_plot)  # 打印散点图到PDF
    } else {
      message("散点图生成失败或返回NULL。")
    }
    dev.off()  # 关闭PDF设备
    
    # 6.2 绘制森林图（Forest Plot）
    single_snp_results <- mr_singlesnp(combined_data)  # 获取单个SNP的MR结果
    pdf(file = paste0(out_prefix, "_forest_plot.pdf"), width = 6.5, height = 5)
    forest_plot <- mr_forest_plot(single_snp_results)  # 生成森林图对象
    if (!is.null(forest_plot)) {
      print(forest_plot)  # 打印森林图到PDF
    } else {
      message("森林图生成失败或返回NULL。")
    }
    dev.off()  # 关闭PDF设备
    
    # 6.3 绘制漏斗图（Funnel Plot）
    pdf(file = paste0(out_prefix, "_funnel_plot.pdf"), width = 6.5, height = 6)
    funnel_plot <- mr_funnel_plot(singlesnp_results = single_snp_results)  # 生成漏斗图对象
    if (!is.null(funnel_plot)) {
      print(funnel_plot)  # 打印漏斗图到PDF
    } else {
      message("漏斗图生成失败或返回NULL。")
    }
    dev.off()  # 关闭PDF设备
    
    # 6.4 进行留一法敏感性分析（Leave-One-Out）并绘图
    pdf(file = paste0(out_prefix, "_leaveoneout_plot.pdf"), width = 6.5, height = 5)
    leaveoneout_results <- mr_leaveoneout(combined_data)  # 计算留一法结果
    leaveoneout_plot <- mr_leaveoneout_plot(leaveoneout_results = leaveoneout_results)  # 生成留一法图对象
    if (!is.null(leaveoneout_plot)) {
      print(leaveoneout_plot)  # 打印留一法图到PDF
    } else {
      message("留一法敏感性图生成失败或返回NULL。")
    }
    dev.off()  # 关闭PDF设备
    
    # 6.5 绘制GWAS全局结果图：曼哈顿图和QQ图
    gwas_data <- read.csv("exposure.F.csv", stringsAsFactors = FALSE)  # 读取用于GWAS全局展示的暴露数据
    
    # 曼哈顿图
    pdf(file = paste0(out_prefix, "_manhattan_plot.pdf"), width = 8, height = 6)
    manhattan(gwas_data,
              chr = "chr.exposure",       # 指定染色体列
              bp = "pos.exposure",        # 指定位点信息列
              p = "pval.exposure",        # 指定p值列
              main = "Exposure GWAS Manhattan Plot")  # 设置图标题
    dev.off()  # 关闭PDF设备
    
    # QQ图
    pdf(file = paste0(out_prefix, "_qq_plot.pdf"), width = 8, height = 6)
    qq(gwas_data$pval.exposure, main = "Exposure GWAS QQ Plot")  # 绘制QQ图
    dev.off()  # 关闭PDF设备
    
    # ------------------------------ #
    # 7. 扩展MR方法与方向性检验
    # ------------------------------ #
    # 7.1 使用MR Egger、加权中位数和加权模式进行扩展MR分析
    extended_mr <- mr(combined_data, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    write.csv(extended_mr, file = paste0(out_prefix, "_MR_extended.csv"), row.names = FALSE)
    
    # 7.2 进行Steiger方向性检验
    steiger_test <- steiger_filtering(combined_data)
    write.csv(steiger_test, file = paste0(out_prefix, "_MR_steiger.csv"), row.names = FALSE)
    
    # ------------------------------ #
    # 8. Radial MR 分析及图形生成
    # ------------------------------ #
    # 8.1 格式化数据用于Radial分析
    radial_data <- format_radial(
      combined_data$beta.exposure,   # 暴露beta值
      combined_data$beta.outcome,    # 结局beta值
      combined_data$se.exposure,       # 暴露标准误
      combined_data$se.outcome,        # 结局标准误
      combined_data$SNP               # SNP标识
    )
    
    # 8.2 进行IVW Radial MR回归
    radial_results <- ivw_radial(radial_data, alpha = 0.05, weights = 3)
    
    # 8.3 生成Radial图并保存
    radial_prefix <- "radial_analysis"  # 定义Radial图输出前缀
    pdf(file = paste0(radial_prefix, "_radial_plot.pdf"), width = 7, height = 6.5)
    radial_plot_obj <- plot_radial(radial_results)  # 生成Radial图对象
    if (!is.null(radial_plot_obj)) {
      print(radial_plot_obj)  # 打印Radial图到PDF（适用于ggplot对象）
    } else {
      message("Radial图生成失败或返回NULL。")
    }
    dev.off()  # 关闭PDF设备
    
    # 8.4 保存Radial分析中识别的异常值（outliers）
    radial_outliers <- radial_results$outliers
    if (!is.null(radial_outliers) && nrow(radial_outliers) > 0) {
      write.csv(radial_outliers, file = "radial_outliers.csv", row.names = FALSE)
    } else {
      message("Radial分析中未检测到异常值或异常值数据为空。")
    }
    
    # 8.5 保存Radial分析中的SNP变异数据信息
    variant_info <- radial_results$data
    if (!is.null(variant_info) && nrow(variant_info) > 0) {
      write.csv(variant_info, file = "radial_variant_data.csv", row.names = FALSE)
    } else {
      message("Radial分析中未检测到变异数据或数据为空。")
    }
    
    # ------------------------------ #
    # 9. 更新进度并输出当前文件处理完成信息
    # ------------------------------ #
    setTxtProgressBar(pb, file_idx)  # 更新进度条
    message(sprintf("步骤完成：已处理文件 %s (%d/%d)", current_file, file_idx, num_files))
  }
  close(pb)  # 处理完所有文件后关闭进度条
} else {
  # 如果目录中没有找到 .gz 文件，则输出提示信息
  message("目录中未找到 .gz 文件，请检查工作目录。")
}

# ------------------------------ #
# 10. 全部处理完成后的提示信息
# ------------------------------ #
message("所有文件均已成功处理。")
