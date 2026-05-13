# =============================================================================
# 编号       : R033
# 脚本名     : gwas为暴露，芬兰为结局。MR分析.R
# 分类       : 09_孟德尔随机化_GWAS处理
# 项目来源   : MR_GWAS暴露_芬兰结局
# 用途       : 将 GWAS 暴露与 FinnGen 结局数据进行孟德尔随机化分析并输出敏感性分析图表。
# 结果图     : 散点图；森林图；漏斗图；Leave-one-out敏感性图；曼哈顿图；QQ图；Radial MR图
# 主要 R 包  : gwasglue; qqman; RadialMR; TwoSampleMR; VariantAnnotation
# 整理时间   : 2026-05-10
# =============================================================================
# 加载所需的R包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
library(qqman)
library(RadialMR)

# 设置工作目录
setwd("H:\\常用分析生信\\35.孟德尔随机化分析，暴露因素GWAS，结局芬兰\\新建文件夹 (2)")

# 获取当前工作目录下所有的 .gz 文件名
file_list <- list.files(pattern = "\\.gz$")
total_files <- length(file_list)

# 检查是否有文件
if (total_files > 0) {
  for (i in seq_along(file_list)) {
    file <- file_list[i]
    
    # 显示当前处理的文件及进度
    message(sprintf("正在处理文件 %d/%d: %s", i, total_files, file))
    
    # 生成输出文件名前缀，例如将 "filename.h.tsv.gz" 转换为 "filename"
    output_prefix <- gsub("\\.gz$", "", file)
    
    # 读取暴露数据
    exposure_dat <- read_exposure_data(
      filename = "exposure.F.csv",  # 暴露数据文件名
      sep = ",",                          # 文件分隔符
      snp_col = "SNP",                    # SNP列名
      beta_col = "beta.exposure",         # beta值列名
      se_col = "se.exposure",             # 标准误列名
      pval_col = "pval.exposure",         # p值列名
      effect_allele_col = "effect_allele.exposure", # 效应等位基因列名
      other_allele_col = "other_allele.exposure",   # 其他等位基因列名
      eaf_col = "eaf.exposure",           # 效应等位基因频率列名
      phenotype_col = "exposure",         # 表型列名
      samplesize_col = "samplesize.exposure", # 样本量列名
      chr_col = "chr.exposure",           # 染色体列名
      pos_col = "pos.exposure",           # 位置列名
      clump = FALSE                       # 是否进行clumping
    )
    
    # 读取整理好的结局数据
    outcome_data <- read_outcome_data(
      snps = exposure_dat$SNP,          # 使用暴露数据中的SNP列表
      filename = file,                  # 当前结局数据文件
      sep = "\t",                       # 指定分隔符
      snp_col = "rsids",                # 指定SNP列名
      beta_col = "beta",                # 指定beta值列名
      se_col = "sebeta",                # 指定标准误列名
      effect_allele_col = "alt",        # 指定效应等位基因列名
      other_allele_col = "ref",         # 指定其他等位基因列名
      pval_col = "pval",                # 指定p值列名
      eaf_col = "af_alt"                # 效应等位基因频率列名
    )
    
    # 将暴露数据和结局数据合并
    outcome_data$outcome <- "lumbar disc prolapse"  # 设置结局名称
    dat <- harmonise_data(exposure_dat, outcome_data)  # 合并数据
    
    # 输出用于孟德尔随机化的工具变量
    outTab <- dat[dat$mr_keep == "TRUE", ]
    write.csv(outTab, file = paste0(output_prefix, "_table.SNP.csv"), row.names = FALSE)
    
    # 孟德尔随机化分析
    mrResult <- mr(dat)
    
    # 对结果进行OR值的计算
    mrTab <- generate_odds_ratios(mrResult)
    write.csv(mrTab, file = paste0(output_prefix, "_table.MRresult.csv"), row.names = FALSE)
    
    # 异质性分析
    heterTab <- mr_heterogeneity(dat)
    write.csv(heterTab, file = paste0(output_prefix, "_table.heterogeneity.csv"), row.names = FALSE)
    
    # 多效性检验
    pleioTab <- mr_pleiotropy_test(dat)
    write.csv(pleioTab, file = paste0(output_prefix, "_table.pleiotropy.csv"), row.names = FALSE)
    
    # 绘制散点图并保存为PDF文件
    pdf(file = paste0(output_prefix, "_scatter_plot.pdf"), width = 7, height = 6.5)
    p1 <- mr_scatter_plot(mrResult, dat)
    print(p1)
    dev.off()
    
    # 绘制森林图并保存为PDF文件
    res_single <- mr_singlesnp(dat)  # 得到每个工具变量对结局的影响
    pdf(file = paste0(output_prefix, "_forest.pdf"), width = 6.5, height = 5)
    p2 <- mr_forest_plot(res_single)
    print(p2)
    dev.off()
    
    # 绘制漏斗图并保存为PDF文件
    pdf(file = paste0(output_prefix, "_funnel_plot.pdf"), width = 6.5, height = 6)
    p3 <- mr_funnel_plot(singlesnp_results = res_single)
    print(p3)
    dev.off()
    
    # 留一法敏感性分析并保存为PDF文件
    pdf(file = paste0(output_prefix, "_leaveoneout.pdf"), width = 6.5, height = 5)
    p4 <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    print(p4)
    dev.off()
    
    ### GWAS全局结果展示（曼哈顿图和QQ图）
    # 假设暴露数据文件 exposure.F.csv 包含列：chr.exposure, pos.exposure, pval.exposure
    exposure_gwas <- read.csv("exposure.F.csv", stringsAsFactors = FALSE)
    
    # 绘制曼哈顿图
    pdf(file = paste0(output_prefix, "_manhattan_plot.pdf"), width = 8, height = 6)
    manhattan(exposure_gwas, 
              chr = "chr.exposure", 
              bp = "pos.exposure", 
              p = "pval.exposure", 
              main = "Manhattan Plot for Exposure GWAS")
    dev.off()
    
    # 绘制QQ图
    pdf(file = paste0(output_prefix, "_qq_plot.pdf"), width = 8, height = 6)
    qq(exposure_gwas$pval.exposure, main = "QQ Plot for Exposure GWAS")
    dev.off()
    
    ### 扩展的MR方法分析
    # (1) MR Egger、加权中位数、加权模式
    mr_extended <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    write.csv(mr_extended, file = paste0(output_prefix, "_table.MR_extended.csv"), row.names = FALSE)
    
    # (2) Steiger方向性检验
    steiger <- steiger_filtering(dat)
    write.csv(steiger, file = paste0(output_prefix, "_table.MR_steiger.csv"), row.names = FALSE)
    ### Radial MR分析
    
    # 1. 格式化数据（注意直接按顺序传入参数，避免使用未定义的参数名）
    dat_radial <- format_radial(
      dat$beta.exposure,
      dat$beta.outcome,
      dat$se.exposure,
      dat$se.outcome,
      dat$SNP
    )
    
    # 2. 进行 IVW Radial MR 回归分析
    radial_ivw <- ivw_radial(dat_radial, alpha = 0.05, weights = 3)
    
    # 3. 绘制 Radial Plot 并保存为 PDF 文件
    # 请根据实际情况设置输出前缀，例如使用 "radial_analysis" 作为前缀
    output_prefix <- "radial_analysis"
    pdf(file = paste0(output_prefix, "_radial_plot.pdf"), width = 7, height = 6.5)
    print(plot_radial(radial_ivw))
    dev.off()
    
    outliers <- radial_ivw$outliers  # 或根据实际情况进行提取
    write.csv(outliers, file = "radial_ivw_outliers.csv", row.names = FALSE)
    
    variant_data <- radial_ivw$data
    head(variant_data)
    outlier_variants <- radial_ivw$outliers
    head(outlier_variants)
    
    write.csv(radial_ivw$data, file = "radial_ivw_variant_data.csv", row.names = FALSE)
    
    # 显示当前文件处理完成
    message(sprintf("文件处理完成: %s (%d/%d)", file, i, total_files))
  }
} else {
  message("目录中未找到 .gz 文件。")
}





message("所有文件处理完成。")