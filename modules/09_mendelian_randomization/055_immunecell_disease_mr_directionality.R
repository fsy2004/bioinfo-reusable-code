# ==========================================================================
# 脚本名     : 免疫细胞疾病MR方向性检验.R
# 分类       : 09_mendelian_randomization
# 项目来源   : 从压缩包 345.731种免疫细胞和疾病的MR分析，方向性检验.rar 整理
# 原始文件   : 345.731种免疫细胞和疾病的MR分析，方向性检验\免疫细胞和疾病进行MR分析 - GWAS数据库.R
# 用途       : 批量读取免疫细胞暴露与疾病结局GWAS/VCF数据，进行 TwoSampleMR 主分析、异质性、多效性和 Steiger 方向性检验。
# 结果图     : 本脚本主要输出结果表；可接入现有MR森林图/漏斗图/leave-one-out/Radial MR可视化模块
# 非肿瘤消化适配: 适合。可用于免疫细胞与非肿瘤消化系统疾病的因果方向筛选。
# 主要 R 包  : VariantAnnotation; gwasglue; TwoSampleMR; R.utils; ggplot2; RadialMR; dplyr
# 整理日期   : 2026-05-13
# 备注       : 保留bioinfo-reusable-code逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
library(VariantAnnotation)   # 用于处理变异注释数据
library(gwasglue)            # 用于 GWAS 数据转换
library(TwoSampleMR)         # 用于 Two-Sample Mendelian Randomization 分析
library(R.utils)             # 提供工具函数
library(ggplot2)             # 用于高级数据可视化
library(RadialMR)            # 加载用于Radial MR分析的包
library(dplyr)

# =============================================================================
# 参数设置
# =============================================================================

# 设置工作目录
work_dir <- "H:\\常用分析生信\\345.731种免疫细胞和疾病的MR分析，方向性检验"
setwd(work_dir)

# 定义输入文件路径
exposure_file <- "immune_exposure.csv"                    # 暴露数据文件
outcome_vcf_files <- list.files(pattern = "\\.vcf.gz$")   # 结局VCF文件列表
outcome_file <- outcome_vcf_files[1]                      # 当前分析的结局文件

# 从结局文件名提取疾病名称
disease_name <- tools::file_path_sans_ext(outcome_file)

# 输出目录
output_dir <- "results_output1"

# Steiger方向性检验所需的结局样本量
outcome_n <- 423258

# =============================================================================
# 数据读取与预处理
# =============================================================================

# 读取原始暴露数据，建立id.exposure到exposure名称的映射关系
exposure_raw <- read.csv(exposure_file, stringsAsFactors = FALSE)
exposure_id_mapping <- unique(exposure_raw[, c("id.exposure", "exposure")])

# 读取暴露数据（TwoSampleMR格式）
exposure_dat <- read_exposure_data(
  filename = exposure_file,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  pval_col = "pval.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  phenotype_col = "id.exposure",
  id_col = "id",
  samplesize_col = "samplesize.exposure",
  chr_col = "chr.exposure",
  pos_col = "pos.exposure",
  clump = FALSE
)

# 读取结局VCF文件并转换为TwoSampleMR格式
outcome_vcf <- readVcf(outcome_vcf_files)
outcome_dat_raw <- gwasvcf_to_TwoSampleMR(vcf = outcome_vcf, type = "outcome")

# 提取与暴露SNP匹配的结局数据
outcome_merged <- merge(exposure_dat, outcome_dat_raw, by.x = "SNP", by.y = "SNP")
write.csv(outcome_merged[, -(2:ncol(exposure_dat))],
          file = "outcome_instruments.csv",
          row.names = FALSE)

# 读取结局数据（TwoSampleMR格式）
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "outcome_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "pval.outcome",
  eaf_col = "eaf.outcome"
)

# 添加结局样本量（Steiger方向性检验必需）
outcome_dat$samplesize.outcome <- outcome_n

# =============================================================================
# MR分析循环
# =============================================================================

# 获取所有唯一暴露ID
unique_exposures <- unique(exposure_dat$exposure)
n_exposures <- length(unique_exposures)

# 初始化进度条
progress_bar <- txtProgressBar(min = 0, max = n_exposures, style = 3)
step_count <- 1

# 定义文件名清理函数
clean_filename <- function(name) {
  name <- gsub("[^[:alnum:] ]", " ", name)   # 替换非字母数字字符为空格
  name <- gsub("\\s+", " ", name)            # 合并多个空格
  name <- trimws(name)                        # 去除首尾空格
  return(name)
}

# 初始化结果存储数据框
results_or <- data.frame()           # OR结果
results_heterogeneity <- data.frame() # 异质性检验结果
results_pleiotropy <- data.frame()    # 多效性检验结果
results_steiger <- data.frame()       # Steiger方向性检验结果

# 主循环：遍历每个暴露进行MR分析
for (i in seq_along(unique_exposures)) {
  current_exposure <- unique_exposures[i]
  current_exposure_clean <- clean_filename(current_exposure)

  cat("Step", step_count, ": Processing exposure", current_exposure,
      "(Progress:", i, "/", n_exposures, ")\n")
  step_count <- step_count + 1

  # 筛选当前暴露的数据子集
  exposure_subset <- exposure_dat[exposure_dat$exposure == current_exposure, ]

  # 检查数据是否为空

  if (nrow(exposure_subset) == 0) {
    warning(paste("Warning: Exposure", current_exposure, "has no data. Skipping!"))
    setTxtProgressBar(progress_bar, i)
    next
  }

  # 数据协调（对齐等位基因方向）
  outcome_dat$outcome <- disease_name
  harmonised_dat <- harmonise_data(exposure_subset, outcome_dat)

  # 检查协调后数据是否为空
  if (nrow(harmonised_dat) == 0) {
    warning(paste("Warning: Harmonised data for exposure", current_exposure, "is empty. Skipping!"))
    setTxtProgressBar(progress_bar, i)
    next
  }

  # 检查SNP数量是否足够
  if (nrow(harmonised_dat) < 2) {
    warning(paste("Warning: Not enough SNPs for MR analysis of", current_exposure, ". Skipping!"))
    setTxtProgressBar(progress_bar, i)
    next
  }

  # 根据结局p值过滤（排除可能的反向因果，保留p > 5e-06）
  filtered_dat <- harmonised_dat[harmonised_dat$pval.outcome > 5e-06, ]
  if (nrow(filtered_dat) < 1) {
    warning(paste("Warning: No data left after filtering for exposure", current_exposure, ". Skipping!"))
    setTxtProgressBar(progress_bar, i)
    next
  }

  # ---------------------------------------------------------------------------
  # MR分析
  # ---------------------------------------------------------------------------
  mr_result <- tryCatch({
    mr(filtered_dat)
  }, error = function(e) {
    warning(paste("MR分析出错，暴露", current_exposure, "跳过：", e$message))
    return(NULL)
  })

  if (is.null(mr_result)) {
    setTxtProgressBar(progress_bar, i)
    next
  }

  # 生成OR结果
  or_result <- tryCatch({
    generate_odds_ratios(mr_result)
  }, error = function(e) {
    warning(paste("生成OR结果出错，暴露", current_exposure, "跳过：", e$message))
    return(NULL)
  })

  if (!is.null(or_result)) {
    results_or <- rbind(results_or, or_result)
  }

  # 异质性检验
  heterogeneity_result <- tryCatch({
    mr_heterogeneity(filtered_dat)
  }, error = function(e) {
    warning(paste("异质性检验出错，暴露", current_exposure, "跳过：", e$message))
    return(NULL)
  })

  if (!is.null(heterogeneity_result)) {
    results_heterogeneity <- rbind(results_heterogeneity, heterogeneity_result)
  }

  # 多效性检验（MR-Egger截距检验）
  pleiotropy_result <- tryCatch({
    mr_pleiotropy_test(filtered_dat)
  }, error = function(e) {
    warning(paste("多效性检验出错，暴露", current_exposure, "跳过：", e$message))
    return(NULL)
  })

  if (!is.null(pleiotropy_result)) {
    results_pleiotropy <- rbind(results_pleiotropy, pleiotropy_result)
  }

  # Steiger方向性检验
  # 原理：比较SNP解释暴露vs结局的方差(R²)
  # correct_causal_direction = TRUE 表示因果方向正确（暴露→结局）
  steiger_result <- tryCatch({
    directionality_test(filtered_dat)
  }, error = function(e) {
    warning(paste("Steiger方向性检验出错，暴露", current_exposure, "跳过：", e$message))
    return(NULL)
  })

  if (!is.null(steiger_result)) {
    results_steiger <- rbind(results_steiger, steiger_result)
  }

  # 更新进度条
  setTxtProgressBar(progress_bar, i)
}

# 关闭进度条
close(progress_bar)

# =============================================================================
# 结果整理与保存
# =============================================================================

# 为OR结果添加暴露名称
if (nrow(results_or) > 0) {
  results_or_annotated <- merge(
    results_or,
    exposure_id_mapping,
    by.x = "exposure",
    by.y = "id.exposure",
    all.x = TRUE
  )

  # 重命名列
  colnames(results_or_annotated)[colnames(results_or_annotated) == "exposure"] <- "id.exposure"
  colnames(results_or_annotated)[colnames(results_or_annotated) == "exposure.y"] <- "exposure_name"

  # 调整列顺序，将exposure_name放到最后
  col_names <- colnames(results_or_annotated)
  if ("exposure_name" %in% col_names) {
    exposure_name_idx <- which(col_names == "exposure_name")
    other_cols <- col_names[-exposure_name_idx]
    results_or_annotated <- results_or_annotated[, c(other_cols, "exposure_name")]
  }

  results_or_final <- results_or_annotated
} else {
  results_or_final <- results_or
}

# 保存OR结果
write.csv(results_or_final, file = "all_odds_ratios.csv", row.names = FALSE)

# 保存异质性检验结果
write.csv(results_heterogeneity, file = "all_heterogeneity_results.csv", row.names = FALSE)

# 保存多效性检验结果
write.csv(results_pleiotropy, file = "all_pleiotropy_results.csv", row.names = FALSE)

# 保存Steiger方向性检验结果
if (nrow(results_steiger) > 0) {
  # 添加细胞名称
  results_steiger_annotated <- merge(
    results_steiger,
    exposure_id_mapping,
    by.x = "exposure",
    by.y = "id.exposure",
    all.x = TRUE
  )

  # 重命名列
  colnames(results_steiger_annotated)[colnames(results_steiger_annotated) == "exposure"] <- "id.exposure"
  colnames(results_steiger_annotated)[colnames(results_steiger_annotated) == "exposure.y"] <- "cell_name"

  write.csv(results_steiger_annotated, file = "all_steiger_results.csv", row.names = FALSE)
  cat("Steiger方向性检验结果已保存至: all_steiger_results.csv\n")

  # 统计方向正确的比例
  n_correct <- sum(results_steiger_annotated$correct_causal_direction == TRUE, na.rm = TRUE)
  n_total <- sum(!is.na(results_steiger_annotated$correct_causal_direction))
  cat("方向性检验：", n_correct, "/", n_total, "个暴露的因果方向正确\n")
} else {
  write.csv(results_steiger, file = "all_steiger_results.csv", row.names = FALSE)
  cat("警告: 没有Steiger方向性检验结果\n")
}

cat("Processing and analysis for all exposures are complete! Results are saved in folder:", output_dir, "\n")
