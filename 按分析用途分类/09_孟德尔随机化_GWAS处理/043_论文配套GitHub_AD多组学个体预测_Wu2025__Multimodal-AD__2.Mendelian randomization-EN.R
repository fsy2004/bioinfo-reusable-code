# =============================================================================
# 编号       : R043
# 脚本名     : 2.Mendelian randomization-EN.R
# 分类       : 09_孟德尔随机化_GWAS处理
# 项目来源   : 论文配套GitHub_AD多组学个体预测_Wu2025
# 用途       : AD GWAS VCF→TwoSampleMR 数据格式转换 + 孟德尔随机化主分析（IVW、MR-Egger、加权中位数、敏感性分析、LOO、漏斗、曼哈顿、QQ）。
# 结果图     : 散点图；森林图；漏斗图；Leave-one-out敏感性图；曼哈顿图；QQ图
# 主要 R 包  : VariantAnnotation; gwasglue; TwoSampleMR; KYNMCZRJ; dplyr
# 整理时间   : 2026-05-10
# =============================================================================
#######VCF文件转换############
# 加载必要库
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(gwasglue)
  library(TwoSampleMR)
  library(KYNMCZRJ)
  library(dplyr)
})

#----------- 第二部分：用户配置 ------------
# 重要参数配置（根据实际情况修改！）
config <- list(
  work_dir = "D:/AD_GWAS",     # 工作路径（使用正斜杠）
  vcf_file = "ebi-a-GCST90027158.vcf",                       # VCF文件名
  output_basename = "my_analysis",                   # 输出文件前缀
  manhattan_threshold = 5e-8,                        # 曼哈顿图显著阈值
  manhattan_cex = 0.5,                               # 点大小（0.1-1.0）
  plot_width = 14,                                   # 图像宽度（英寸）
  plot_height = 8                                    # 图像高度（英寸）
)

#----------- 第三部分：函数定义 ------------
# 安全文件写入函数（自动创建目录）
safe_write <- function(data, filename, ...) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  message("\n正在保存文件: ", filename)
  write.table(data, file = filename, ...)
}

#----------- 第四部分：主程序 ------------
tryCatch({
  # 设置工作目录
  setwd(config$work_dir)
  message("\n当前工作目录: ", getwd())
  
  # 检查输入文件是否存在
  if (!file.exists(config$vcf_file)) {
    stop("错误：VCF文件不存在！请检查路径: ", config$vcf_file)
  }
  
  # 读取VCF文件（显示进度条）
  message("\n⏳ 正在读取VCF文件...")
  vcf_data <- readVcf(config$vcf_file)
  
  # 转换为TwoSampleMR格式
  message("🔄 转换数据格式...")
  mr_data <- gwasvcf_to_TwoSampleMR(vcf_data)
  
  
  
  
  #----------- 第五部分：可视化 ------------
  # 准备绘图数据
  plot_data <- mr_data %>%
    select(SNP, chr.exposure, pos.exposure, pval.exposure) %>%
    rename(SNP = SNP, CHR = chr.exposure, BP = pos.exposure, pvalue = pval.exposure)
  
  # 生成曼哈顿图（线性布局）
  message("\n🎨 生成曼哈顿图...")
  CMplot(plot_data,
         plot.type = "m",
         LOG10 = TRUE,
         threshold = config$manhattan_threshold,
         threshold.col = "red",
         threshold.lwd = 2,
         amplify = TRUE,  # 自动放大显著点
         cex = config$manhattan_cex,
         ylim = c(0, 60),  # 调整Y轴范围
         file = "pdf",
         file.name = paste0(config$output_basename, "_Manhattan"),
         width = config$plot_width,
         height = config$plot_height)
  
  # 生成曼哈顿圈图
  message("🎨 生成圈型曼哈顿图...")
  CMplot(plot_data,
         plot.type = "c",
         cir.chr.h = 1.5,  # 染色体标签高度
         cir.legend.cex = 0.8,
         file = "pdf",
         file.name = paste0(config$output_basename, "_Circular"),
         width = 10,
         height = 10)
  
  message("\n✅ 所有分析已完成！输出文件前缀: ", config$output_basename)
  
}, error = function(e) {
  message("\n❌ 运行出错: ", e$message)
  if(grepl("VariantAnnotation", e$message)) {
    message("提示：请检查是否已安装VariantAnnotation包：BiocManager::install('VariantAnnotation')")
  }
})

# 安装并加载 dplyr 包
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# 确认 mr_data 是一个有效的数据框并显示其列名
print(class(mr_data))
print(colnames(mr_data))

# 删除指定的列
mr_data_clean <- dplyr::select(mr_data, -c(exposure, mr_keep.exposure, 
                                           pval_origin.exposure, id.exposure, 
                                           ncase.exposure, ncontrol.exposure))

# 检查结果
print(mr_data_clean)


######正确定义列名映射#####
custom_colnames <- c(
  CHR = "chr.exposure",
  BP = "pos.exposure",
  other_allele = "other_allele.exposure",
  effect_allele = "effect_allele.exposure",
  beta = "beta.exposure",
  se = "se.exposure",
  pval = "pval.exposure",
  eaf = "eaf.exposure",
  samplesize = "samplesize.exposure"
)

missing_cols <- setdiff(unname(custom_colnames), colnames(mr_data_clean))
if (length(missing_cols) > 0) {
  warning("以下原始列不存在，无法重命名：\n", paste(missing_cols, collapse = ", "))
}
library(dplyr)


######验证结果#########
print(colnames(mr_data_renamed))
######保存结果########
output_filename <- "mr_data_processed.txt"  # 自定义文件名
write.table(
  mr_data_renamed,
  file = output_filename,
  sep = "\t",          # 制表符分隔
  row.names = FALSE,   # 不保存行名
  quote = FALSE,       # 禁用引号包裹字符
  na = "NA"            # 缺失值标记为NA
)
# 验证文件保存
message("\n✅ 数据已保存至: ", normalizePath(output_filename))


library(KYNMCZRJ)
library(TwoSampleMR)

# 设置工作路径和读取文件
FileNames <- list.files(getwd(), pattern = ".csv")
exp_dat_ids <- FileNames
exps <- FileNames

# 读取并处理结局数据
outcome_file <- "AD_GCST90027158.txt"
out <- fread(outcome_file, header = TRUE)
out$trait <- 'AD'
outcomeid <- out
head(outcomeid)

# 创建存放结果的文件夹
output_dir <- "mendelian_test"
dir.create(path = output_dir)

# 定义辅助函数
get_f_noeaf <- function(dat, F_value = 10) {
  if(is.null(dat$beta.exposure[1]) || is.na(dat$beta.exposure[1])) {
    print("数据不包含beta，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$se.exposure[1]) || is.na(dat$se.exposure[1])) {
    print("数据不包含se，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$samplesize.exposure[1]) || is.na(dat$samplesize.exposure[1])) {
    print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)
  }
  
  R2 <- (dat$beta.exposure^2) / ((dat$se.exposure^2 * dat$samplesize.exposure) + dat$beta.exposure^2)
  F <- (dat$samplesize.exposure - 2) * R2 / (1 - R2)
  dat$R2 <- R2
  dat$F <- F
  dat <- subset(dat, F > F_value)
  return(dat)
}

get_f <- function(dat, F_value = 10) {
  log <- is.na(dat$eaf.exposure)
  log <- unique(log)
  if(length(log) == 1 && log == TRUE) {
    print("数据不包含eaf，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$beta.exposure[1]) || is.na(dat$beta.exposure[1])) {
    print("数据不包含beta，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$se.exposure[1]) || is.na(dat$se.exposure[1])) {
    print("数据不包含se，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$samplesize.exposure[1]) || is.na(dat$samplesize.exposure[1])) {
    print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)
  }
  
  if("FALSE" %in% log) {
    R2 <- (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) / 
      ((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) + 
         (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$se.exposure^2) * dat$samplesize.exposure))
    F <- (dat$samplesize.exposure - 2) * R2 / (1 - R2)
    dat$R2 <- R2
    dat$F <- F
    dat <- subset(dat, F > F_value)
    return(dat)
  }
}

steiger_test <- function(dat) {
  dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
  dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
  res_steiger <- mr_steiger(
    p_exp = dat$pval.exposure,
    p_out = dat$pval.outcome,
    n_exp = dat$samplesize.exposure,
    n_out = dat$samplesize.outcome,
    r_exp = dat$r.exposure,
    r_out = dat$r.outcome
  )
  res_steiger <- directionality_test(dat)
  
  return(res_steiger)
}

results_binary <- function(N, alpha, R2xz, K, OR, epower) {
  threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
  f.value <- 1 + N * R2xz / (1 - R2xz)
  
  if (is.na(epower)) {
    b_MR <- K * (OR / (1 + K * (OR - 1)) - 1)
    v_MR <- (K * (1 - K) - b_MR^2) / (N * R2xz)
    NCP <- b_MR^2 / v_MR
    
    power <- 1 - pchisq(threschi, 1, NCP)
    data.frame(Parameter = c("Power", "NCP", "F-statistic"), Value = c(power, NCP, f.value), Description = c("", "Non-Centrality-Parameter", "The strength of the instrument"))    
  } else {
    z1 <- qnorm(1 - alpha / 2)
    z2 <- qnorm(epower)
    Z <- (z1 + z2)^2
    
    b_01 <- K * (OR / (1 + K * (OR - 1)) - 1)
    f <- K * (1 - K) - b_01^2
    N1 <- Z * f / (b_01^2 * R2xz)
    N1 <- ceiling(N1)
    data.frame(Parameter = "Sample Size", Value = N1)
  }
}

choose_MR <- function(dat = dat) {
  res_hete <- NULL  
  if (nrow(dat) < 3) {
    res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
  } else {
    res_hete <- mr_heterogeneity(dat)
    if (res_hete$Q_pval[2] < 0.05) {
      res <- mr(dat, method_list = c(
        "mr_egger_regression", "mr_weighted_median", "mr_ivw_mre", "mr_weighted_mode", "mr_simple_mode"
      ))
    } else {
      res <- mr(dat, method_list = c(
        "mr_egger_regression", "mr_weighted_median", "mr_ivw_fe", "mr_weighted_mode", "mr_simple_mode"
      ))
    }
  }
  AAA <- list(res_hete, res)
  return(list(AAA))
}

# 主函数：循环处理每个暴露数据文件
process_exposure_data <- function(exp_dat_id, exp, outcomeid) {
  d3 <- try(fread(paste0(getwd(), "/", exp_dat_id), fill = TRUE), silent = TRUE)
  d3 <- subset(d3, d3$pval < 1e-5)
  
  if (nrow(d3) == 0) return(NULL)
  
  d3 <- format_data(as.data.frame(d3), type = "exposure")
  
  d4 <- ld_clump(
    clump_kb = 10000,
    clump_r2 = 0.01,
    pop = "EUR",
    dplyr::tibble(rsid = d3$SNP, pval = d3$pval.exposure, id = d3$id.exposure),
    plink_bin = "D:/AD_GWAS/参考文件/本地plink/plink_win64_20230116/plink.exe",
    bfile = "D:/AD_GWAS/参考文件/本地plink/EUR/EUR"
  )
  
  exp_data <- subset(d3, SNP %in% d4$rsid)
  if (nrow(exp_data) == 0) return(NULL)
  
  outcome_dat <- merge(exp_data, outcomeid, by.x = "SNP", by.y = "SNP")
  if (nrow(outcome_dat) == 0) return(NULL)
  
  write.csv(outcome_dat, file = "d.csv")
  out_data <- read_outcome_data(
    snps = exp_data$SNP,
    filename = "d.csv",
    sep = ","
  )
  
  out_data <- subset(out_data, pval.outcome > 5e-8)
  
  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_data,
    outcome_dat = out_data
  )
  
  dat <- subset(dat, mr_keep == TRUE)
  dat <- get_f(dat, F_value = 10)
  
  res <- choose_MR(dat = dat)
  res1 <- generate_odds_ratios(res[[1]][[2]])
  
  res1$estimate <- paste0(
    format(round(res1$or, 2), nsmall = 2), " (", 
    format(round(res1$or_lci95, 2), nsmall = 2), "-",
    format(round(res1$or_uci95, 2), nsmall = 2), ")"
  )
  
  openxlsx::write.xlsx(dat, file = paste0(output_dir, "/", exp, "-dat.xlsx"), rowNames = FALSE)
  openxlsx::write.xlsx(res1, paste0(output_dir, "/", exp, "-res.xlsx"))
  
  res_steiger <- steiger_test(dat)
  
  N <- dat$samplesize.outcome[1]
  alpha <- 0.05
  R2xz <- sum(dat$R2)
  K <- (dat$ncase.outcome[1] / dat$ncontrol.outcome[1])
  OR <- if (nrow(dat) == 1) {
    res1 %>% filter(method == "Wald ratio") %>% pull(or)
  } else {
    res1 %>% filter(grepl("Inverse variance weighted", method)) %>% pull(or)
  }
  
  epower <- NA
  power <- results_binary(N, alpha, R2xz, K, OR, epower)
  
  res3 <- res1[, -c(10:14)]
  res4 <- tidyr::pivot_wider(
    res3, names_from = "method", names_vary = "slowest",
    values_from = c("b", "se", "pval", "estimate")
  )
  
  colnames(res4) <- gsub("\\(.*\\)", "", colnames(res4))
  
  res_steiger2 <- dplyr::select(res_steiger, correct_causal_direction, steiger_pval)
  power2 <- tidyr::pivot_wider(
    power, names_from = "Parameter", names_vary = "slowest",
    values_from = c("Value", "Description")
  )[, 1]
  
  res_ALL <- cbind(res4, res_steiger2, power2)
  res_ALL$F <- mean(dat$F, na.rm = TRUE)
  res_ALL$R2 <- sum(dat$R2)
  
  if (nrow(dat) <= 2) {
    write.csv(res_ALL, file = paste0(output_dir, "/", exp, "1.csv"), row.names = FALSE)
  } else {
    res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    res_leaveone <- mr_leaveoneout(dat)
    
    p1 <- mr_scatter_plot(res[[1]][[2]], dat)
    pdf(paste0(output_dir, "/", exp, "_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    
    res_single <- mr_singlesnp(dat, all_method)
    p2 <- mr_forest_plot(res_single)
    pdf(paste0(output_dir, "/", exp, "_forest.pdf"))
    print(p2[[1]])
    dev.off()
    
    p3 <- mr_funnel_plot(res_single)
    pdf(paste0(output_dir, "/", exp, "_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    
    res_loo <- mr_leaveoneout(dat)
    pdf(paste0(output_dir, "/", exp, "_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
    
    res_hete <- purrr::map(.x = seq_along(res), .f = ~res[[.x]][[1]])
    res_hete <- do.call(rbind, res_hete)
    res_hete2 <- tidyr::pivot_wider(
      res_hete, names_from = "method", names_vary = "slowest",
      values_from = c("Q", "Q_df", "Q_pval")
    )[, 4:6]
    
    res_plei2 <- dplyr::select(res_plei, egger_intercept, se, pval)
    
    res_ALL <- cbind(res_ALL, res_hete2, res_plei2)
    write.csv(res_ALL, file = paste0(output_dir, "/", exp, ".csv"), row.names = FALSE)
  }
}

# 循环处理所有暴露数据文件
for (qaq in 1:length(exp_dat_ids)) {
  process_exposure_data(exp_dat_ids[qaq], exps[qaq], outcomeid)
}

# 获取所有以 .csv1.csv 和 .csv.csv 结尾的文件
csv_files <- list.files(output_dir, pattern = "\\.csv1\\.csv$|\\.csv\\.csv$", full.names = TRUE)

# 读取第一个文件
combined_df <- read.csv(csv_files[1], stringsAsFactors = FALSE)

# 处理 global_test_p 列
if (!is.null(combined_df$global_test_p)) {
  combined_df$global_test_p <- as.character(combined_df$global_test_p)
}

# 循环读取剩余的文件并合并
for (i in 2:length(csv_files)) {
  temp_df <- read.csv(csv_files[i], stringsAsFactors = FALSE)
  
  if (!is.null(temp_df$global_test_p)) {
    temp_df$global_test_p <- as.character(temp_df$global_test_p)
    temp_df$global_test_p[temp_df$global_test_p == "<"] <- ""
  }
  
  combined_df <- bind_rows(combined_df, temp_df)
}

# 写入最终合并的CSV文件
write.csv(combined_df, "imm_cell-ad-1.csv", row.names = FALSE)


######森林图##########
# 加载必要的包
library(forestploter)
library(grid)

# 读取数据（请替换为您的文件路径）
dt <- read.csv("imm_cell-ad-1.CSV", check.names = FALSE)

# 重命名列名（根据实际列名调整）
names(dt) <- gsub("\\.$", "", names(dt))  # 去除列名末尾的句号（如果存在）

# 计算OR和置信区间
dt$OR <- exp(dt$b_Inverse.variance.weighted)
dt$lower <- exp(dt$b_Inverse.variance.weighted - 1.96 * dt$se_Inverse.variance.weighted)
dt$upper <- exp(dt$b_Inverse.variance.weighted + 1.96 * dt$se_Inverse.variance.weighted)

# 确保 b_Inverse.variance.weighted 和 se_Inverse.variance.weighted 为数值类型
dt$b_Inverse.variance.weighted <- as.numeric(gsub("[^0-9.-]", "", dt$b_Inverse.variance.weighted))
dt$se_Inverse.variance.weighted <- as.numeric(gsub("[^0-9.-]", "", dt$se_Inverse.variance.weighted))

# 检查是否存在 NA 值
#dt <- dt[complete.cases(dt), ]  # 去除含有 NA 的行

# 计算 OR 和置信区间
dt$OR <- exp(dt$b_Inverse.variance.weighted)  # 计算 OR
dt$lower <- exp(dt$b_Inverse.variance.weighted - 1.96 * dt$se_Inverse.variance.weighted)  # 下限
dt$upper <- exp(dt$b_Inverse.variance.weighted + 1.96 * dt$se_Inverse.variance.weighted)  # 上限

# 创建展示文本列
dt$`OR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", dt$OR, dt$lower, dt$upper)
dt$`p-value` <- ifelse(dt$pval_Inverse.variance.weighted < 0.05, "<0.001",
                       sprintf("%.3f", dt$pval_Inverse.variance.weighted))

# 创建空白列用于绘制森林图
dt$` ` <- paste(rep(" ", 20), collapse = " ")

# 选择需要展示的列
show_cols <- c("exposure", "nsnp", "OR (95% CI)", "p-value", " ")

# 设置主题样式
tm <- forest_theme(
  base_size = 10,
  ci_pch = 15,
  ci_col = "#1B9E77",
  ci_lwd = 1.5,
  ci_Theight = 0.2,
  refline_lwd = 1,
  refline_lty = "dashed",
  refline_col = "grey20",
  summary_fill = "#7570B3",
  summary_col = "#7570B3",
  footnote_cex = 0.8
)

# 生成森林图
p <- forest(dt[, show_cols],
            est = dt$OR,
            lower = dt$lower,
            upper = dt$upper,
            ci_column = 5,          # 置信区间在第五列
            ref_line = 1,           # 参考线在OR=1处
            x_trans = "log",        # X轴对数转换
            xlim = c(0.5, 2),       # X轴范围
            ticks_at = c(0.5, 1, 1.5, 2), # 刻度位置
            arrow_lab = c("Protective", "Risk"),
            footnote = "Inverse Variance Weighted Method",
            theme = tm)

# 调整边距和保存为高分辨率图片（通过调整mar参数来增加上下边距）
png("MR_forestplot.png", width = 3000, height = 4000, res = 300, 
    units = "px", bg = "white")
par(mar = c(5, 5, 4, 4))  # 设置边距，底部、左侧、顶部、右侧

# 打印并保存图形
print(p)

# 关闭设备
dev.off()



# 保存为 PDF 文件
pdf("免疫细胞MR_forestplot.pdf", width = 12, height = 32)  # 设定PDF画幅宽高
par(mar = c(5, 5, 4, 4))  # 设置边距，底部、左侧、顶部、右侧

# 打印并保存图形
print(p)

# 关闭设备
dev.off()  # 关闭设备以保存文件

# 首先，计算 OR 和置信区间



###########火山图##########
#引用包
library(dplyr)
library(ggplot2)
library(ggrepel)
# 读取数据（请确保路径正确）
dt <- read.csv("me_cell-ad.csv", check.names = FALSE)

# 重命名列名（根据实际列名调整，去掉末尾的句号）
names(dt) <- gsub("\\.$", "", names(dt))  # 去除列名末尾的句号（如果存在）

# 打印列名以检查重复
print(names(dt))

# 确保 b_Inverse.variance.weighted 和 se_Inverse.variance.weighted 为数值类型
dt$b_Inverse.variance.weighted <- as.numeric(gsub("[^0-9.-]", "", dt$b_Inverse.variance.weighted))
dt$se_Inverse.variance.weighted <- as.numeric(gsub("[^0-9.-]", "", dt$se_Inverse.variance.weighted))

# 计算 OR 和置信区间
dt$OR <- exp(dt$b_Inverse.variance.weighted)
dt$lower <- exp(dt$b_Inverse.variance.weighted - 1.96 * dt$se_Inverse.variance.weighted)
dt$upper <- exp(dt$b_Inverse.variance.weighted + 1.96 * dt$se_Inverse.variance.weighted)

# 创建展示文本列
dt$`OR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", dt$OR, dt$lower, dt$upper)
dt$`p-value` <- ifelse(dt$pval_Inverse.variance.weighted < 0.001, "<0.001",
                       sprintf("%.3f", dt$pval_Inverse.variance.weighted))

# 计算 log2(OR) 用于 x 轴（以 1 为中心）
dt$log2_OR <- log2(dt$OR)

# 确保 p-value 存在有效的数值
dt$pvalue <- as.numeric(dt$pval_Inverse.variance.weighted)

# 检查并修复重复列
library(stringr)

# 检查重复列名
duplicated_cols <- names(dt)[duplicated(names(dt))]

if (length(duplicated_cols) > 0) {
  cat("发现重复列名: ", unique(duplicated_cols), "\n")
  
  # 为重复列名添加后缀
  for (col in unique(duplicated_cols)) {
    # 确定该列的所有实例
    indices <- which(names(dt) == col)
    
    # 如果有多个，重命名后面的列
    if (length(indices) > 1) {
      for (i in 2:length(indices)) {
        names(dt)[indices[i]] <- paste0(col, "_dup", i - 1)
      }
    }
  }
}

# 打印新列名以确认
print(names(dt))

# 确保 b_Inverse.variance.weighted 和 se_Inverse.variance.weighted 为数值类型
dt$b_Inverse.variance.weighted <- as.numeric(gsub("[^0-9.-]", "", dt$b_Inverse.variance.weighted))
dt$se_Inverse.variance.weighted <- as.numeric(gsub("[^0-9.-]", "", dt$se_Inverse.variance.weighted))

# 检查是否有 NA 值
cat("b_Inverse.variance.weighted 中 NA 的数量:", sum(is.na(dt$b_Inverse.variance.weighted)), "\n")
cat("se_Inverse.variance.weighted 中 NA 的数量:", sum(is.na(dt$se_Inverse.variance.weighted)), "\n")

# 计算 OR 和置信区间
dt$OR <- exp(dt$b_Inverse.variance.weighted)
dt$lower <- exp(dt$b_Inverse.variance.weighted - 1.96 * dt$se_Inverse.variance.weighted)
dt$upper <- exp(dt$b_Inverse.variance.weighted + 1.96 * dt$se_Inverse.variance.weighted)

# 创建展示文本列
dt$`OR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", dt$OR, dt$lower, dt$upper)
dt$`p-value` <- ifelse(dt$pval_Inverse.variance.weighted < 0.001, "<0.001",
                       sprintf("%.3f", dt$pval_Inverse.variance.weighted))

# 计算 log2(OR) 用于 x 轴（以 1 为中心）
dt$log2_OR <- log2(dt$OR)


library(ggplot2)
library(ggrepel)

# 创建一个新的颜色列，根据条件进行赋值
dt$color <- ifelse(dt$log2_OR > 0 & dt$pvalue < 0.05, "#FA8260", 
                   ifelse(dt$log2_OR < 0 & dt$pvalue < 0.05, "#4D8FD1", "#ABABA6"))

# 绘制火山图
p <- ggplot(dt, aes(x = log2_OR, y = -log10(pvalue), color = color, size = -log10(pvalue))) +
  geom_point() + 
  xlab("log2(OR)") + ylab("-log10(P-value)") + 
  xlim(-0.3, 0.3) +
  scale_color_identity() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", alpha = 0.5) +
  scale_size_continuous(range = c(0.1, 4)) + 
  labs(title = "Volcano Plot (IVW Method)") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 12),
        axis.title = element_text(size = 12))

# 在图形中标注显著基因数据
showData <- dt[dt$pvalue < 0.05, ]  # 筛选显著基因数据用于标注
p1 <- p + geom_text_repel(data = showData,
                          box.padding = 0.2, point.padding = 0.2, min.segment.length = 0.01,
                          size = 3, aes(label = exposure),
                          max.overlaps = 20)
p1


# 输出图形
pdf(file = "volcano_OR_IVW_formatted2.pdf", width = 8, height = 6)
print(p)
dev.off()

# 输出图形
pdf(file = "volcano_OR_IVW1.pdf", width = 7, height = 5.5)
print(p1)
dev.off()


#######代谢物富集分析#######
library(ggplot2)
library(readr)

# 读取数据
dt <- read.csv("msea_ora_result.csv", check.names = FALSE)

# 检查数据结构
print(head(dt))


# 确保 FDR 列为数值类型
dt$FDR <- as.numeric(dt$FDR)

# 按 FDR 排序并选择前五个结果
top5 <- dt[order(dt$FDR), ][1:10, ]
# 打印前五个结果以确认
print(top5)
# 载入所需的库
library(ggplot2)
library(readr)

# 读取数据
dt <- read.csv("msea_ora_result.csv", check.names = FALSE)

# 确保 FDR 列为数值类型
dt$FDR <- as.numeric(dt$FDR)

# 按 FDR 排序并选择前五个结果
top5 <- dt[order(dt$FDR), ][1:10, ]

# 绘制气泡图
p <- ggplot(top5, aes(x = -log10(FDR), y = reorder(Term, FDR))) +
  geom_point(aes(size = Enrichment_Ratio, color = FDR), alpha = 0.7) +
  scale_color_gradientn(colors = c("#A50026","#D73027", "#F46D43","#FDAE61" ), 
                        values = scales::rescale(c(0, 0.5, 1, 1.5)), 
                        guide = "colourbar") +  # 设置颜色渐变
  scale_size(range = c(1, 10), name = "Enrichment Ratio") +  # 调整气泡大小
  xlab("-log10(FDR)") +  # x 轴标签
  ylab("Biological Terms") +  # y 轴标签
  labs(title = "Genetic Susceptibility: Overview of Enriched Metabolite Sets")#+
#theme_minimal() +  # 使用简洁主题
#theme(plot.title = element_text(size = 16, face = "bold", color = "#d23927"),
#      axis.title.x = element_text(size = 14),
#     axis.title.y = element_text(size = 14),
#     legend.title = element_text(size = 12))
p
# 输出图形
pdf(file = "危险-bubble_plot_genetic_susceptibility6.pdf", width = 8, height = 6)
print(p)
dev.off()


library(ggplot2)
library(readr)

# 读取数据
dt <- read.csv("msea_ora_result.csv", check.names = FALSE)

# 检查数据结构
print(head(dt))


# 确保 FDR 列为数值类型
dt$FDR <- as.numeric(dt$FDR)

# 按 FDR 排序并选择前五个结果
top5 <- dt[order(dt$FDR), ][1:10, ]
# 打印前五个结果以确认
print(top5)
# 载入所需的库
library(ggplot2)
library(readr)

# 读取数据
dt <- read.csv("msea_ora_result.csv", check.names = FALSE)

# 确保 FDR 列为数值类型
dt$FDR <- as.numeric(dt$FDR)

# 按 FDR 排序并选择前五个结果
top5 <- dt[order(dt$FDR), ][1:10, ]

# 绘制气泡图
p <- ggplot(top5, aes(x = -log10(FDR), y = reorder(Term, FDR))) +
  geom_point(aes(size = Enrichment_Ratio, color = FDR), alpha = 0.7) +
  scale_color_gradientn(colors = c("#313695", "#4575B4","#74ADD1", "#ABD9E9" ), 
                        values = scales::rescale(c(0, 0.5, 1, 1.5)), 
                        guide = "colourbar") +  # 设置颜色渐变
  scale_size(range = c(1, 10), name = "Enrichment Ratio") +  # 调整气泡大小
  xlab("-log10(FDR)") +  # x 轴标签
  ylab("Biological Terms") +  # y 轴标签
  labs(title = "Genetic Susceptibility: Overview of Enriched Metabolite Sets") #+
#theme_minimal() +  # 使用简洁主题
#theme(plot.title = element_text(size = 16, face = "bold", color = "#313695"),
#      axis.title.x = element_text(size = 14),
#     axis.title.y = element_text(size = 14),
#    legend.title = element_text(size = 12))
p 


# 输出图形
pdf(file = "保护bubble_plot_genetic_susceptibility6.pdf", width = 8, height = 6)
print(p)
dev.off()



############提取免疫细胞#############  exposureID.txt和input.txt
setwd("D:\\AD_GWAS\\免疫细胞\\保护")
# 读取 exposureID.txt 文件
rt = read.table("exposureID.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# 将数据转换为矩阵，不进行数值转换
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = as.matrix(exp)

# 获取自噬基因表达量
gene = read.table("input.txt", header=FALSE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)

# 找到两个文件中的共同基因
sameGene = intersect(as.vector(gene[, 1]), rownames(data))

# 提取共同基因的表达量
geneExp = data[sameGene, ]
print(geneExp)
#输出结果
# 假设你已经在之前的步骤中提取了 geneExp

# 设置要保存的文件名称和路径
output_file_path <- "D:\\AD_GWAS\\免疫细胞\\保护\\gene_expression_output.txt"

# 将 geneExp 保存为 txt 格式
write.table(geneExp, file = output_file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# 读取 imm_cell_ad.txt 文件
imm_cell_data = read.table("imm_cell_ad.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# 将 geneExp 转换为数据框，并将行名作为一列
geneExp_df = as.data.frame(geneExp)
geneExp_df$Gene = rownames(geneExp_df)

# 合并数据框，使用 Gene 列作为键
# 假设 imm_cell_data 的第一列为基因名处理
merged_data = merge(geneExp_df, imm_cell_data, by.x = "Gene", by.y = imm_cell_data[,1], all = TRUE)

# 检查合并后的数据
print(merged_data)


library(ggplot2)
library(dplyr)

# 读取原始数据
df <- read.csv("imm_cell-ad-2已分类 - 危险因素.csv", check.names = FALSE)
library(dplyr)
library(ggplot2)

# 读取数据
df <- read.csv("imm_cell-ad-2已分类 - 危险因素.csv", check.names = FALSE)

# 1. 汇总：按 id 合并（去掉 Ratio）
df_grouped <- df %>%
  group_by(id) %>%
  summarise(
    Count = n(),                      # 相同 id 的数量
    Pval  = mean(pval, na.rm = TRUE)  # 平均 P 值
  ) %>%
  arrange(Pval) %>%                   # 可读性：小 P 值在上
  mutate(id = factor(id, levels = id))

# 2. 绘图：X=Count, Y=id，大小=Count，颜色=Pval
p <- ggplot(df_grouped, aes(x = Count, y = id, size = Count, color = Pval)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(name = "Count", range = c(3, 12)) +  # 调整大小范围
  scale_color_gradient(name = "P-value", low = "red", high = "orange") +
  labs(x = "Count", y = "") +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

# 3. 导出 PDF
pdf(file = "危险因素_bubble_plot_Count_Pval.pdf", width = 8, height = 6)
print(p)
dev.off()


#######免疫保护因素#########

setwd('D:\\AD_GWAS\\免疫细胞\\保护')
library(dplyr)
library(ggplot2)

# 读取数据
df <- read.csv("imm_cell-ad-2已分类保护因素.csv", check.names = FALSE)

# 1. 汇总：按 id 合并（不计算 Ratio）
df_grouped <- df %>%
  group_by(id) %>%
  summarise(
    Count = n(),                      # 相同 id 的数量
    Pval  = mean(pval, na.rm = TRUE)  # 平均 P 值
  )

# 2. Y 轴按 P 值排序（从小到大）
df_grouped <- df_grouped %>%
  arrange(Pval) %>%
  mutate(id = factor(id, levels = id))

# 3. 绘制气泡图
#    - X 轴：Count
#    - Y 轴：id
#    - 点大小：Count
#    - 颜色：Pval（蓝色渐变）
p <- ggplot(df_grouped, aes(x = Count, y = id, size = Count, color = Pval)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(name = "Count", range = c(3, 12)) +  # 可根据需要调整点大小范围
  scale_color_gradientn(
    name = "P-value",
    colours = c("#313695", "#4575B4", "#74ADD1")
  ) +
  labs(x = "Count", y = "") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

p
# 输出图形
pdf(file = "保护因素免疫细胞bubble_plot_6.pdf", width = 5, height = 4)
print(p)
dev.off()