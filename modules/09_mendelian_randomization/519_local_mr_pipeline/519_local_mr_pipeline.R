# =============================================================================
# 编号   : 519
# 脚本名 : 纯本地孟德尔随机化(MR)全流程 — 零 OpenGWAS API
# 分类   : 09_mendelian_randomization
# 用途   : 用本地下载的 GWAS summary stats 跑完整 MR:
#          取工具变量(p<5e-8 + F>10) → 本地 LD clumping → harmonise →
#          IVW/Egger/加权中位数/加权众数 + MR-PRESSO + 异质性/多效性/Steiger/留一,
#          全程不连 OpenGWAS API(规避 token 失效 / 限流 / 宕机)。
# 关键   : 本地 LD clumping = ieugwasr::ld_clump(plink_bin=, bfile=1000G EUR),不走 API。
# 依赖   : TwoSampleMR · ieugwasr · (plinkbinr 或 genetics.binaRies 提供 plink) · MRPRESSO · ggplot2
# 运行   : Rscript 519_local_mr_pipeline.R                         # 合成示例(无 bfile 则跳 clumping)
#          Rscript 519_local_mr_pipeline.R --exposure exp.csv --outcome out.tsv.gz --bfile D:/ref/EUR
# 输入   : exposure/outcome = 本地 GWAS summary 文件,列含 SNP/effect_allele/other_allele/eaf/beta/se/pval
#          (列名可在下方 read_*_data 调整);bfile = 1000G EUR plink 前缀(EUR.bed/bim/fam 去扩展名)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(TwoSampleMR); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  exposure = file.path(DDAT, "exposure_gwas.csv"),
  outcome  = file.path(DDAT, "outcome_gwas.csv"),
  bfile    = Sys.getenv("MR_LD_BFILE", unset = ""),   # 留空=跳过本地 LD clumping;服务器设此环境变量指向 1000G EUR plink 前缀
  outdir   = file.path(SCRIPT_DIR, "results"),
  p_thresh = 5e-8, F_thresh = 10, clump_r2 = 0.001, clump_kb = 10000))
for (k in c("p_thresh","F_thresh","clump_r2","clump_kb")) args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## 本地 plink 二进制(plinkbinr 或 genetics.binaRies,均离线,不连网) -------------
get_local_plink <- function() {
  if (requireNamespace("plinkbinr", quietly = TRUE))
    return(tryCatch(plinkbinr::get_plink_exe(), error = function(e) NA_character_))
  if (requireNamespace("genetics.binaRies", quietly = TRUE))
    return(tryCatch(genetics.binaRies::get_plink_binary(), error = function(e) NA_character_))
  p <- Sys.getenv("MR_PLINK_BIN", unset = NA_character_); if (!is.na(p) && nzchar(p)) p else NA_character_
}

## ---- 0. 合成示例(若无输入文件) ------------------------------------------
# exposure: 200 SNP,其中前 25 为强工具(大效应,p<5e-8),且对 outcome 有真因果(beta_out=0.30*beta_exp)
if (!(file.exists(args$exposure) && file.exists(args$outcome))) {
  set.seed(42); n <- 200
  snp <- sprintf("rs%07d", sample(1e6:9e6, n))
  ea  <- sample(c("A","C","G","T"), n, TRUE)
  oa  <- vapply(ea, function(x) sample(setdiff(c("A","C","G","T"), x), 1), character(1))
  eaf <- runif(n, 0.05, 0.5)
  is_iv    <- c(rep(TRUE, 25), rep(FALSE, n - 25))
  beta_exp <- ifelse(is_iv, rnorm(n, 0, 0.12), rnorm(n, 0, 0.02))
  se_exp   <- runif(n, 0.012, 0.020); pval_exp <- 2 * pnorm(-abs(beta_exp / se_exp))
  causal   <- 0.30
  se_out   <- runif(n, 0.012, 0.020)
  beta_out <- causal * beta_exp + rnorm(n, 0, 0.010); pval_out <- 2 * pnorm(-abs(beta_out / se_out))
  write.csv(data.frame(SNP = snp, effect_allele = ea, other_allele = oa, eaf = eaf,
    beta = beta_exp, se = se_exp, pval = pval_exp, samplesize = 50000, Phenotype = "ExposureX",
    chr = sample(1:22, n, TRUE), pos = sample(1e6:2.4e8, n)), args$exposure, row.names = FALSE)
  write.csv(data.frame(SNP = snp, effect_allele = ea, other_allele = oa, eaf = eaf,
    beta = beta_out, se = se_out, pval = pval_out, samplesize = 80000, Phenotype = "OutcomeY"),
    args$outcome, row.names = FALSE)
  cat(sprintf("[gen] 合成 GWAS: exposure %d SNP(25 真工具, 真因果=%.2f) + outcome\n", n, causal))
}

## ---- 1. 读 exposure + 筛工具变量(p<5e-8) + F 统计量 -----------------------
cat("Step 1: 读 exposure → 筛 p<5e-8 → 算 F...\n")
exp_all <- read_exposure_data(args$exposure, sep = ",", snp_col = "SNP", beta_col = "beta",
  se_col = "se", pval_col = "pval", effect_allele_col = "effect_allele",
  other_allele_col = "other_allele", eaf_col = "eaf", samplesize_col = "samplesize",
  phenotype_col = "Phenotype", chr_col = "chr", pos_col = "pos")
exp_iv <- exp_all[exp_all$pval.exposure < args$p_thresh, ]
exp_iv$F_stat <- (exp_iv$beta.exposure / exp_iv$se.exposure)^2
exp_iv <- exp_iv[exp_iv$F_stat > args$F_thresh, ]
cat(sprintf("  工具变量: %d 个 (p<%.0e & F>%g; 平均 F=%.1f)\n",
  nrow(exp_iv), args$p_thresh, args$F_thresh, mean(exp_iv$F_stat)))

## ---- 2. 本地 LD clumping(不走 OpenGWAS API) ------------------------------
cat("Step 2: 本地 LD clumping...\n")
plink_bin <- get_local_plink(); bed_ok <- file.exists(paste0(args$bfile, ".bed"))
if (!is.na(plink_bin) && bed_ok) {
  cl <- ieugwasr::ld_clump(
    dat = dplyr::tibble(rsid = exp_iv$SNP, pval = exp_iv$pval.exposure, id = exp_iv$exposure),
    plink_bin = plink_bin, bfile = args$bfile, clump_r2 = args$clump_r2, clump_kb = args$clump_kb)
  exp_iv <- exp_iv[exp_iv$SNP %in% cl$rsid, ]
  cat(sprintf("  本地 clumping 后保留 %d 个独立工具 (plink=%s, r2<%g, %gkb)\n",
    nrow(exp_iv), basename(plink_bin), args$clump_r2, args$clump_kb))
} else {
  cat("  ⚠ 未检测到本地 plink 或 bfile → 跳过 clumping。\n",
      "    真实数据务必配置(见 README): 装 plinkbinr/genetics.binaRies + 下 1000G EUR,\n",
      "    并用 --bfile 指向 EUR 前缀。合成示例已是独立 SNP,可继续。\n", sep = "")
}

## ---- 3. 读 outcome(仅提取工具 SNP)+ harmonise ---------------------------
cat("Step 3: 读 outcome → harmonise...\n")
out_dat <- read_outcome_data(snps = exp_iv$SNP, filename = args$outcome, sep = ",",
  snp_col = "SNP", beta_col = "beta", se_col = "se", pval_col = "pval",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  eaf_col = "eaf", samplesize_col = "samplesize")
out_dat$outcome <- "OutcomeY"
dat <- harmonise_data(exp_iv, out_dat); dat <- dat[dat$mr_keep, ]
write.csv(dat, file.path(args$outdir, "harmonised_instruments.csv"), row.names = FALSE)
cat(sprintf("  harmonise 后用于 MR 的工具: %d 个\n", nrow(dat)))

## ---- 4. MR 主分析 + 敏感性 ------------------------------------------------
cat("Step 4: MR(IVW/Egger/WM/WMode)+ 异质性/多效性/Steiger/MR-PRESSO...\n")
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
res_or <- generate_odds_ratios(res); write.csv(res_or, file.path(args$outdir, "MR_results.csv"), row.names = FALSE)
het  <- mr_heterogeneity(dat);     write.csv(het,  file.path(args$outdir, "MR_heterogeneity.csv"), row.names = FALSE)
plei <- mr_pleiotropy_test(dat);   write.csv(plei, file.path(args$outdir, "MR_pleiotropy.csv"),   row.names = FALSE)
stg  <- tryCatch(directionality_test(dat), error = function(e) NULL)
if (!is.null(stg)) write.csv(stg, file.path(args$outdir, "MR_steiger.csv"), row.names = FALSE)
presso_ok <- FALSE
try({
  if (requireNamespace("MRPRESSO", quietly = TRUE) && nrow(dat) >= 4) {
    pr <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure", data = as.data.frame(dat),
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000, SignifThreshold = 0.05)
    capture.output(print(pr), file = file.path(args$outdir, "MR_PRESSO.txt")); presso_ok <- TRUE
  }
}, silent = TRUE)
ivw <- res_or[res_or$method == "Inverse variance weighted", ]
cat(sprintf("  IVW OR=%.3f (%.3f-%.3f, p=%.2e) · Egger 截距 p=%.2f · 异质性 Q p=%.2e · PRESSO=%s\n",
  ivw$or, ivw$or_lci95, ivw$or_uci95, ivw$pval, plei$pval,
  het$Q_pval[het$method == "Inverse variance weighted"], presso_ok))

## ---- 5. 论文级敏感性图(scatter/forest/funnel/leave-one-out) ---------------
cat("Step 5: 出图...\n")
sng <- mr_singlesnp(dat); loo <- mr_leaveoneout(dat)
savep <- function(p, f, w = 6.5, h = 6) {
  ggsave(paste0(f, ".pdf"), p, width = w, height = h)
  ggsave(paste0(f, ".png"), p, width = w, height = h, dpi = 300) }
try(savep(mr_scatter_plot(res, dat)[[1]],        file.path(ASSETS, "MR_scatter")),     silent = TRUE)
try(savep(mr_forest_plot(sng)[[1]],              file.path(ASSETS, "MR_forest")),      silent = TRUE)
try(savep(mr_funnel_plot(sng)[[1]],              file.path(ASSETS, "MR_funnel")),      silent = TRUE)
try(savep(mr_leaveoneout_plot(loo)[[1]],         file.path(ASSETS, "MR_leaveoneout")), silent = TRUE)
cat("完成。结果表见", normalizePath(args$outdir), ";图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖版本快照(铁律6)
