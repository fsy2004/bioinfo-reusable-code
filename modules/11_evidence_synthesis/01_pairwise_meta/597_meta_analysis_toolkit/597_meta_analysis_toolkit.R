# =============================================================================
# 编号       : 597
# 脚本名     : 证据合成 / 系统评价 meta 分析(meta-analysis-toolkit 函数库演示)
# 分类       : 11_evidence_synthesis
# 用途       : 对 2x2 试验计数做完整 meta 分析:效应量 → 随机效应合并(REML+
#              Knapp-Hartung)→ 异质性 → 亚组 / meta 回归 → 发表偏倚 → 影响
#              诊断 → 森林图 / 漏斗图。统计全部由随模块捆绑的函数库
#              toolkit/R 完成(wrap metafor / meta / netmeta / mada /
#              bayesmeta / robvis / metasens),本脚本只做 I/O 与展示。
# 方法/包    : R · metafor(核心), meta / metadat / mada / bayesmeta / robvis
# 结果图     : 597_forest;597_funnel;597_metareg_bcg_ablat;597_influence_*
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 597_meta_analysis_toolkit.R
# 运行(自己): Rscript 597_meta_analysis_toolkit.R --input data/your.csv --outdir results/run1
# 输入规格 : CSV 需含列 author, year, tpos, tneg, cpos, cneg(2x2 计数);
#            可选列 ablat(连续调节变量)、alloc(分类调节变量)。数据格式与
#            metadat::dat.bcg 一致(见 example_data/bcg_trials.csv)。
# 整理日期 : 2026-08-14(自 fsy2004/meta-analysis-toolkit @ 3af629d 并入,
#            来源见本模块 README 与 toolkit/README.md)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "bcg_trials.csv"),
                      outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

## ---- 加载捆绑的 meta-analysis-toolkit 函数库 -------------------------------
TOOLKIT_R <- file.path(SCRIPT_DIR, "toolkit", "R")
if (!dir.exists(TOOLKIT_R)) stop("未找到 toolkit/R —— 函数库缺失")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a[1L])) b else a
for (f in sort(list.files(TOOLKIT_R, pattern = "\\.R$", full.names = TRUE))) source(f)

set.seed(597)   # 可复现(本流程为确定性估计,种子用于 trim&fill 等随机化步骤)

cat("Step 1/6: 读取 2x2 计数数据 + 计算效应量 + 随机效应合并...\n")
d <- read_table_smart(args$input)
stopifnot(all(c("tpos", "tneg", "cpos", "cneg") %in% names(d)))
fit_or <- ma_pairwise(d, measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                      slab = ~paste(author, year))
print(fit_or)
write.csv(ma_summary_row(fit_or, label = "BCG vaccine (OR, REML+KH)"),
          file.path(args$outdir, "597_pooled_summary.csv"), row.names = FALSE)

cat("Step 2/6: 异质性 —— 亚组(alloc)+ meta 回归(ablat)...\n")
sg <- ma_subgroup(fit_or, by = "alloc")
write.csv(sg$table, file.path(args$outdir, "597_subgroup_summary.csv"), row.names = FALSE)
mr <- ma_metareg(fit_or, mods = ~ablat, out = file.path(args$outdir, "597_metareg_bcg_ablat.pdf"))
write.csv(as.data.frame(mr$coef), file.path(args$outdir, "597_metareg_summary.csv"), row.names = FALSE)

cat("Step 3/6: 发表偏倚(Egger / rank / trim&fill / PET-PEESE)...\n")
pb <- ma_pubbias(fit_or, min_k = 10)
pb_row <- data.frame(
  k = pb$k, pooled = pb$pooled,
  egger.t = if (pb$egger$ok) pb$egger$stat else NA,
  egger.p = if (pb$egger$ok) pb$egger$p else NA,
  rank.tau = if (pb$rank$ok) pb$rank$tau else NA,
  rank.p = if (pb$rank$ok) pb$rank$p else NA,
  tf.k0 = if (pb$trimfill$ok) pb$trimfill$k0 else NA,
  tf.est = if (pb$trimfill$ok) pb$trimfill$est else NA,
  pp.model = if (pb$petpeese$ok) pb$petpeese$model else NA,
  pp.est = if (pb$petpeese$ok) pb$petpeese$corrected else NA)
write.csv(pb_row, file.path(args$outdir, "597_pubbias_summary.csv"), row.names = FALSE)
ma_funnel(fit_or, out = file.path(args$outdir, "597_funnel.pdf"))

cat("Step 4/6: 影响诊断(leave-one-out / 诊断量)...\n")
ma_influence(fit_or, out_prefix = file.path(args$outdir, "597_influence"))

cat("Step 5/6: 输出矢量图(森林图 / 漏斗图 / meta 回归图 → results/)...\n")
ma_forest(fit_or, out = file.path(args$outdir, "597_forest.pdf"))

cat("Step 6/6: 渲染 300dpi 展示图(→ assets/ 供 README 引用)...\n")
## ma_forest()/ma_funnel() 固定输出矢量 PDF(见 results/);这里用同一 metafor 对象
## 渲染 PNG 预览,数值与 PDF 完全一致。
png(file.path(ASSETS, "597_forest.png"), width = 9, height = max(5, 0.30 * fit_or$re$k + 3),
    units = "in", res = 300, type = "cairo")
metafor::forest(fit_or$re, atransf = exp, refline = 0, showweights = TRUE, addpred = TRUE,
                header = TRUE, xlab = "Odds Ratio (log scale)")
dev.off()
png(file.path(ASSETS, "597_funnel.png"), width = 7, height = 6, units = "in", res = 300,
    type = "cairo")
metafor::funnel(metafor::trimfill(fit_or$re), level = c(90, 95, 99),
                shade = c("white", "gray85", "gray70"), legend = TRUE,
                back = "white", refline = 0)
dev.off()

cat(sprintf("完成。合并 OR = %.3f (95%% CI %.3f–%.3f, 预测区间 %.3f–%.3f, p=%.3g), I² = %.1f%%, τ² = %.4f\n",
            exp(as.numeric(fit_or$re$b)), exp(fit_or$re$ci.lb), exp(fit_or$re$ci.ub),
            exp(fit_or$pred$pi.lb), exp(fit_or$pred$pi.ub), fit_or$re$pval,
            fit_or$re$I2, fit_or$re$tau2))
cat("输出: results/597_pooled_summary.csv · 597_subgroup_summary.csv ·",
    "597_metareg_summary.csv · 597_pubbias_summary.csv ·",
    "597_forest.pdf · 597_funnel.pdf · 597_metareg_bcg_ablat.pdf · 597_influence_* ·",
    "assets/597_forest.png · 597_funnel.png\n")

## 依赖快照(铁律 6):记录运行环境,便于复现
writeLines(capture.output(sessionInfo()), file.path(args$outdir, "597_sessionInfo.txt"))
