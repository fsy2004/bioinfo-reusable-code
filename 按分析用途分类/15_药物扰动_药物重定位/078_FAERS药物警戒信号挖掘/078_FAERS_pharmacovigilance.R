# =============================================================================
# 编号       : 078
# 脚本名     : FAERS 药物警戒信号挖掘 (turnkey + 顶刊图)
# 分类       : 15_药物扰动_药物重定位
# 用途       : 对药物-不良事件报告做四种不相称性分析(ROR/PRR/BCPNN-IC/EBGM),识别
#              安全信号,输出 ROR 森林图与信号热图。
# 方法/包    : base R(四算法)+ ggplot2;主题 theme_pub.R
# 结果图     : ROR_forest(信号森林图);Signal_heatmap(药×事件信号)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 078_FAERS_pharmacovigilance.R
# 运行(自己): Rscript 078_FAERS_pharmacovigilance.R --input data/cases.csv
# 输入规格 : 二选一 —— ① 原始报告行(列 case_id, drug, event);② 预计算四格计数(列 drug,event,n11,n10,n01,n00)。
# 整理日期 : 2026-06-23(turnkey 重构;保留四算法,补 theme_pub 森林/热图)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages(library(ggplot2)))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "drug_event_cases.csv"), outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/3: 读取数据 + 构建四格表...\n")
x <- read_table_smart(args$input)
if (all(c("n11", "n10", "n01", "n00") %in% names(x))) { ct <- x } else {
  stopifnot(all(c("drug", "event") %in% names(x))); cc <- intersect(c("case_id", "case", "id"), names(x))[1]
  if (is.na(cc)) { x$case_id <- seq_len(nrow(x)); cc <- "case_id" }
  x <- unique(x[, c(cc, "drug", "event")]); allc <- unique(x[[cc]])
  ct <- do.call(rbind, lapply(sort(unique(x$drug)), function(d) {
    cd <- unique(x[[cc]][x$drug == d])
    do.call(rbind, lapply(sort(unique(x$event)), function(e) {
      ce <- unique(x[[cc]][x$event == e])
      data.frame(drug = d, event = e, n11 = length(intersect(cd, ce)), n10 = length(setdiff(cd, ce)),
                 n01 = length(setdiff(ce, cd)), n00 = length(setdiff(allc, union(cd, ce)))) })) }))
}

cat("Step 2/3: 四算法不相称性...\n")
a <- ct$n11 + .5; b <- ct$n10 + .5; c <- ct$n01 + .5; d <- ct$n00 + .5; n <- a + b + c + d
ct$ROR <- (a * d) / (b * c); rl <- log(ct$ROR); rse <- sqrt(1/a + 1/b + 1/c + 1/d)
ct$ROR025 <- exp(rl - 1.96 * rse); ct$ROR975 <- exp(rl + 1.96 * rse)
ct$PRR <- (a / (a + b)) / (c / (c + d)); ct$chi2 <- n * (a * d - b * c)^2 / ((a + b) * (c + d) * (a + c) * (b + d))
ct$IC <- log2(a / (((a + b) * (a + c)) / n)); ct$IC025 <- ct$IC - 1.96 * sqrt(1 / pmax(a, 1)); ct$EBGM <- 2^ct$IC
ct$signal <- ct$n11 >= 3 & ct$ROR025 > 1 & ct$PRR >= 2 & ct$IC025 > 0   # 三准则共识信号
ct$pair <- paste(ct$drug, ct$event, sep = " — ")
write.csv(ct, file.path(args$outdir, "signals.csv"), row.names = FALSE)
cat("  共识信号:", sum(ct$signal), "/", nrow(ct), "个药-事件对\n")

cat("Step 3/3: ROR 森林图 + 信号热图...\n")
top <- ct[order(-ct$ROR), ]; top <- head(top[top$n11 >= 3, ], 15); top$pair <- factor(top$pair, levels = rev(top$pair))
p_for <- ggplot(top, aes(ROR, pair, colour = signal)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(aes(xmin = ROR025, xmax = ROR975), height = .3) + geom_point(size = 2.6) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#E64B35"), labels = c("ns", "signal"), name = NULL) +
  scale_x_log10() + labs(title = "Pharmacovigilance ROR (95% CI)", x = "ROR (log scale)", y = NULL) +
  theme_pub(base_size = 10, border = TRUE)
save_fig(p_for, file.path(ASSETS, "ROR_forest"), 7, 5.5); save_fig(p_for, file.path(args$outdir, "ROR_forest"), 7, 5.5)

p_hm <- ggplot(ct, aes(event, drug, fill = log2(ROR))) + geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = ifelse(signal, "★", "")), colour = "white", size = 4) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, name = "log2 ROR") +
  labs(title = "Signal heatmap (★ = signal)", x = NULL, y = NULL) +
  theme_pub(base_size = 11, border = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_hm, file.path(ASSETS, "Signal_heatmap"), 6.5, 4); save_fig(p_hm, file.path(args$outdir, "Signal_heatmap"), 6.5, 4)
cat("完成。药物警戒图/表见", normalizePath(args$outdir), "\n")
