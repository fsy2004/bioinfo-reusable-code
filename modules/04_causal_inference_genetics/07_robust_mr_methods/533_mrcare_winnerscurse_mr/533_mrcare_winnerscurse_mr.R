# =============================================================================
# 编号   : 533
# 脚本名 : Winner-curse-free 稳健 MR (CARE / RIVW 概念演示 + 诚实基线)
# 分类   : 09_mendelian_randomization
# 用途   : 演示「赢家诅咒(winner's curse)」如何让两样本 MR 因果效应被系统性虚高,
#          以及 CARE(JASA 2026, ChongWuLab)的 RIVW(re-randomized IVW)
#          内生化选择偏倚校正的思想。
#          - 真实工具:MRcare::mr_care() / RIVW()(装上时自动调用,见 README 安装命令);
#          - 降级路径(本机/无 MRcare 时):用已装的 TwoSampleMR 基础估计子(IVW/Egger/WM)
#            在合成数据上跑通,复现并校正赢家诅咒,出全部展示图。
# ★诚实基线 : 朴素 IVW/Egger 直接用「发现样本」筛过 p 的 SNP 效应(含赢家诅咒)→ SNP-暴露效应
#             被虚高(分母膨胀),致 IVW 因果估计被系统性「偏向零(衰减)」;对照「独立验证样本」
#             效应(无诅咒)= 真值参照。展示选择偏倚把因果估计拉偏了多少(本演示约 -17%)。
# 依赖   : TwoSampleMR(诚实基线/降级路径,已装) · ggplot2 · (可选 MRcare = 真实 CARE/RIVW)
# 运行   : Rscript 533_mrcare_winnerscurse_mr.R
#          Rscript 533_mrcare_winnerscurse_mr.R --n_snp 400 --causal 0.20 --p_sel 5e-5
# 输入   : 无需外部输入;脚本内合成两样本 MR summary(synthetic, demo only):
#          注入赢家诅咒——工具按「发现样本」p 值筛选,发现样本效应即被高估。
#          换真实数据时,装 MRcare 后走 mr_care()(见 README),本脚本主体演示其纠偏逻辑。
# =============================================================================

## ---- 定位共享框架 theme_pub.R(顶刊图主题/配色/save_fig)---------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  n_snp   = 600,        # 候选 SNP 数(synthetic)
  causal  = 0.15,       # 真因果效应 beta(参照真值)
  p_sel   = 5e-5,       # 工具选择阈值(CARE 默认 5e-5,故意比 5e-8 宽 → 诅咒更明显)
  outdir  = file.path(SCRIPT_DIR, "results")))
for (k in c("n_snp","causal","p_sel")) args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## 0. 真实工具 MRcare 入口(装上则用真 CARE/RIVW;未装则降级到诚实基线路径)
##    —— 不臆造 API:函数名/参数取自官方文档 https://chongwulab.github.io/MRcare/
##       (mr_care() 主函数 / RIVW() 估计子 / preprocess_gwas_data())。
## =============================================================================
MRCARE_OK <- FALSE
try({ suppressWarnings(suppressMessages(library(MRcare))); MRCARE_OK <- TRUE }, silent = TRUE)
cat(sprintf("[env] MRcare 可用: %s  (FALSE = 走诚实基线降级路径,见 README 顶部 🟡 安装命令)\n",
            MRCARE_OK))

## =============================================================================
## 1. 合成两样本 MR summary,显式注入赢家诅咒(synthetic, demo only)
##    机制:同一批 SNP 在两个【独立】GWAS 子样本里各测一次 SNP-暴露效应:
##      - 发现样本(discovery): 选择 + 估计都用它 → 选中的 SNP 效应被「赢家诅咒」高估;
##      - 验证样本(validation): 独立重测同一 SNP → 无诅咒,作真值参照(CARE 用再随机化模拟它)。
##    SNP-结局效应由真因果 causal × 真 beta 生成(真 beta 对两样本相同)。
## =============================================================================
cat("Step 1: 合成两样本 MR(注入 winner's curse)...\n")
n <- as.integer(args$n_snp)
# 关键设定:多数 SNP 为弱/近零真效应,少数为真工具。当用「发现样本 p」筛选时,
# 弱真效应 SNP 只有正好被噪声推大的那次才会过阈值 → 选中集的发现 beta 系统性高估
# (= 赢家诅咒)。这正是 CARE/RIVW 要内生化校正的偏倚来源。
frac_strong <- 0.08                                    # 仅 8% 为强真工具
is_strong   <- runif(n) < frac_strong
# 弱 SNP 的真效应略低于选择阈值附近,使其大多靠「噪声推大」才入选 → 诅咒最重
true_b   <- ifelse(is_strong, rnorm(n, 0, 0.050), rnorm(n, 0, 0.018))
se_exp   <- runif(n, 0.011, 0.015)                     # 发现样本测量 se
se_val   <- runif(n, 0.011, 0.015)                     # 验证样本测量 se(独立)
b_disc   <- true_b + rnorm(n, 0, 1) * se_exp           # 发现样本观测 beta(含噪)
b_val    <- true_b + rnorm(n, 0, 1) * se_val           # 验证样本观测 beta(独立噪声)
p_disc   <- 2 * pnorm(-abs(b_disc / se_exp))           # 发现样本 p(用于选工具)
# 结局:由真 beta 决定(真因果 causal),与「观测到哪个样本」无关
se_out   <- runif(n, 0.010, 0.018)
b_out    <- args$causal * true_b + rnorm(n, 0, 1) * se_out

snp <- sprintf("rs%07d", sample(1e6:9e6, n))
dat_all <- data.frame(SNP = snp,
  b_disc = b_disc, se_disc = se_exp, p_disc = p_disc,   # 发现样本(被诅咒)
  b_val  = b_val,  se_val  = se_val,                    # 验证样本(无诅咒,真值参照)
  b_out  = b_out,  se_out  = se_out, true_b = true_b)
write.csv(dat_all, file.path(DDAT, "synthetic_two_sample_mr.csv"), row.names = FALSE)

## 按发现样本 p 值筛工具(这一步就是赢家诅咒的来源)----------------------------
sel <- dat_all[dat_all$p_disc < args$p_sel, ]
cat(sprintf("  候选 %d SNP → p_disc<%.0e 选中 %d 个工具\n", n, args$p_sel, nrow(sel)))
if (nrow(sel) < 5) stop("选中工具过少;请增大 --n_snp 或放宽 --p_sel")

## 量化诅咒幅度:选中 SNP 的发现 beta vs 验证 beta(同一真值,差额=诅咒膨胀)----
infl <- mean(abs(sel$b_disc)) / mean(abs(sel$b_val))
cat(sprintf("  ★ 选中工具效应膨胀比 |b_disc|/|b_val| = %.2f (>1 即赢家诅咒)\n", infl))

## =============================================================================
## 2. 诚实基线 MR:朴素(被诅咒) vs 无诅咒参照
##    估计子来自 TwoSampleMR 底层函数(mr_ivw / mr_egger_regression /
##    mr_weighted_median),输入原始 beta/se 向量,真实 API,非臆造。
## =============================================================================
cat("Step 2: 诚实基线 — 朴素 IVW/Egger/WM(含诅咒) vs 验证样本参照(无诅咒)...\n")
suppressWarnings(suppressMessages(library(TwoSampleMR)))
PAR <- default_parameters()

run_panel <- function(b_exp, se_exp, b_out, se_out, label) {
  iv <- mr_ivw(b_exp, b_out, se_exp, se_out, PAR)
  eg <- mr_egger_regression(b_exp, b_out, se_exp, se_out, PAR)
  wm <- mr_weighted_median(b_exp, b_out, se_exp, se_out, PAR)
  rbind(
    data.frame(scenario = label, method = "IVW",             b = iv$b, se = iv$se, pval = iv$pval),
    data.frame(scenario = label, method = "MR-Egger",        b = eg$b, se = eg$se, pval = eg$pval),
    data.frame(scenario = label, method = "Weighted median", b = wm$b, se = wm$se, pval = wm$pval))
}

## (A) 朴素/被诅咒:选择 + 估计都用发现样本效应(标准两样本 MR 的常见错误)
res_naive <- run_panel(sel$b_disc, sel$se_disc, sel$b_out, sel$se_out, "Naive (winner's curse)")
## (B) 无诅咒参照:对【同一批选中工具】改用独立验证样本效应(CARE 再随机化所逼近的目标)
res_corr  <- run_panel(sel$b_val,  sel$se_val,  sel$b_out, sel$se_out, "Curse-corrected (CARE-style)")

## (C) 真实 CARE/RIVW(仅当 MRcare 装上;否则跳过,不臆造结果)------------------
res_care <- NULL
if (MRCARE_OK) {
  try({
    # 真实 API(官方文档): RIVW(...) 再随机化 IVW。装上时按其 vignette 入参调用。
    # 此处仅在包真实可用时执行,避免编造返回值。
    rv <- MRcare::RIVW(
      b_exp = sel$b_disc, se_exp = sel$se_disc,
      b_out = sel$b_out,  se_out = sel$se_out,
      p_threshold = args$p_sel)
    res_care <- data.frame(scenario = "MRcare RIVW (real)", method = "RIVW",
                           b = rv$estimate, se = rv$standard_error, pval = rv$p_value)
  }, silent = TRUE)
}

res <- rbind(res_naive, res_corr, if (!is.null(res_care)) res_care)
res$lci <- res$b - 1.96 * res$se
res$uci <- res$b + 1.96 * res$se
res$true_b <- args$causal
write.csv(res, file.path(args$outdir, "MR_results_winnerscurse.csv"), row.names = FALSE)

ivw_naive <- res_naive$b[res_naive$method == "IVW"]
ivw_corr  <- res_corr$b[res_corr$method == "IVW"]
cat(sprintf("  ★ 诚实基线实测: 真因果=%.3f | 朴素 IVW=%.3f (偏倚 %+.1f%%) | 校正后 IVW=%.3f (偏倚 %+.1f%%)\n",
  args$causal, ivw_naive, 100 * (ivw_naive - args$causal) / args$causal,
  ivw_corr,  100 * (ivw_corr  - args$causal) / args$causal))

## =============================================================================
## 3. 顶刊级展示图(全部非条形图;每图独立成文件,save_fig 出 PDF+PNG)
## =============================================================================
cat("Step 3: 出图(lollipop / forest / 校正斜率散点 / 诅咒膨胀 dumbbell)...\n")
PAL <- pal_pub(3, "npg")
names(PAL) <- c("Naive (winner's curse)", "Curse-corrected (CARE-style)", "MRcare RIVW (real)")

## (图1) 各方法因果估计 lollipop —— 朴素 vs 校正,横虚线=真值 -------------------
res$lab <- paste0(res$method, "\n", sub(" .*", "", res$scenario))
p1 <- ggplot(res, aes(x = reorder(interaction(method, scenario), b), y = b, colour = scenario)) +
  geom_hline(yintercept = args$causal, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.7, y = args$causal, label = "True causal effect",
           hjust = 0, vjust = -0.6, size = 3, colour = "grey30") +
  geom_linerange(aes(ymin = 0, ymax = b), linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 3.2) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.7) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +
  scale_colour_manual(values = PAL, name = NULL) +
  coord_flip() +
  labs(title = "Causal estimates: naive vs winner's-curse-corrected",
       subtitle = "Curse inflates SNP-exposure -> naive estimate biased toward the null",
       x = NULL, y = "Causal effect estimate (95% CI)") +
  theme_pub(base_size = 11, legend = "bottom")
save_fig(p1, file.path(ASSETS, "fig1_estimates_lollipop"), width = 7, height = 5)

## (图2) CARE-style vs IVW/Egger forest —— 经典森林图,真值参照线 ----------------
res$ord <- factor(paste0(res$scenario, " · ", res$method),
                  levels = rev(unique(paste0(res$scenario, " · ", res$method))))
p2 <- ggplot(res, aes(x = b, y = ord, colour = scenario)) +
  geom_vline(xintercept = args$causal, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
  geom_errorbar(aes(xmin = lci, xmax = uci), orientation = "y", width = 0.22, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = PAL, name = NULL) +
  labs(title = "Forest plot: CARE-style correction vs IVW / MR-Egger",
       subtitle = "Dashed line = true causal effect; naive estimates sit short of it (attenuated)",
       x = "Causal effect estimate (95% CI)", y = NULL) +
  theme_pub(base_size = 11, legend = "bottom")
save_fig(p2, file.path(ASSETS, "fig2_forest_care_vs_classic"), width = 7.2, height = 5)

## (图3) SNP-暴露 vs SNP-结局 散点 + 校正斜率 ----------------------------------
##   两组点:被诅咒(发现 beta,膨胀) vs 校正(验证 beta);叠各自 IVW 斜率(过原点)。
sc <- rbind(
  data.frame(bx = sel$b_disc, by = sel$b_out, scenario = "Naive (winner's curse)"),
  data.frame(bx = sel$b_val,  by = sel$b_out, scenario = "Curse-corrected (CARE-style)"))
slopes <- data.frame(
  scenario = c("Naive (winner's curse)", "Curse-corrected (CARE-style)"),
  slope = c(ivw_naive, ivw_corr))
p3 <- ggplot(sc, aes(bx, by, colour = scenario)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point(alpha = 0.5, size = 1.6) +
  geom_abline(data = slopes, aes(slope = slope, intercept = 0, colour = scenario),
              linewidth = 1, show.legend = FALSE) +
  geom_abline(slope = args$causal, intercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = PAL, name = NULL) +
  labs(title = "SNP-exposure vs SNP-outcome with fitted IVW slope",
       subtitle = "Curse inflates SNP-exposure effects (x-axis) -> shallower, biased slope",
       x = "SNP-exposure effect", y = "SNP-outcome effect") +
  theme_pub(base_size = 11, legend = "bottom")
save_fig(p3, file.path(ASSETS, "fig3_scatter_corrected_slope"), width = 6.8, height = 5.4)

## (图4) 赢家诅咒膨胀 dumbbell —— 每个工具 验证 beta → 发现 beta 的拉伸 ----------
dd <- sel[order(abs(sel$b_val)), ]
dd$idx <- seq_len(nrow(dd))
top <- tail(dd, min(30, nrow(dd)))   # 取效应最大的 30 个工具,避免过密
p4 <- ggplot(top) +
  geom_segment(aes(x = abs(b_val), xend = abs(b_disc), y = idx, yend = idx),
               colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = abs(b_val),  y = idx, colour = "Validation (no curse)"),  size = 2) +
  geom_point(aes(x = abs(b_disc), y = idx, colour = "Discovery (curse)"), size = 2) +
  scale_colour_manual(values = c("Validation (no curse)" = unname(PAL[2]),
                                 "Discovery (curse)" = unname(PAL[1])), name = NULL) +
  labs(title = "Winner's curse per instrument (|effect| inflation)",
       subtitle = sprintf("Discovery-sample effects are inflated x%.2f on average", infl),
       x = "|SNP-exposure effect|", y = "Instrument (sorted)") +
  theme_pub(base_size = 11, legend = "bottom") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_fig(p4, file.path(ASSETS, "fig4_winnerscurse_dumbbell"), width = 6.6, height = 5.6)

cat("完成。结果表见", normalizePath(args$outdir), ";展示图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖快照(铁律6)
