# =============================================================================
# 编号   : 557
# 脚本名 : sccomp 单细胞细胞组成差异检验(compositional differential abundance, DA)
# 分类   : 08_singlecell_spatial_trajectory
# 用途   : 检验各【细胞类型比例】是否随条件(condition)变化。sccomp 用
#          【sum-constrained beta-binomial】贝叶斯模型,在样本层级直接建模"计数
#          落入各细胞类型的比例",同时:① 尊重组成性(各类型比例和=1,非独立);
#          ② 联建【差异变异性】(formula_variability,某类型在某条件下更/欠分散);
#          ③ 离群稳健 + 经验贝叶斯收缩 → 小样本(每组 ~5)也稳健,假阳性受控。
#          输出每类型的 logit 效应 c_effect + 可信区间 + c_FDR(后验 H0 概率的 FDR)。
# ★诚实基线(内置对照,核心卖点):
#          (1) naive 比例 t-test —— 把各类型比例当【独立连续量】做 Welch t 检验。
#              这是组成性数据的【经典误用】:比例被迫和=1,当少数类型真变,其余
#              "稳定"类型的比例被动反向偏移 → naive 法把它们也判显著(假阳性级联)。
#              合成示例仅 2 类真变,naive 却判 ~7/8 显著 → 直观暴露陷阱。
#          (2) naive 卡方 —— 把"细胞总计数"汇总成 类型×条件 列联表做卡方;
#              忽略【样本层级变异】(把多个样本的细胞当独立观测=伪重复 pseudoreplication)
#              → 因 n=细胞数(上万)极大而过度显著。
#          (3) voomCLR(2026 对照法,CLR + limma)—— 中心对数比变换 + limma 经验贝叶斯
#              + 偏倚校正(applyBiasCorrection/topTableBC)。是组成性感知的【频率派】
#              对照;比 naive 大幅改善,但 CLR 的几何均值参照随真变类型漂移,仍会在
#              个别稳定类型上残留弱假阳性 → 与 sccomp 的样本级贝叶斯并列,体现方法谱。
#          → 三种基线与 sccomp 并排出图:展示"组成性建模 + 样本级层级"如何同时
#            压住伪重复(卡方)与组成性反弹(naive t)两类假阳性。
# ★工具接地(已最小实跑确认,见脚本注释):
#          sccomp 1.10.0 / voomCLR 0.99.41 / ggplot2 4.0.x。真实 API:
#          sccomp_estimate(.data, formula_composition=~condition, formula_variability=~1,
#            .sample=, .cell_group=, .abundance=, inference_method="pathfinder")
#          → sccomp_test(est, test_composition_above_logit_fold_change=) 返回 tibble:
#            cell_group/parameter/factor/c_lower/c_effect/c_upper/c_pH0/c_FDR/v_*。
#          sccomp_proportional_fold_change(est, ~condition, from=, to=) → 可解释比例倍数。
#          voomCLR(counts, design) → EList → lmFit → eBayes → applyBiasCorrection → topTableBC。
#          ★坑1:sccomp 把编译好的 Stan 模型缓存到【构建机器】的绝对路径
#            (C:\Users\biocbuild\...\.sccomp_models\),本机不可写 → 估计直接报
#            "cannot open the connection"。修复:估计前用 assignInNamespace() 把
#            包内常量 sccomp_stan_models_cache_dir 改指到可写目录(见 .fix_sccomp_cache)。
#          ★坑2:sccomp 自带 sccomp_boxplot()/plot_1D_intervals() 在 ggplot2 4.x 下
#            失效(内部用了已弃用的 S3 theme 构造 → "must be an <S7_object>")。本脚本
#            不依赖其绘图方法,全部从 sccomp_test 结果【自绘顶刊图】(也符合框架 §3:
#            不用包默认简陋图)。
# 依赖   : sccomp · voomCLR · limma · edgeR · ggplot2 · ggbeeswarm · dplyr · tidyr
#          (sccomp 需 cmdstanr + cmdstan 工具链;首次运行会编译 Stan 模型,稍慢)
# 运行   : Rscript 557_sccomp_composition_da.R                      # 合成示例(脚本内生成)
#          Rscript 557_sccomp_composition_da.R --input data/counts.csv --fdr 0.05 --lfc 0.2
# 输入   : long-format CSV(synthetic demo only,脚本自动生成到 example_data/):
#          列 sample,cell_group,count,condition —— 每行 = 某样本中某细胞类型的细胞计数。
#          换真实数据:从 Seurat/SCE 的 colData 用 table(sample, cell_type) 汇总即可
#          (见 README);condition 为样本级协变量,每个 sample 唯一对应一个 condition。
# =============================================================================

## ---- 定位共享顶刊主题库 _framework/theme_pub.R(同 519/558 样板)----------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
set.seed(42)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggbeeswarm)
}))

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input  = file.path(DDAT, "composition_counts.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  fdr    = 0.05,    # 显著阈值(sccomp c_FDR / 基线 adj.P)
  lfc    = 0.10,    # sccomp 检验的 logit 效应阈值(test_composition_above_logit_fold_change)
  cores  = 1))      # Stan 并行核数(小数据 1 即可)
for (k in c("fdr","lfc","cores")) args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 修复 sccomp 模型缓存到不可写构建路径的 bug(★坑1)----------------------
# sccomp 把编译好的 Stan 模型 .rds 缓存到包内常量 sccomp_stan_models_cache_dir,
# 该常量在【构建机器】上冻结为 C:\Users\biocbuild\...\.sccomp_models\<ver>,本机不存在
# 且不可写 → 估计阶段 gzfile/saveRDS 报 "cannot open the connection"。
# 修复:把该 namespace 常量改指到当前用户可写目录(~/.sccomp_models 或临时目录)。
.fix_sccomp_cache <- function() {
  ok <- requireNamespace("sccomp", quietly = TRUE)
  if (!ok) return(invisible(FALSE))
  cache <- tryCatch(normalizePath(file.path(path.expand("~"), ".sccomp_models",
                    as.character(utils::packageVersion("sccomp"))), mustWork = FALSE),
                    error = function(e) file.path(tempdir(), "sccomp_models"))
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  if (!file.access(cache, mode = 2) == 0) {           # 仍不可写 → 退回临时目录
    cache <- file.path(tempdir(), "sccomp_models"); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  }
  tryCatch(assignInNamespace("sccomp_stan_models_cache_dir", cache, ns = "sccomp"),
           error = function(e) invisible(FALSE))
  cache
}

# =============================================================================
# 0. 合成 long-format 计数(synthetic demo only)
#    设计:每条件 5 个样本(生物学重复),K=8 种细胞类型,每样本约 3000 个细胞。
#    真信号:CellType1 在 treat 中比例 ×1.9 富集;CellType5 在 treat 中 ×0.3 耗竭;
#    其余 6 类基线比例不变。但因比例【和=1】,这 6 类的【观测比例】仍会被动漂移
#    → 正是 naive t-test 误判它们的根源(组成性陷阱)。
# =============================================================================
gen_synthetic <- function(path, G_per = 5, K = 8, cells_per_sample = 3000) {
  samples <- c(sprintf("ctrl_%d", seq_len(G_per)), sprintf("treat_%d", seq_len(G_per)))
  cond    <- rep(c("ctrl", "treat"), each = G_per)
  base_p  <- c(0.30, 0.20, 0.15, 0.12, 0.08, 0.07, 0.05, 0.03)[seq_len(K)]
  truth   <- rep("unchanged", K); truth[1] <- "up (x1.9)"; truth[5] <- "down (x0.3)"
  mk <- function(cnd) {                                   # 一个样本的 K 维计数 ~ Multinomial
    p <- base_p
    if (cnd == "treat") { p[1] <- p[1] * 1.9; p[5] <- p[5] * 0.3 }
    p <- p / sum(p)
    # 加少量样本间过散(每样本 p 抖动)→ 让样本级模型有用武之地,贴近真实
    p <- p * exp(rnorm(K, 0, 0.08)); p <- p / sum(p)
    as.numeric(rmultinom(1, cells_per_sample, p))
  }
  mat <- vapply(seq_along(samples), function(i) mk(cond[i]), numeric(K))
  rownames(mat) <- paste0("CellType", seq_len(K)); colnames(mat) <- samples
  long <- as.data.frame(as.table(mat)); colnames(long) <- c("cell_group", "sample", "count")
  long$condition <- cond[match(long$sample, samples)]
  long$count <- as.integer(long$count)
  long <- long[, c("sample", "cell_group", "count", "condition")]
  attr(long, "truth") <- data.frame(cell_group = rownames(mat), truth = truth)
  write.csv(long, path, row.names = FALSE)
  # 旁置真值表(synthetic ground truth):便于复跑/读 CSV 时仍能做 sanity-check 与图注
  write.csv(data.frame(cell_group = rownames(mat), truth = truth),
            file.path(dirname(path), "composition_truth.csv"), row.names = FALSE)
  long
}

cat("Step 0: 准备 long-format 计数...\n")
if (!file.exists(args$input)) {
  long <- gen_synthetic(args$input)
  cat(sprintf("  [gen] 合成: %d 样本(ctrl×5/treat×5) × %d 类型 → %s\n",
              length(unique(long$sample)), length(unique(long$cell_group)), basename(args$input)))
} else {
  long <- read.csv(args$input, stringsAsFactors = FALSE)
  long$count <- as.integer(long$count)
  cat(sprintf("  [read] %s: %d 行, %d 样本, %d 类型\n", basename(args$input),
              nrow(long), length(unique(long$sample)), length(unique(long$cell_group))))
}
# 真值表(若有;真实数据无 truth 时为 NULL)——仅用于 sanity-check 与图注。
# 优先用合成时附在 long 上的属性;否则尝试读旁置的 composition_truth.csv(复跑/读 CSV 场景)。
truth_tab <- attr(long, "truth")
if (is.null(truth_tab)) {
  tpath <- file.path(dirname(args$input), "composition_truth.csv")
  if (file.exists(tpath)) truth_tab <- read.csv(tpath, stringsAsFactors = FALSE)
}
samp_cond <- unique(long[, c("sample", "condition")])
cond_levels <- c("ctrl", "treat")
if (!all(cond_levels %in% long$condition)) cond_levels <- sort(unique(long$condition))

# 宽矩阵(cell_group × sample)——voomCLR / 卡方 / 比例计算共用
wide <- long %>% select(sample, cell_group, count) %>%
  pivot_wider(names_from = sample, values_from = count) %>% as.data.frame()
rownames(wide) <- wide$cell_group; wide$cell_group <- NULL
mat <- as.matrix(wide)
samples <- colnames(mat)
cond_vec <- samp_cond$condition[match(samples, samp_cond$sample)]
prop <- sweep(mat, 2, colSums(mat), "/")                 # 各样本内比例(列和=1)

# =============================================================================
# 1. sccomp —— 贝叶斯 sum-constrained beta-binomial(主方法)
#    formula_composition=~condition:比例随 condition 变;
#    formula_variability=~condition:同时建模"离散度随条件变"(差异变异性)。
# =============================================================================
cat("Step 1: sccomp 估计(beta-binomial, pathfinder)...\n")
cache <- .fix_sccomp_cache()
cat(sprintf("  [cache fix] sccomp 模型缓存目录 → %s\n", cache))
suppressWarnings(suppressMessages(library(sccomp)))

est <- sccomp_estimate(
  long,
  formula_composition = ~ condition,    # 比例 ~ condition
  formula_variability = ~ condition,    # ★差异变异性:离散度也随 condition(sccomp 特色)
  .sample = sample, .cell_group = cell_group, .abundance = count,
  cores = args$cores, inference_method = "pathfinder", verbose = FALSE)

res <- sccomp_test(est, test_composition_above_logit_fold_change = args$lfc)
# 只取 condition 这个因子的行(剔除 Intercept;factor 列标了因子名,Intercept 行为 NA)
sccomp_cond <- res %>% filter(!is.na(factor) & factor == "condition") %>%
  transmute(cell_group, c_effect, c_lower, c_upper, c_pH0, c_FDR,
            v_effect, v_FDR, sig = c_FDR < args$fdr)
# 可解释的比例倍数(logit 效应 → proportion fold change)
pfc <- tryCatch(
  sccomp_proportional_fold_change(est, formula_composition = ~ condition,
                                  from = cond_levels[1], to = cond_levels[2]),
  error = function(e) NULL)
if (!is.null(pfc)) sccomp_cond <- sccomp_cond %>%
  left_join(pfc %>% select(cell_group, proportion_fold_change, statement), by = "cell_group")
write.csv(sccomp_cond, file.path(args$outdir, "sccomp_composition_DA.csv"), row.names = FALSE)
cat(sprintf("  sccomp 显著(c_FDR<%.2f): %d/%d 类型 → %s\n",
            args$fdr, sum(sccomp_cond$sig), nrow(sccomp_cond),
            paste(sccomp_cond$cell_group[sccomp_cond$sig], collapse = ", ")))
# sanity-check:真变类型应被检出,真不变类型不应被检出(验证管道有效)
if (!is.null(truth_tab)) {
  m <- sccomp_cond %>% left_join(truth_tab, by = "cell_group")
  tp <- sum(m$sig & m$truth != "unchanged"); fp <- sum(m$sig & m$truth == "unchanged")
  cat(sprintf("  [sanity] 真变=%d → sccomp 检出真变 TP=%d, 误检不变 FP=%d\n",
              sum(truth_tab$truth != "unchanged"), tp, fp))
}

# =============================================================================
# 2. ★诚实基线对照
# =============================================================================
cat("Step 2: 诚实基线(naive t / naive 卡方 / voomCLR)...\n")

## (2a) naive 比例 Welch t-test —— 把各类型比例当独立连续量 ----------------------
naive_t <- do.call(rbind, lapply(rownames(prop), function(g) {
  pc <- prop[g, cond_vec == cond_levels[1]]; pt <- prop[g, cond_vec == cond_levels[2]]
  tt <- tryCatch(t.test(pt, pc), error = function(e) list(statistic = NA, p.value = NA))
  data.frame(cell_group = g,
             diff = mean(pt) - mean(pc),
             t = unname(tt$statistic %||% NA), p = tt$p.value %||% NA)
}))
naive_t$adj.P <- p.adjust(naive_t$p, "BH"); naive_t$sig <- naive_t$adj.P < args$fdr
write.csv(naive_t, file.path(args$outdir, "baseline_naive_ttest.csv"), row.names = FALSE)

## (2b) naive 卡方 —— 类型×条件 细胞总计数列联表(忽略样本=伪重复)-------------
ct <- t(rowsum(t(mat), group = cond_vec))                # K × 2 条件 细胞总计数
naive_chisq <- do.call(rbind, lapply(rownames(ct), function(g) {
  a <- ct[g, cond_levels[2]]; b <- ct[g, cond_levels[1]]            # 该类型 treat/ctrl
  c2 <- sum(ct[, cond_levels[2]]) - a; d2 <- sum(ct[, cond_levels[1]]) - b  # 其余类型
  tab <- matrix(c(a, c2, b, d2), 2)
  ch <- suppressWarnings(chisq.test(tab))
  data.frame(cell_group = g, chisq = unname(ch$statistic), p = ch$p.value)
}))
naive_chisq$adj.P <- p.adjust(naive_chisq$p, "BH"); naive_chisq$sig <- naive_chisq$adj.P < args$fdr
write.csv(naive_chisq, file.path(args$outdir, "baseline_naive_chisq.csv"), row.names = FALSE)

## (2c) voomCLR(2026 对照)—— CLR + limma + 偏倚校正 ---------------------------
voom_tt <- NULL
try({
  suppressWarnings(suppressMessages({ library(voomCLR); library(limma) }))
  cond_f <- factor(cond_vec, levels = cond_levels)
  design <- model.matrix(~ cond_f)
  v   <- voomCLR(counts = mat, design = design)          # CLR + 方差权重(features=类型行)
  fit <- eBayes(lmFit(v, design))
  fit_bc <- tryCatch(applyBiasCorrection(fit), error = function(e) fit)  # 偏倚校正(失败则退回)
  tt  <- tryCatch(topTableBC(fit_bc, coef = 2, number = Inf, sort.by = "none"),
                  error = function(e) topTable(fit, coef = 2, number = Inf, sort.by = "none"))
  voom_tt <- data.frame(cell_group = rownames(tt), logFC = tt$logFC,
                        p = tt$P.Value, adj.P = tt$adj.P.Val)
  voom_tt$sig <- voom_tt$adj.P < args$fdr
  write.csv(voom_tt, file.path(args$outdir, "baseline_voomCLR.csv"), row.names = FALSE)
}, silent = TRUE)

# 汇总四法的"显著标记"(真值并入)----------------------------------------------
summ <- data.frame(cell_group = rownames(prop)) %>%
  left_join(sccomp_cond %>% transmute(cell_group, sccomp = sig), by = "cell_group") %>%
  left_join(naive_t     %>% transmute(cell_group, naive_t = sig), by = "cell_group") %>%
  left_join(naive_chisq %>% transmute(cell_group, naive_chisq = sig), by = "cell_group")
if (!is.null(voom_tt)) summ <- summ %>% left_join(voom_tt %>% transmute(cell_group, voomCLR = sig), by = "cell_group")
if (!is.null(truth_tab)) summ <- summ %>% left_join(truth_tab, by = "cell_group")
write.csv(summ, file.path(args$outdir, "method_comparison_significance.csv"), row.names = FALSE)
cat("  各法显著数: ",
    sprintf("sccomp=%d · naive_t=%d · naive_chisq=%d%s",
            sum(summ$sccomp, na.rm = TRUE), sum(summ$naive_t, na.rm = TRUE),
            sum(summ$naive_chisq, na.rm = TRUE),
            if (!is.null(voom_tt)) sprintf(" · voomCLR=%d", sum(summ$voomCLR, na.rm = TRUE)) else ""), "\n")
if (!is.null(truth_tab)) {
  fp_naive <- sum(summ$naive_t & summ$truth == "unchanged", na.rm = TRUE)
  fp_chi   <- sum(summ$naive_chisq & summ$truth == "unchanged", na.rm = TRUE)
  fp_scc   <- sum(summ$sccomp & summ$truth == "unchanged", na.rm = TRUE)
  cat(sprintf("  [对照核心] 不变类型上的假阳性 FP: naive_t=%d, naive_chisq=%d, sccomp=%d (真变仅 %d 类)\n",
              fp_naive, fp_chi, fp_scc, sum(truth_tab$truth != "unchanged")))
}

# =============================================================================
# 3. 出图(全部顶刊级;禁平凡条形图;每图独立成文件 PDF+PNG)
# =============================================================================
cat("Step 3: 出图...\n")
pal <- pal_pub(name = "npg")
W <- 7.2

## 图1:组成 boxplot + 每样本点(代替堆叠条形)--------------------------------
# 各类型在两条件下的【样本级比例】箱线 + 每个样本一个点(本设计为非配对两组,
# 故不画跨条件连线以免暗示伪配对);真变类型由 sccomp 显著性加 "*" 标到 facet。
prop_long <- as.data.frame(as.table(prop))
colnames(prop_long) <- c("cell_group", "sample", "proportion")
prop_long$condition <- factor(cond_vec[match(prop_long$sample, samples)], levels = cond_levels)
prop_long$cell_group <- factor(prop_long$cell_group, levels = rownames(prop))
# 给 sccomp 显著类型加 "*" 标记到 facet 标签
sig_set <- sccomp_cond$cell_group[sccomp_cond$sig]
lab_map <- setNames(ifelse(rownames(prop) %in% sig_set,
                           paste0(rownames(prop), " *"), rownames(prop)), rownames(prop))
prop_long$facet <- factor(lab_map[as.character(prop_long$cell_group)],
                          levels = lab_map[rownames(prop)])
p_box <- ggplot(prop_long, aes(condition, proportion)) +
  geom_boxplot(aes(fill = condition), width = 0.55, outlier.shape = NA, alpha = 0.85,
               colour = "grey20", linewidth = 0.4) +
  ggbeeswarm::geom_quasirandom(colour = "grey25", size = 1.4, alpha = 0.85, width = 0.16) +
  facet_wrap(~ facet, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(pal[4], pal[1]), name = "Condition") +
  labs(title = "Per-sample cell-type composition",
       subtitle = "Box = across samples; points = individual samples; * = significant by sccomp (beta-binomial)",
       x = NULL, y = "Proportion within sample") +
  theme_pub(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 0),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_box, file.path(ASSETS, "fig1_composition_boxplot"), width = W, height = 5.4)

## 图2:credible-effect lollipop —— sccomp c_effect ± 可信区间,色=方向 ---------
# y 轴按效应排序;棒=0→效应,点=c_effect,误差线=95% 可信区间;显著者实心描边。
le <- sccomp_cond %>%
  mutate(cell_group = factor(cell_group, levels = cell_group[order(c_effect)]),
         dir = ifelse(c_effect > 0, "Enriched in treat", "Depleted in treat"),
         dir = ifelse(sig, dir, "n.s."))
dir_cols <- c("Enriched in treat" = pal[1], "Depleted in treat" = pal[4], "n.s." = "grey70")
p_lol <- ggplot(le, aes(c_effect, cell_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = c_effect, yend = cell_group, colour = dir), linewidth = 0.7) +
  geom_errorbarh(aes(xmin = c_lower, xmax = c_upper, colour = dir), height = 0.22, linewidth = 0.5) +
  geom_point(aes(colour = dir, shape = sig, size = sig)) +
  scale_colour_manual(values = dir_cols, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.4, `TRUE` = 3.6), guide = "none") +
  labs(title = "sccomp credible composition effects",
       subtitle = sprintf("Bayesian beta-binomial logit effect ± 95%% credible interval (c_FDR<%.2g solid)", args$fdr),
       x = expression("Composition effect ("*log[e]~"fold-change, treat vs ctrl)"),
       y = NULL) +
  theme_pub(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_lol, file.path(ASSETS, "fig2_sccomp_credible_lollipop"), width = 7.0, height = 5.0)

## 图3:比例 raincloud —— 各类型样本级比例的 半提琴+箱+抖动点(按条件分色)-----
# raincloud = half-violin(密度)+ boxplot(摘要)+ jitter(原始样本点),代替条形。
rc <- prop_long
p_rain <- ggplot(rc, aes(cell_group, proportion, fill = condition, colour = condition)) +
  geom_violin(width = 0.85, alpha = 0.35, colour = NA, scale = "width",
              position = position_dodge(width = 0.8), trim = TRUE) +
  geom_boxplot(width = 0.16, alpha = 0.9, outlier.shape = NA, linewidth = 0.4,
               position = position_dodge(width = 0.8), colour = "grey20") +
  ggbeeswarm::geom_quasirandom(dodge.width = 0.8, size = 1.0, alpha = 0.8,
                               colour = "grey25", width = 0.10) +
  scale_fill_manual(values = c(pal[4], pal[1]), name = "Condition") +
  scale_colour_manual(values = c(pal[4], pal[1]), guide = "none") +
  labs(title = "Cell-type proportion distributions (raincloud)",
       subtitle = "Half-violin density + box summary + per-sample points, split by condition",
       x = NULL, y = "Proportion within sample") +
  theme_pub(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_rain, file.path(ASSETS, "fig3_proportion_raincloud"), width = W, height = 4.8)

## 图4:★诚实基线对照 —— 四法显著性 dot-matrix(谁判谁显著)--------------------
# 行=细胞类型(标真值),列=方法,点=是否显著(实心红/空心灰),
# 直观看到 naive 法在不变类型上"亮成一片"(假阳性级联),sccomp/voomCLR 干净。
meths <- c("sccomp", "voomCLR", "naive_t", "naive_chisq")
meths <- meths[meths %in% colnames(summ)]
heat <- summ %>% select(cell_group, all_of(meths)) %>%
  pivot_longer(-cell_group, names_to = "method", values_to = "sig")
heat$method <- factor(heat$method, levels = meths,
                      labels = c(sccomp = "sccomp\n(beta-binom)", voomCLR = "voomCLR\n(CLR+limma)",
                                 naive_t = "naive t-test\n(proportions)",
                                 naive_chisq = "naive chi-sq\n(pooled cells)")[meths])
heat$sig <- factor(ifelse(is.na(heat$sig), FALSE, heat$sig), levels = c(FALSE, TRUE))
# 行序:真变类型置顶(若有真值),并把真值并入 y 标签
if (!is.null(truth_tab)) {
  ord <- truth_tab$cell_group[order(truth_tab$truth == "unchanged", truth_tab$cell_group)]
  ylab <- setNames(ifelse(truth_tab$truth == "unchanged",
                          truth_tab$cell_group,
                          paste0(truth_tab$cell_group, "  (TRUE Δ)")), truth_tab$cell_group)
} else { ord <- sort(unique(heat$cell_group)); ylab <- setNames(ord, ord) }
heat$cell_group <- factor(heat$cell_group, levels = rev(ord))
p_cmp <- ggplot(heat, aes(method, cell_group)) +
  geom_point(aes(fill = sig, size = sig), shape = 21, colour = "grey30", stroke = 0.5) +
  scale_fill_manual(values = c(`FALSE` = "grey92", `TRUE` = pal[1]),
                    labels = c("n.s.", paste0("sig (FDR<", args$fdr, ")")), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 3.2, `TRUE` = 5.2), guide = "none") +
  scale_y_discrete(labels = ylab) +
  labs(title = "Honest baseline: who calls what significant?",
       subtitle = "Only 2 types truly change; naive t adds compositional FPs, naive chi-sq lights up all (pseudoreplication)",
       x = NULL, y = NULL) +
  theme_pub(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_cmp, file.path(ASSETS, "fig4_baseline_method_comparison"), width = 8.0, height = 4.8)

cat("完成。结果表见", normalizePath(args$outdir), ";图见 assets/\n")

# 依赖版本快照(铁律6)
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
