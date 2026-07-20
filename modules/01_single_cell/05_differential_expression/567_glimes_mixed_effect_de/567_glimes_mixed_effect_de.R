# =============================================================================
# 567 · GLIMES —— 单细胞广义线性混合效应差异表达 (donor random effect)
# -----------------------------------------------------------------------------
# 上游方法 : Wu CH, Zhou X, Chen M. Exploring and mitigating shortcomings in
#            single-cell differential expression analysis with a new statistical
#            paradigm. Genome Biol 2025;26(1):58.
#            PMID 40098192 · doi:10.1186/s13059-025-03525-6   (已核实,见 README)
# 上游仓库 : https://github.com/C-HW/GLIMES  (仅 GitHub,不在 CRAN/Bioconductor)
# 上游 API : 已实读克隆到本地的上游源码 upstream-sources/567_GLIMES/R/DE_methods.R
#            与 NAMESPACE(导出 poisson_glmm_DE:45 / binomial_glmm_DE:144 /
#            identifyDEGs:255 / simple_mean_DE:12),函数名、参数名与默认值以该源码为准。
#
# 本模块做什么:
#   在「多供体 (donor) × 两组条件」的单细胞计数上比较三条 DE 路线,
#   核心问题是 **供体层面的伪重复 (pseudoreplication)**:
#     A. naive_ttest  —— 把细胞当独立重复直接做 t 检验 → 一类错误膨胀
#                        (与上游 simple_mean_DE 同属"细胞级朴素对照",但**不是复刻**:
#                         上游的 p 值用固定自由度 2*n_cells-1 手算,本模块用 Welch
#                         t.test 的 p 值。详见 README「API 与诚实边界」。)
#     B. pseudobulk   —— 每供体求和成 pseudobulk 再做供体级 t 检验(标准对照)
#     C. poisson_glmm —— GLIMES 的模型:原始 UMI 计数 + 供体随机截距,
#                        MASS::glmmPQL(family = poisson)
#   合成数据自带 ground truth,因此可以直接量出三者的 ROC-AUC、空基因假阳率
#   与 p 值均匀性。B/C 两条基线本机零依赖即可跑(MASS/nlme 随 R 分发)。
#
# ★ 诚实声明:本机未安装 GLIMES 包(不在 CRAN/Bioc,需 devtools::install_github)。
#   - 加 --use-glimes 时脚本会守卫式尝试 library(GLIMES) 并调用**官方**
#     poisson_glmm_DE()/binomial_glmm_DE(),失败则打印真实安装命令后跳过。
#   - 默认路径下的 poisson_glmm 是**按上游源码逐行复刻的等价拟合**
#     (同样是 MASS::glmmPQL(count ~ comparison, random = list(replicates = ~1),
#      family = poisson) + BH 校正),已在 README 中标注为「复刻,非官方包」。
#     未安装官方包时的官方结果一栏 = 未验证,不做等价性断言。
#
# 用法:
#   Rscript 567_glimes_mixed_effect_de.R
#   Rscript 567_glimes_mixed_effect_de.R --counts data/my_counts.csv \
#           --meta data/my_meta.csv --group condition --donor donor --outdir results/run1
#   Rscript 567_glimes_mixed_effect_de.R --use-glimes      # 若已装官方包
# =============================================================================

suppressWarnings(suppressMessages({
  library(MASS)      # glmmPQL(随 R 分发)
  library(nlme)      # glmmPQL 后端
  library(ggplot2)
}))

set.seed(2026)   # 固定随机种子

# ---- 定位脚本目录 + 载入全库统一顶刊主题 ------------------------------------
.script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", a[m[1]]))))
  getwd()
}
ROOT <- .script_dir()
source(file.path(ROOT, "..", "..", "..", "_framework", "theme_pub.R"))

# ---- 参数区(默认全部指向 example_data/,支持 --key value 覆盖)--------------
args <- bio_args(list(
  counts     = file.path(ROOT, "example_data", "counts.csv"),
  meta       = file.path(ROOT, "example_data", "cell_meta.csv"),
  truth      = file.path(ROOT, "example_data", "truth.csv"),  # 可选:合成数据的金标准
  outdir     = file.path(ROOT, "results"),
  assets     = file.path(ROOT, "assets"),
  group      = "condition",     # 两水平的比较变量 (对应上游 comparison)
  donor      = "donor",         # 供体/重复变量   (对应上游 replicates)
  freq_expressed = "0.05",      # 基因检出率下限 (对应上游 freq_expressed)
  `use-glimes` = FALSE          # 守卫式调用官方 GLIMES 包
))
FREQ <- as.numeric(args$freq_expressed)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$assets, recursive = TRUE, showWarnings = FALSE)

read_csv_c <- function(p) {
  if (!file.exists(p)) stop(sprintf("输入文件不存在: %s", p))
  utils::read.csv(p, header = TRUE, comment.char = "#", check.names = FALSE,
                  stringsAsFactors = FALSE)
}

# =============================== Step 1 读数据 ================================
cat("Step 1 · 读取计数矩阵与细胞元数据\n")
cnt_df <- read_csv_c(args$counts)
meta   <- read_csv_c(args$meta)

gene_ids <- as.character(cnt_df[[1]])
counts   <- as.matrix(cnt_df[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- gene_ids

stopifnot(all(c(args$group, args$donor) %in% colnames(meta)))
cell_col <- colnames(meta)[1]
meta <- meta[match(colnames(counts), meta[[cell_col]]), , drop = FALSE]
if (anyNA(meta[[cell_col]])) stop("counts 的细胞列名与 meta 第一列对不上")

grp   <- factor(meta[[args$group]])
donor <- factor(meta[[args$donor]])
if (nlevels(grp) != 2) stop(sprintf("--group '%s' 必须恰好两个水平,当前 %d 个",
                                    args$group, nlevels(grp)))
cat(sprintf("  基因 %d · 细胞 %d · 供体 %d · 组别 %s\n",
            nrow(counts), ncol(counts), nlevels(donor),
            paste(levels(grp), collapse = " vs ")))

# 供体-组别对应表:检查每组供体数(伪重复问题的真实自由度)
donor_tab <- unique(data.frame(donor = as.character(donor), group = as.character(grp)))
n_donor_per_grp <- table(donor_tab$group)
cat("  每组供体数 (真实统计自由度):",
    paste(sprintf("%s=%d", names(n_donor_per_grp), as.integer(n_donor_per_grp)),
          collapse = ", "), "\n")

# 检出率:仅方法 C (GLMM) 会据此跳过低表达基因(与上游 freq_expressed 语义一致:
# 检出率 <= 阈值 的基因不拟合,状态标为 lowly expressed / zero mean)。
# 方法 A/B 两条基线不做此过滤,以保持与"朴素做法"的对照真实。
det_rate <- rowMeans(counts > 0)
cat(sprintf("  检出率 > %.2f 的基因(仅方法 C 会拟合这些): %d / %d\n",
            FREQ, sum(det_rate > FREQ), nrow(counts)))

# =========================== Step 2 方法 A 朴素细胞级 t 检验 ==================
# 对原始计数按细胞分组做两样本 Welch t 检验:"把细胞当独立重复"的朴素对照。
# ★ 不是 GLIMES::simple_mean_DE 的复刻 —— 上游那个函数取 t.test 的统计量后,
#   用 2*pt(-|t|, df = length(cellgroup1)+length(cellgroup2)-1) 手算 p 值
#   (两个入参是全长逻辑向量,故 df = 2*n_cells-1),与 Welch 自由度不同。
#   本模块只需要一条"朴素细胞级基线",直接用 t.test 的 p 值,不冒充上游实现。
cat("Step 2 · 方法 A:naive cell-level t-test(把细胞当重复 = 伪重复)\n")
g1 <- grp == levels(grp)[1]; g2 <- grp == levels(grp)[2]
naive <- data.frame(gene = gene_ids, t = NA_real_, pval = NA_real_,
                    log2FC = NA_real_, stringsAsFactors = FALSE)
for (i in seq_len(nrow(counts))) {
  x <- counts[i, g1]; y <- counts[i, g2]
  tt <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
  if (!is.null(tt)) { naive$t[i] <- unname(tt$statistic); naive$pval[i] <- tt$p.value }
  naive$log2FC[i] <- log2((mean(y) + 1e-8) / (mean(x) + 1e-8))
}
naive$BH <- stats::p.adjust(naive$pval, method = "BH")

# =========================== Step 3 方法 B pseudobulk ========================
# 每供体把细胞计数加和 → log-CPM → 供体级两样本 t 检验(n = 每组供体数)。
cat("Step 3 · 方法 B:pseudobulk(供体求和 → 供体级 t 检验)\n")
pb <- t(rowsum(t(counts), group = as.character(donor)))         # gene × donor
lib <- colSums(pb)
pb_cpm <- log2(t(t(pb) / lib) * 1e6 + 1)
pb_grp <- factor(donor_tab$group[match(colnames(pb), donor_tab$donor)])
pbA <- pb_grp == levels(grp)[1]; pbB <- pb_grp == levels(grp)[2]
pseudobulk <- data.frame(gene = gene_ids, pval = NA_real_, log2FC = NA_real_,
                         stringsAsFactors = FALSE)
for (i in seq_len(nrow(pb_cpm))) {
  x <- pb_cpm[i, pbA]; y <- pb_cpm[i, pbB]
  tt <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
  if (!is.null(tt)) pseudobulk$pval[i] <- tt$p.value
  pseudobulk$log2FC[i] <- mean(y) - mean(x)
}
pseudobulk$BH <- stats::p.adjust(pseudobulk$pval, method = "BH")

# =========================== Step 4 方法 C Poisson-GLMM ======================
# GLIMES 的核心模型。拟合调用逐行复刻自上游源码 R/DE_methods.R:
#   MASS::glmmPQL(count ~ comparison, random = list(replicates = ~1),
#                 family = poisson, data = countdf)
# 输出列名也对齐上游(genes / mu / beta_comparison / log2FC / sigma_square / status /
# pval / BH / log2mean / log2meandiff,见上游 R/DE_methods.R:74-76)。
# ★ 这是复刻,不是官方包。
cat("Step 4 · 方法 C:Poisson-GLMM + 供体随机截距 (GLIMES 模型,glmmPQL)\n")
glmm_fit <- function(counts, comparison, replicates, freq_expressed = 0.05) {
  countdf <- data.frame(comparison = as.factor(comparison),
                        replicates = as.factor(replicates))
  random_effects_list <- list(replicates = ~1)
  df <- data.frame(genes = rownames(counts), mu = NA_real_, beta_comparison = NA_real_,
                   log2FC = NA_real_, sigma_square = NA_real_, status = "done",
                   pval = NA_real_, BH = NA_real_, log2mean = NA_real_,
                   log2meandiff = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(counts))) {
    if (i %% max(1, round(nrow(counts) / 5)) == 0)
      cat(sprintf("    进度 %d/%d 基因\n", i, nrow(counts)))
    countdf$count <- as.numeric(round(pmax(counts[i, ], 0)))
    gm_mean <- stats::aggregate(count ~ comparison, data = countdf, FUN = mean, na.rm = TRUE)
    gm_mean <- gm_mean[order(gm_mean$comparison), ]
    m1 <- gm_mean[1, 2]; m2 <- gm_mean[2, 2]
    df$log2mean[i]     <- log2(m1 * m2) / 2
    df$log2meandiff[i] <- log2(abs(m1 - m2))
    dr <- mean(countdf$count != 0, na.rm = TRUE)
    if (dr <= freq_expressed) {
      df$status[i] <- if (dr == 0) "zero mean" else "lowly expressed"; next
    }
    gm <- tryCatch(suppressWarnings(summary(MASS::glmmPQL(
            count ~ comparison, random = random_effects_list,
            family = stats::poisson, data = countdf, verbose = FALSE))),
          error = function(e) NULL)
    if (is.null(gm)) { df$status[i] <- "not converge"; next }
    df$pval[i]            <- gm$tTable[2, "p-value"]
    df$sigma_square[i]    <- gm$sigma^2
    df$mu[i]              <- gm$coefficients$fixed[1]
    df$beta_comparison[i] <- gm$coefficients$fixed[2]
  }
  df$log2FC <- log2(exp(df$beta_comparison))
  df$BH     <- stats::p.adjust(df$pval, method = "BH")
  df
}
glmm <- glmm_fit(counts, grp, donor, freq_expressed = FREQ)
cat("  拟合状态:", paste(sprintf("%s=%d", names(table(glmm$status)),
                                  as.integer(table(glmm$status))), collapse = ", "), "\n")

# ------- identifyDEGs:上游的"新判定标准"(复刻自源码,含默认阈值) ----------
identify_degs <- function(adj_pval, log2FC, log2mean = NA, log2meandiff = -Inf,
                          pvalcutoff = 0.05, log2FCcutoff = log2(1.5),
                          log2meancutoff = -2.25, log2meandiffcutoff = -1,
                          newcriteria = TRUE) {
  d <- if (newcriteria)
    adj_pval < pvalcutoff & abs(log2FC) > log2FCcutoff &
      (log2mean > log2meancutoff | log2meandiff > log2meandiffcutoff)
  else adj_pval < pvalcutoff & abs(log2FC) > log2FCcutoff
  ifelse(is.na(adj_pval), NA, d)
}
glmm$DEG_new <- identify_degs(glmm$BH, glmm$log2FC, glmm$log2mean, glmm$log2meandiff,
                              newcriteria = TRUE)
glmm$DEG_old <- identify_degs(glmm$BH, glmm$log2FC, newcriteria = FALSE)

# =========================== Step 5 官方 GLIMES 包(守卫式)===================
official <- NULL
if (isTRUE(args$`use-glimes`) || identical(args$`use-glimes`, "TRUE")) {
  cat("Step 5 · 尝试调用官方 GLIMES 包\n")
  # 上游 poisson_glmm_DE 内部直接调用 SummarizedExperiment::colData(),故一并检查
  ok <- requireNamespace("GLIMES", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE) &&
        requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("S4Vectors", quietly = TRUE)
  if (!ok) {
    cat("  [跳过] 未检出 GLIMES / SingleCellExperiment / SummarizedExperiment。真实安装命令:\n",
        '    install.packages("devtools")\n',
        '    devtools::install_github("C-HW/GLIMES")\n',
        "    # 依赖 (Bioconductor): SummarizedExperiment, edgeR, MAST\n",
        "  官方函数签名(实读上游源码,未在本机验证运行):\n",
        "    poisson_glmm_DE(sce, comparison, replicates, exp_batch = NULL,",
        " other_fixed = NULL, freq_expressed = 0.05)\n",
        "    binomial_glmm_DE(sce, comparison, replicates, exp_batch = NULL,",
        " other_fixed = NULL, freq_expressed = 0.05)\n", sep = "")
  } else {
    # 上游按列名取 colData(sce)[, comparison],并按 sce@assays@data$counts[i, ] 取计数,
    # 故 assay 必须命名为 "counts",colData 行序必须与 counts 列序一致(上面已 match 对齐)。
    cd <- S4Vectors::DataFrame(meta)
    rownames(cd) <- colnames(counts)
    sce <- SingleCellExperiment::SingleCellExperiment(
             assays = list(counts = counts), colData = cd)
    official <- GLIMES::poisson_glmm_DE(sce, comparison = args$group,
                                        replicates = args$donor,
                                        freq_expressed = FREQ)
    utils::write.csv(official, file.path(args$outdir, "glimes_official_poisson.csv"),
                     row.names = FALSE)
    cat("  官方结果已写出 glimes_official_poisson.csv\n")
  }
} else {
  cat("Step 5 · 跳过官方 GLIMES 包(加 --use-glimes 启用;本机未安装)\n")
}

# =========================== Step 6 汇总 + 评估 ==============================
cat("Step 6 · 汇总三法结果并评估\n")
res <- data.frame(
  gene            = gene_ids,
  naive_pval      = naive$pval,      naive_BH      = naive$BH,      naive_log2FC      = naive$log2FC,
  pseudobulk_pval = pseudobulk$pval, pseudobulk_BH = pseudobulk$BH, pseudobulk_log2FC = pseudobulk$log2FC,
  glmm_pval       = glmm$pval,       glmm_BH       = glmm$BH,       glmm_log2FC       = glmm$log2FC,
  glmm_sigma_sq   = glmm$sigma_square, glmm_status  = glmm$status,
  glmm_DEG_new    = glmm$DEG_new,    glmm_DEG_old  = glmm$DEG_old,
  detection_rate  = det_rate, stringsAsFactors = FALSE)

has_truth <- file.exists(args$truth)
if (has_truth) {
  tr <- read_csv_c(args$truth)
  res$gene_class <- tr$gene_class[match(res$gene, tr$gene)]
  res$is_de      <- tr$is_de[match(res$gene, tr$gene)]
}
utils::write.csv(res, file.path(args$outdir, "de_results_all_methods.csv"), row.names = FALSE)

METHODS <- c(naive = "naive_pval", pseudobulk = "pseudobulk_pval", glmm = "glmm_pval")
MLAB <- c(naive = "Naive cell-level t-test", pseudobulk = "Pseudobulk t-test",
          glmm = "Poisson-GLMM (donor RE)")

# Mann-Whitney 式 AUC(不引外部包)
auc_of <- function(score, label) {
  ok <- !is.na(score) & !is.na(label)
  s <- score[ok]; y <- label[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
roc_of <- function(score, label) {
  ok <- !is.na(score) & !is.na(label)
  s <- score[ok]; y <- label[ok]
  o <- order(s, decreasing = TRUE); y <- y[o]
  data.frame(fpr = cumsum(y == 0) / sum(y == 0), tpr = cumsum(y == 1) / sum(y == 1))
}

metrics <- do.call(rbind, lapply(names(METHODS), function(m) {
  p <- res[[METHODS[m]]]; bh <- res[[sub("_pval", "_BH", METHODS[m])]]
  out <- data.frame(method = MLAB[m], n_tested = sum(!is.na(p)),
                    n_sig_BH05 = sum(bh < 0.05, na.rm = TRUE))
  if (has_truth) {
    nullg <- res$is_de == 0
    out$AUC <- auc_of(-log10(p), res$is_de)
    out$FDP_at_BH05 <- {
      sig <- which(bh < 0.05); if (!length(sig)) NA_real_ else mean(res$is_de[sig] == 0)
    }
    out$null_FPR_raw_p05 <- mean(p[nullg] < 0.05, na.rm = TRUE)   # 名义 5%,越接近越好
    out$null_pval_KS_p <- tryCatch(
      suppressWarnings(stats::ks.test(p[nullg & !is.na(p)], "punif")$p.value),
      error = function(e) NA_real_)
    out$power_true_DE_BH05 <- mean(bh[res$is_de == 1] < 0.05, na.rm = TRUE)
  }
  out
}))
utils::write.csv(metrics, file.path(args$outdir, "metrics_summary.csv"), row.names = FALSE)
cat("\n---- 评估汇总 ----\n"); print(metrics, row.names = FALSE); cat("\n")

# =========================== Step 7 出图(顶刊风格,无条形图)=================
cat("Step 7 · 出图\n")
PAL <- setNames(pal_pub(3, "nejm"), MLAB)
FIG <- function(n) file.path(args$assets, n)

# 图 1:空基因 p 值 ECDF vs 均匀分布 —— 一类错误控制(折线)
if (has_truth) {
  ecdf_df <- do.call(rbind, lapply(names(METHODS), function(m) {
    p <- res[[METHODS[m]]][res$is_de == 0]; p <- sort(p[!is.na(p)])
    data.frame(p = unname(p), ecdf = seq_along(p) / length(p),
               method = unname(MLAB[m]))
  }))
  p1 <- ggplot(ecdf_df, aes(p, ecdf, colour = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    geom_step(linewidth = 0.9) +
    scale_colour_manual(values = PAL, name = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Null genes: p-value calibration",
         subtitle = "Dashed line = uniform (correct type-I error). Above it = anti-conservative.",
         x = "Nominal p-value", y = "Empirical CDF") +
    theme_pub(11, legend = c(0.98, 0.02)) +
    theme(legend.justification = c(1, 0))
  save_fig(p1, FIG("fig1_null_pvalue_calibration"), width = 5.6, height = 4.6)
}

# 图 2:ROC 曲线(折线)
if (has_truth) {
  roc_df <- do.call(rbind, lapply(names(METHODS), function(m) {
    r <- roc_of(-log10(res[[METHODS[m]]]), res$is_de)
    a <- auc_of(-log10(res[[METHODS[m]]]), res$is_de)
    r$method <- sprintf("%s (AUC %.3f)", MLAB[m], a); r$key <- MLAB[m]; r
  }))
  rp <- setNames(PAL[unique(roc_df$key)], unique(roc_df$method))
  p2 <- ggplot(roc_df, aes(fpr, tpr, colour = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = rp, name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Recovery of true DE genes", x = "False positive rate", y = "True positive rate") +
    theme_pub(11, legend = c(0.98, 0.02)) +
    theme(legend.justification = c(1, 0))
  save_fig(p2, FIG("fig2_roc_true_de"), width = 5.2, height = 5.0)
}

# 图 3:naive vs GLMM 的 -log10 p 散点,按真值着色
if (has_truth) {
  sc <- data.frame(x = -log10(res$naive_pval), y = -log10(res$glmm_pval),
                   cls = res$gene_class, gene = res$gene)
  sc <- sc[stats::complete.cases(sc[, c("x", "y")]), ]
  cls_pal <- setNames(pal_pub(name = "okabe_ito")[c(6, 5, 3)],
                      c("true_de", "donor_driven_null", "plain_null"))
  p3 <- ggplot(sc, aes(x, y, colour = cls)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey55") +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", colour = "grey55") +
    geom_point(size = 2.1, alpha = 0.85) +
    scale_colour_manual(values = cls_pal, name = "Gene class",
                        labels = c(true_de = "True DE",
                                   donor_driven_null = "Null, donor-driven",
                                   plain_null = "Null, plain")) +
    labs(title = "Where the naive test goes wrong",
         subtitle = "Points right of the diagonal are called by the cell-level test but not by the GLMM",
         x = expression(-log[10]~"p  ·  naive cell-level t-test"),
         y = expression(-log[10]~"p  ·  Poisson-GLMM (donor RE)")) +
    theme_pub(11)
  save_fig(p3, FIG("fig3_naive_vs_glmm_scatter"), width = 6.2, height = 4.8)
}

# 图 4:哑铃图 —— naive 与 GLMM 对同一批"供体驱动空基因"的显著性落差
if (has_truth) {
  cand <- res[res$gene_class == "donor_driven_null" & !is.na(res$naive_pval), ]
  cand <- cand[order(cand$naive_pval), ][seq_len(min(20, nrow(cand))), ]
  cand$gene <- factor(cand$gene, levels = rev(cand$gene))
  dd <- data.frame(gene = rep(cand$gene, 2),
                   val = c(-log10(cand$naive_pval), -log10(cand$glmm_pval)),
                   method = rep(c(MLAB["naive"], MLAB["glmm"]), each = nrow(cand)))
  p4 <- ggplot() +
    geom_segment(data = cand, aes(y = gene, yend = gene,
                                  x = -log10(naive_pval), xend = -log10(glmm_pval)),
                 colour = "grey70", linewidth = 0.7) +
    geom_point(data = dd, aes(val, gene, colour = method), size = 2.6) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", colour = "grey40") +
    scale_colour_manual(values = PAL[c(MLAB["naive"], MLAB["glmm"])], name = NULL) +
    labs(title = "Donor-driven null genes lose their false signal",
         subtitle = "Top 20 null genes ranked by the naive test; dotted line = p 0.05",
         x = expression(-log[10]~"p-value"), y = NULL) +
    theme_pub(10, legend = "top")
  save_fig(p4, FIG("fig4_dumbbell_null_genes"), width = 5.8, height = 6.0)
}

# 图 5:单基因供体级 raincloud/点图 —— 展示"细胞多但供体少"的伪重复
top_fp <- if (has_truth) {
  cd <- res[res$gene_class == "donor_driven_null" & !is.na(res$naive_pval), ]
  cd$gene[which.min(cd$naive_pval)]
} else res$gene[which.min(res$naive_pval)]
gi <- match(top_fp, gene_ids)
cell_df <- data.frame(expr = log2(counts[gi, ] + 1), donor = as.character(donor),
                      group = as.character(grp))
don_df <- stats::aggregate(expr ~ donor + group, data = cell_df, FUN = mean)
p5 <- ggplot(cell_df, aes(donor, expr, colour = group)) +
  geom_violin(aes(fill = group), colour = NA, alpha = 0.20, width = 0.9,
              scale = "width", show.legend = FALSE) +
  geom_jitter(width = 0.16, height = 0, size = 0.9, alpha = 0.45, show.legend = FALSE) +
  geom_point(data = don_df, aes(donor, expr, fill = group), shape = 23, size = 3.4,
             colour = "black", stroke = 0.5) +
  scale_colour_manual(values = pal_pub(2, "nejm")) +
  scale_fill_manual(values = pal_pub(2, "nejm"), name = args$group) +
  labs(title = sprintf("Gene %s: signal lives at the donor level, not the cell level", top_fp),
       subtitle = sprintf("naive p = %.2g   ·   GLMM p = %.2g   (diamonds = donor means)",
                          res$naive_pval[gi], res$glmm_pval[gi]),
       x = "Donor", y = expression(log[2]*"(count + 1)")) +
  theme_pub(10, legend = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p5, FIG("fig5_donor_level_variation"), width = 7.0, height = 4.6)

# 图 6:GLMM 供体随机效应方差 vs 朴素检验显著性(散点)
sig_df <- data.frame(sigma2 = glmm$sigma_square, naive = -log10(naive$pval),
                     cls = if (has_truth) res$gene_class else "all")
sig_df <- sig_df[stats::complete.cases(sig_df), ]
p6 <- ggplot(sig_df, aes(sigma2, naive, colour = cls)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey45") +
  geom_point(size = 2.1, alpha = 0.85) +
  scale_x_continuous(trans = "log10") +
  { if (has_truth) scale_colour_manual(values = pal_pub(name = "okabe_ito")[c(6, 5, 3)],
        breaks = c("true_de", "donor_driven_null", "plain_null"), name = "Gene class",
        labels = c("True DE", "Null, donor-driven", "Null, plain"))
    else scale_colour_manual(values = "grey30", guide = "none") } +
  labs(title = "Donor dispersion drives naive false positives",
       subtitle = "x = GLMM residual dispersion; higher dispersion inflates the cell-level test",
       x = expression("GLMM "*sigma^2*" (log"[10]*" scale)"),
       y = expression(-log[10]~"p  ·  naive cell-level t-test")) +
  theme_pub(11)
save_fig(p6, FIG("fig6_dispersion_vs_naive_significance"), width = 6.2, height = 4.6)

# ---- 依赖版本快照(可复现性)------------------------------------------------
utils::capture.output(utils::sessionInfo(),
                      file = file.path(args$outdir, "sessionInfo.txt"))

cat("\n完成。结果表 →", normalizePath(args$outdir, mustWork = FALSE),
    "\n展示图 →", normalizePath(args$assets, mustWork = FALSE), "\n")
