# =============================================================================
# 594 · ColocBoost — 梯度提升的多性状(multi-trait)共定位
# -----------------------------------------------------------------------------
# 用途   : 一个基因座上有 GWAS + 多个分子 QTL(eQTL/sQTL/pQTL)时,判断哪些性状
#          子集共享同一个因果变异。经典 coloc 只能两两跑且假设"每性状至多 1 个
#          因果变异";ColocBoost 用多任务梯度提升一次性对所有性状联合建模。
#
# 诚实说明(本库硬规矩):
#   · 本机 **未安装** colocboost(不允许装包)。脚本对真包做 requireNamespace 守卫:
#     装了就按官方真实签名真跑并落盘;没装则只跑基线,并打印真实安装命令。
#   · **基线是真包 `coloc::coloc.abf`(本机已装)**,不是我自己重写的近似 —— 即
#     ColocBoost 声称要超越的 pairwise ABF 框架本身。因此无论真包在不在,
#     `Rscript 594_colocboost_colocalization.R` 都能跑完并出全部图。
#
# 真实 API 来源(2026-07-21 复核:直接读 StatFunGen/colocboost v1.0.9 源码,非网页转述):
#   · colocboost()            定义于 R/colocboost.R:143;形参 sumstat=/LD= 见 :144
#   · sumstat 列 z/n/variant  见 vignettes/Input_Data_Format.Rmd:49,55-57
#   · LD 需带 dimnames        见 vignettes/Input_Data_Format.Rmd:60
#   · $cos_summary            R/colocboost.R:107 @return · 赋值 R/colocboost_assemble.R:238
#   · $vcp(带 variant 名)    R/colocboost.R:108 @return · 计算 R/colocboost_output.R:241-243
#   · $cos_details$cos$cos_index  R/colocboost.R:135(官方 roxygen 示例原句)
#   核心调用: colocboost(sumstat = <list of data.frame(z, n, variant)>, LD = <矩阵,带 dimnames>)
#
# 输入   : example_data/sumstat_<trait>.csv (variant,beta,se,z,n,maf,pos) + region_ld.csv
# 输出   : results/ (csv + versions.txt) · 展示图 assets/
# 出图   : dot(区域图) / heatmap(两两 PP.H4) / dumbbell(H3 vs H4) / lollipop(变异级)
#          —— 无条形图(库内硬规矩)
# 随机种子: SEED = 42
# =============================================================================

suppressWarnings(suppressMessages({
  library(ggplot2)
}))

# ---- 定位脚本目录 + 载入统一顶刊主题 ----------------------------------------
.script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", a[m[1]]))))
  getwd()
}
ROOT <- .script_dir()
source(file.path(ROOT, "..", "..", "..", "_framework", "theme_pub.R"))

# ---- 参数区(全部支持 --key value 覆盖)-------------------------------------
opt <- bio_args(list(
  datadir  = file.path(ROOT, "example_data"),
  outdir   = file.path(ROOT, "results"),
  assets   = file.path(ROOT, "assets"),
  focal    = "GWAS",   # 用于变异级图的焦点性状
  seed     = "42"
))
SEED <- as.integer(opt$seed)
set.seed(SEED)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$assets, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c("GWAS", "eQTL_geneA", "sQTL_geneA", "pQTL_geneB")

# =============================================================================
# Step 0 · 合成示例数据(仅当 example_data 缺失时生成;synthetic, for demo only)
# -----------------------------------------------------------------------------
# 设计意图:一个区域内放 **两个独立因果变异**(两个 LD 块各一个),
#   causal_1 (pos 60) 被 GWAS / eQTL_geneA / sQTL_geneA 共享  → 3 性状共定位
#   causal_2 (pos 150) 被 sQTL_geneA / pQTL_geneB 共享        → 另一组 2 性状共定位
# 这正是 pairwise coloc 的痛点:sQTL 有两个信号 → 单因果假设被违反 →
#   GWAS vs sQTL 的 PP.H4 被稀释,而真值上它们确实共享 causal_1。
# =============================================================================
make_example_data <- function(datadir) {
  cat("Step 0 · 生成合成示例数据 (synthetic, for demo only)\n")
  dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
  set.seed(SEED)

  P <- 200L; N <- 3000L
  # 两个 LD 块(块内 AR(1) 强相关,块间独立)—— 模拟真实基因座的 haplotype 结构
  blk <- c(rep(1L, 100L), rep(2L, 100L))
  ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
  Sig <- matrix(0, P, P)
  Sig[blk == 1, blk == 1] <- ar1(sum(blk == 1), 0.93)
  Sig[blk == 2, blk == 2] <- ar1(sum(blk == 2), 0.93)

  # 潜变量 → 两条单倍型阈值化 → 真实 0/1/2 dosage(有真实 MAF,不是连续伪基因型)
  maf <- runif(P, 0.08, 0.45)
  cut <- stats::qnorm(1 - maf)
  L <- chol(Sig + diag(1e-8, P))
  draw_hap <- function() (matrix(rnorm(N * P), N, P) %*% L > matrix(cut, N, P, byrow = TRUE)) * 1L
  X <- draw_hap() + draw_hap()
  colnames(X) <- sprintf("rs%03d", seq_len(P))
  keep <- apply(X, 2, sd) > 0
  X <- X[, keep, drop = FALSE]; maf <- maf[keep]
  P <- ncol(X)

  causal <- c(causal_1 = 60L, causal_2 = 150L)   # 列索引
  # 各性状的真效应(标准化尺度)
  eff <- list(
    GWAS       = c(causal_1 = 0.16, causal_2 = 0.00),
    eQTL_geneA = c(causal_1 = 0.28, causal_2 = 0.00),
    sQTL_geneA = c(causal_1 = 0.22, causal_2 = 0.24),  # 两个独立信号 → 违反单因果假设
    pQTL_geneB = c(causal_1 = 0.00, causal_2 = 0.26)
  )
  Xs <- scale(X)

  for (tr in TRAITS) {
    b <- eff[[tr]]
    y <- Xs[, causal["causal_1"]] * b["causal_1"] +
         Xs[, causal["causal_2"]] * b["causal_2"] + rnorm(N)
    y <- as.numeric(scale(y))                     # 标准化 → sdY = 1(coloc 用得上)
    # 逐变异边际回归(单变量 GWAS 的标准做法)
    st <- t(apply(X, 2, function(g) {
      f <- stats::lm.fit(cbind(1, g), y); r <- f$residuals
      s2 <- sum(r^2) / (N - 2)
      gc_ <- g - mean(g); se <- sqrt(s2 / sum(gc_^2))
      c(beta = unname(f$coefficients[2]), se = se)
    }))
    df <- data.frame(
      variant = colnames(X), pos = seq_len(P), maf = round(maf, 4),
      beta = round(st[, "beta"], 6), se = round(st[, "se"], 6),
      z = round(st[, "beta"] / st[, "se"], 4), n = N, row.names = NULL
    )
    write.csv(df, file.path(datadir, sprintf("sumstat_%s.csv", tr)), row.names = FALSE)
  }
  LD <- stats::cor(X)
  write.csv(round(LD, 4), file.path(datadir, "region_ld.csv"), row.names = TRUE)
  write.csv(data.frame(signal = names(causal), variant = colnames(X)[causal], pos = unname(causal)),
            file.path(datadir, "true_causal.csv"), row.names = FALSE)
  cat(sprintf("   已写入 %d 变异 × %d 性状 + LD 矩阵\n", P, length(TRAITS)))
}

if (!file.exists(file.path(opt$datadir, "region_ld.csv"))) make_example_data(opt$datadir)

# ---- 读入 -------------------------------------------------------------------
cat("Step 1 · 读入区域 summary statistics 与 LD\n")
SS <- lapply(TRAITS, function(tr)
  read.csv(file.path(opt$datadir, sprintf("sumstat_%s.csv", tr)), stringsAsFactors = FALSE))
names(SS) <- TRAITS
LD <- as.matrix(read.csv(file.path(opt$datadir, "region_ld.csv"), row.names = 1, check.names = FALSE))
colnames(LD) <- rownames(LD)
TRUE_CAUSAL <- read.csv(file.path(opt$datadir, "true_causal.csv"), stringsAsFactors = FALSE)
cat(sprintf("   %d 性状 · %d 变异 · 真因果 %s\n", length(SS), nrow(SS[[1]]),
            paste(TRUE_CAUSAL$variant, collapse = ", ")))

# =============================================================================
# Step 2 · 基线:真包 coloc::coloc.abf 的两两共定位(Giambartolomei et al. 2014)
# -----------------------------------------------------------------------------
# 这是 ColocBoost 声称要超越的框架本身,用真包跑,不做近似重写。
# =============================================================================
cat("Step 2 · 基线 pairwise coloc (coloc::coloc.abf)\n")
has_coloc <- requireNamespace("coloc", quietly = TRUE)

mk_ds <- function(df) list(beta = df$beta, varbeta = df$se^2, snp = df$variant,
                           type = "quant", N = df$n[1], sdY = 1, MAF = df$maf)

pairs_df <- t(combn(TRAITS, 2))
base_rows <- list(); snp_pp <- list()
if (has_coloc) {
  for (i in seq_len(nrow(pairs_df))) {
    a <- pairs_df[i, 1]; b <- pairs_df[i, 2]
    r <- suppressWarnings(suppressMessages(
      coloc::coloc.abf(dataset1 = mk_ds(SS[[a]]), dataset2 = mk_ds(SS[[b]]))))
    s <- r$summary
    base_rows[[i]] <- data.frame(trait1 = a, trait2 = b,
      PP.H0 = unname(s["PP.H0.abf"]), PP.H1 = unname(s["PP.H1.abf"]),
      PP.H2 = unname(s["PP.H2.abf"]), PP.H3 = unname(s["PP.H3.abf"]),
      PP.H4 = unname(s["PP.H4.abf"]))
    if (!is.null(r$results) && "SNP.PP.H4" %in% names(r$results))
      snp_pp[[paste(a, b, sep = "|")]] <-
        data.frame(pair = paste(a, b, sep = " vs "), variant = r$results$snp,
                   SNP.PP.H4 = r$results$SNP.PP.H4)
  }
  BASE <- do.call(rbind, base_rows)
} else {
  cat("   ! 未检测到 coloc 包 —— 基线不可跑。安装: install.packages('coloc')\n")
  BASE <- data.frame(trait1 = pairs_df[, 1], trait2 = pairs_df[, 2],
                     PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA)
}
write.csv(BASE, file.path(opt$outdir, "baseline_pairwise_coloc.csv"), row.names = FALSE)
print(BASE[, c("trait1", "trait2", "PP.H3", "PP.H4")], digits = 3)

# =============================================================================
# Step 3 · ColocBoost(守卫式;本机未装则跳过,绝不伪造结果)
# -----------------------------------------------------------------------------
# 官方真实签名(直接读本地克隆源码 R/colocboost.R:143-201,非网页转述):
#   colocboost(X=NULL, Y=NULL, sumstat=NULL, LD=NULL, X_ref=NULL, ...,
#              focal_outcome_idx=NULL, output_level=1, coverage=0.95, ...)
#   逐个默认值出处: focal_outcome_idx=NULL :150 · M=500 :160 · tau=0.01 :162 ·
#   learning_rate_init=0.01 :163 · coverage=0.95 :182 · output_level=1 :201
# sumstat 每个元素需列 `z`(或 `beta`+`sebeta`)、`n`、`variant`;LD 需带 dimnames。
# =============================================================================
cat("Step 3 · ColocBoost(真包路径,守卫)\n")
CB_STATUS <- "NOT-INSTALLED"; CB <- NULL
if (requireNamespace("colocboost", quietly = TRUE)) {
  sumstat_list <- lapply(SS, function(df) data.frame(z = df$z, n = df$n, variant = df$variant))
  names(sumstat_list) <- TRAITS
  CB <- try(colocboost::colocboost(sumstat = sumstat_list, LD = LD), silent = TRUE)
  if (inherits(CB, "try-error")) {
    CB_STATUS <- paste0("INSTALLED-BUT-FAILED: ", trimws(as.character(CB)))
    cat("   ! 真包调用失败:", CB_STATUS, "\n"); CB <- NULL
  } else {
    CB_STATUS <- "REAL-RUN-OK"
    if (!is.null(CB$cos_summary))
      write.csv(CB$cos_summary, file.path(opt$outdir, "colocboost_cos_summary.csv"), row.names = FALSE)
    vcp <- CB$vcp
    if (!is.null(vcp))
      write.csv(data.frame(variant = if (!is.null(names(vcp))) names(vcp) else SS[[1]]$variant,
                           vcp = as.numeric(vcp)),
                file.path(opt$outdir, "colocboost_vcp.csv"), row.names = FALSE)
    cat("   真包运行成功,cos_summary / vcp 已落盘\n")
  }
} else {
  cat("   ⏭ 本机未安装 colocboost —— 只跑基线。真实安装命令:\n",
      "     install.packages('colocboost')                       # CRAN 稳定版\n",
      "     devtools::install_github('StatFunGen/colocboost')     # 开发版\n", sep = "")
}

# =============================================================================
# Step 4 · 出图(全部非条形图)
# =============================================================================
cat("Step 4 · 出图\n")
PAL <- pal_pub(4, "npg")
causal_pos <- TRUE_CAUSAL$pos

# --- 图1 区域多性状散点(dot)------------------------------------------------
reg <- do.call(rbind, lapply(TRAITS, function(tr) {
  d <- SS[[tr]]
  data.frame(trait = tr, pos = d$pos, variant = d$variant,
             logp = -log10(pmax(2 * pnorm(-abs(d$z)), 1e-300)))
}))
reg$trait <- factor(reg$trait, levels = TRAITS)
reg$is_causal <- reg$variant %in% TRUE_CAUSAL$variant
p1 <- ggplot(reg, aes(pos, logp)) +
  geom_vline(xintercept = causal_pos, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_point(aes(colour = trait), size = 1.5, alpha = 0.75) +
  geom_point(data = subset(reg, is_causal), shape = 8, size = 3, colour = "black", stroke = 0.7) +
  facet_wrap(~trait, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = PAL, guide = "none") +
  labs(x = "Variant position in region (index)", y = expression(-log[10](italic(P))),
       title = "Regional association across four traits",
       subtitle = "Dashed lines = simulated causal variants (star). Two independent signals in the locus.") +
  theme_pub(base_size = 11)
save_fig(p1, file.path(opt$assets, "region_multitrait_dots"), width = 9, height = 6.2)

# --- 图2 两两 PP.H4 热图(heatmap)-------------------------------------------
hm <- rbind(
  data.frame(t1 = BASE$trait1, t2 = BASE$trait2, v = BASE$PP.H4),
  data.frame(t1 = BASE$trait2, t2 = BASE$trait1, v = BASE$PP.H4))
hm$t1 <- factor(hm$t1, levels = TRAITS); hm$t2 <- factor(hm$t2, levels = rev(TRAITS))
p2 <- ggplot(hm, aes(t1, t2, fill = v)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", v)),
            colour = ifelse(hm$v > 0.6, "white", "black"), size = 3.6) +
  scale_fill_viridis_c(option = "D", limits = c(0, 1), name = "PP.H4") +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "Pairwise colocalization posterior (coloc.abf baseline)",
       subtitle = "PP.H4 = shared single causal variant. Six independent pairwise tests.") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.line = element_blank(),
        axis.ticks = element_blank())
save_fig(p2, file.path(opt$assets, "pairwise_pph4_heatmap"), width = 6.6, height = 5.6)

# --- 图3 H3 vs H4 哑铃图(dumbbell)------------------------------------------
db <- BASE
db$pair <- paste(db$trait1, "vs", db$trait2)
db$pair <- factor(db$pair, levels = db$pair[order(db$PP.H4)])
p3 <- ggplot(db) +
  geom_segment(aes(y = pair, yend = pair, x = PP.H3, xend = PP.H4),
               colour = "grey65", linewidth = 1.1) +
  geom_point(aes(y = pair, x = PP.H3, colour = "PP.H3 (distinct variants)"), size = 3.4) +
  geom_point(aes(y = pair, x = PP.H4, colour = "PP.H4 (shared variant)"), size = 3.4) +
  scale_colour_manual(values = c("PP.H3 (distinct variants)" = PAL[4],
                                 "PP.H4 (shared variant)" = PAL[1]), name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Posterior probability", y = NULL,
       title = "Pairwise coloc: competing hypotheses per trait pair",
       subtitle = "Pairs involving sQTL_geneA carry two independent signals, violating the one-causal-variant assumption.") +
  theme_pub(base_size = 11, legend = "bottom")
save_fig(p3, file.path(opt$assets, "pph_dumbbell"), width = 8.4, height = 5.0)

# --- 图4 变异级共定位概率棒棒糖(lollipop)------------------------------------
# 基线量 = coloc 的 SNP.PP.H4(真包输出列);若 ColocBoost 跑通则叠加其 VCP。
lolli_src <- NULL
if (length(snp_pp)) {
  fk <- grep(opt$focal, names(snp_pp), fixed = TRUE)
  key <- if (length(fk)) names(snp_pp)[fk[1]] else names(snp_pp)[1]
  d <- snp_pp[[key]]
  d <- d[order(-d$SNP.PP.H4), ][1:15, ]
  d$metric <- sprintf("coloc SNP.PP.H4 (%s)", d$pair[1])
  lolli_src <- data.frame(variant = d$variant, value = d$SNP.PP.H4, metric = d$metric)
}
if (!is.null(CB) && !is.null(CB$vcp)) {
  v <- as.numeric(CB$vcp)
  nm <- if (!is.null(names(CB$vcp))) names(CB$vcp) else SS[[1]]$variant
  o <- order(-v)[1:15]
  lolli_src <- rbind(lolli_src,
    data.frame(variant = nm[o], value = v[o], metric = "ColocBoost VCP (real package)"))
}
if (!is.null(lolli_src)) {
  lolli_src$variant <- factor(lolli_src$variant,
    levels = unique(lolli_src$variant[order(lolli_src$value)]))
  lolli_src$is_causal <- as.character(lolli_src$variant) %in% TRUE_CAUSAL$variant
  p4 <- ggplot(lolli_src, aes(value, variant)) +
    geom_segment(aes(x = 0, xend = value, yend = variant), colour = "grey72", linewidth = 0.7) +
    geom_point(aes(colour = is_causal), size = 3.2) +
    scale_colour_manual(values = c(`FALSE` = PAL[2], `TRUE` = PAL[1]),
                        labels = c("other variant", "true causal variant"), name = NULL) +
    facet_wrap(~metric, ncol = 2, scales = "free_y") +
    labs(x = "Variant-level colocalization probability", y = NULL,
         title = "Which variant drives the shared signal?",
         subtitle = "Top 15 variants. Red = simulated causal variant.") +
    theme_pub(base_size = 10, legend = "bottom")
  save_fig(p4, file.path(opt$assets, "variant_level_lollipop"),
           width = ifelse(length(unique(lolli_src$metric)) > 1, 9, 5.6), height = 5.6)
}

# =============================================================================
# Step 5 · 依赖快照与小结
# =============================================================================
si <- utils::sessionInfo()
writeLines(c(
  paste0("R: ", si$R.version$version.string),
  paste0("coloc: ", if (has_coloc) as.character(utils::packageVersion("coloc")) else "NOT-INSTALLED"),
  paste0("colocboost: ", CB_STATUS),
  paste0("ggplot2: ", as.character(utils::packageVersion("ggplot2"))),
  paste0("seed: ", SEED)
), file.path(opt$outdir, "versions.txt"))

cat("\n完成。\n")
cat("  results/ ->", paste(list.files(opt$outdir), collapse = ", "), "\n")
cat("  assets/  ->", paste(grep("\\.png$", list.files(opt$assets), value = TRUE), collapse = ", "), "\n")
cat("  ColocBoost 状态:", CB_STATUS, "\n")
