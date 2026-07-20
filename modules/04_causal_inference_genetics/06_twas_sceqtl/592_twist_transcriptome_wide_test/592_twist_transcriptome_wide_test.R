# =============================================================================
# 592 · TWiST — 细胞状态(拟时序)层面的转录组关联分析
# -----------------------------------------------------------------------------
# 上游方法: TWiST ("TWAS in pseudoTime"), Qi G, Lila E, Ji Z, Shojaie A,
#           Battle A, Sun W. Cell Genomics 2026 Jan 14.
#           doi:10.1016/j.xgen.2025.101060 · PMID 41187759
#           repo: https://github.com/gqi/TWiST
#
# 本模块做两件事:
#  (A) 【基线,本机零依赖可跑】经典"静态 TWAS burden 检验" vs "拟时序分辨扫描"。
#      这是 TWiST 想要取代的朴素对照:传统 TWAS 把一个细胞类型压成一个权重向量,
#      沿拟时序反向的基因-性状效应会互相抵消。本基线用同样的 B-spline 权重矩阵
#      Wmat 在多个拟时序点上分别做 TWAS,直接把"静态会漏掉什么"量化出来。
#      基线不是 TWiST,也不冒充 TWiST 的似然比检验。
#  (B) 【守卫式封装】真正的 TWiST 三联检验(global / dynamic / nonlinear),
#      需要 R 包 TWiST + plink2R + fda + grpreg。未安装时优雅退出并打印真实安装
#      命令与已核实的函数签名(签名逐字取自上游 man/*.Rd 与 example_data/example.R)。
#
# 零改动即跑:  Rscript 592_twist_transcriptome_wide_test.R
# 换数据:      Rscript 592_twist_transcriptome_wide_test.R --outdir results/run1
# =============================================================================

suppressPackageStartupMessages({
  library(splines)   # base R,提供 bs() B-spline 基,与上游 twist_train_model 同一套
  library(ggplot2)
})

# ---- 定位与框架 -------------------------------------------------------------
HERE <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
})
source(file.path(HERE, "..", "..", "..", "_framework", "theme_pub.R"))

args <- bio_args(list(
  datadir = file.path(HERE, "example_data"),
  outdir  = file.path(HERE, "results"),
  npt     = "21",     # 拟时序扫描网格点数
  seed    = "592"
))
DATADIR <- args$datadir
OUTDIR  <- args$outdir
ASSETS  <- file.path(HERE, "assets")
NPT     <- as.integer(args$npt)
SEED    <- as.integer(args$seed)
RUN_TWIST <- isTRUE(args$`run-twist`)
REGEN     <- isTRUE(args$`regen-example`)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

# TWiST stage-1 默认基:内部节点 c(.25,.5,.75)、三次 B-spline、含截距 -> 7 个基函数
KNOTS_W  <- c(0.25, 0.50, 0.75)
DEGREE_W <- 3
bs_w <- function(t) bs(t, knots = KNOTS_W, degree = DEGREE_W,
                       Boundary.knots = c(0, 1), intercept = TRUE)
NBASIS <- ncol(bs_w(0.5))

# =============================================================================
# Step 0 · 合成示例数据(格式对齐上游真实输入)
# =============================================================================
gen_example <- function(datadir, ngene = 60, nsnp = 50, ngwas = 6.9e4,
                        rho = 0.7, eff = 0.045) {
  dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
  cat("Step 0: 生成合成示例数据 (synthetic, for demo only)\n")

  # --- LD 矩阵:单位点 AR(1),保证半正定 ---
  idx <- seq_len(nsnp)
  R <- rho ^ abs(outer(idx, idx, "-"))
  cholR <- chol(R)
  snp_id <- sprintf("rs%06d", 900000 + idx)
  alleles <- c("A", "C", "G", "T")
  a1 <- sample(alleles, nsnp, TRUE)
  a2 <- vapply(a1, function(x) sample(setdiff(alleles, x), 1), character(1))

  # --- 每个基因的情景:null / static / dynamic ---
  scen <- rep(c("null", "static", "dynamic"), length.out = ngene)
  gene_id <- sprintf("GENE%03d", seq_len(ngene))

  tgrid <- seq(0.001, 0.999, length.out = 400)   # 数值积分网格
  B <- bs_w(tgrid)                               # 400 x NBASIS

  wmat_rows <- list(); z_all <- numeric(0)
  for (g in seq_len(ngene)) {
    # eQTL 权重矩阵 Wmat: nsnp x NBASIS,稀疏(5 个 cis-eQTL)
    W <- matrix(0, nsnp, NBASIS)
    causal <- sort(sample(idx, 5))
    W[causal, ] <- matrix(rnorm(5 * NBASIS, 0, 0.35), 5, NBASIS)
    # 让权重沿拟时序平滑(相邻基函数相关),更像真实 stage-1 输出
    W[causal, ] <- t(apply(W[causal, , drop = FALSE], 1, function(r) cumsum(r) / sqrt(seq_along(r))))

    # 基因-性状效应曲线 beta(t)
    beta_t <- switch(scen[g],
      null    = rep(0, length(tgrid)),
      static  = rep(eff, length(tgrid)),
      # 沿拟时序变号 -> 静态 TWAS 的平均权重与真实效应方向近乎正交,信号被抵消。
      # 幅度取 6x 使其在单个拟时序点上的 burden 信号与 static 情景可比,
      # 这样两条臂的差别只来自"要不要沿细胞状态分辨",而不是来自效应量大小。
      dynamic = eff * 6.0 * (2 * tgrid - 1)
    )
    # 性状上的 SNP 效应: a = Wmat %*% ∫ b(t)beta(t) dt
    cvec <- as.vector(crossprod(B, beta_t) * mean(diff(tgrid)))
    a_snp <- as.vector(W %*% cvec)

    # GWAS 边际 z:信号经 LD 传播 + 相关噪声
    z <- sqrt(ngwas) * as.vector(R %*% a_snp) + as.vector(crossprod(cholR, rnorm(nsnp)))
    z_all <- rbind(z_all, z)

    wmat_rows[[g]] <- data.frame(
      ID = gene_id[g], SNP = snp_id,
      setNames(as.data.frame(W), sprintf("basis%d", seq_len(NBASIS))),
      stringsAsFactors = FALSE
    )
  }

  # 示例里每个基因占一段独立的 cis 区间(SNP 名带基因前缀),因此各自一套 sumstats 分片
  sumstat <- do.call(rbind, lapply(seq_len(ngene), function(g) {
    data.frame(SNP = paste0(gene_id[g], ":", snp_id), A1 = a1, A2 = a2,
               Z = z_all[g, ], stringsAsFactors = FALSE)
  }))
  wgt <- do.call(rbind, wmat_rows)
  wgt$SNP <- paste0(wgt$ID, ":", wgt$SNP)

  wgtlist <- data.frame(
    ID = gene_id, CHR = 6L,
    P0 = 1e6 + seq_len(ngene) * 3e5,
    P1 = 1e6 + seq_len(ngene) * 3e5 + 2e4,
    tss = 1e6 + seq_len(ngene) * 3e5,
    scenario = scen, stringsAsFactors = FALSE
  )

  hdr <- "# synthetic, for demo only — 592 TWiST module\n"
  writeLines(c(hdr, ""), file.path(datadir, "README_synthetic.txt"))
  write.table(sumstat, file.path(datadir, "gwas_sumstats_synth.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(wgt, file.path(datadir, "twist_weights_Wmat_synth.csv"), row.names = FALSE)
  write.csv(wgtlist, file.path(datadir, "wgtlist_synth.csv"), row.names = FALSE)
  write.csv(as.data.frame(R), file.path(datadir, "ld_reference_synth.csv"), row.names = FALSE)
  writeLines(as.character(ngwas), file.path(datadir, "ngwas_synth.txt"))
  invisible(TRUE)
}

need <- c("gwas_sumstats_synth.txt", "twist_weights_Wmat_synth.csv",
          "wgtlist_synth.csv", "ld_reference_synth.csv", "ngwas_synth.txt")
if (REGEN || !all(file.exists(file.path(DATADIR, need)))) gen_example(DATADIR)

# =============================================================================
# Step 1 · 读入
# =============================================================================
cat("Step 1: 读入 GWAS 汇总统计 / Wmat 权重 / LD 参考\n")
sumstat <- read.delim(file.path(DATADIR, "gwas_sumstats_synth.txt"), stringsAsFactors = FALSE)
wgt     <- read.csv(file.path(DATADIR, "twist_weights_Wmat_synth.csv"), stringsAsFactors = FALSE)
wgtlist <- read.csv(file.path(DATADIR, "wgtlist_synth.csv"), stringsAsFactors = FALSE)
R       <- as.matrix(read.csv(file.path(DATADIR, "ld_reference_synth.csv")))
ngwas   <- as.numeric(readLines(file.path(DATADIR, "ngwas_synth.txt"))[1])
basis_cols <- grep("^basis", names(wgt), value = TRUE)
cat(sprintf("       genes=%d  SNPs/gene=%d  basis=%d  n_eff=%.0f\n",
            nrow(wgtlist), nrow(R), length(basis_cols), ngwas))

# =============================================================================
# Step 2 · 基线:静态 TWAS burden 检验 + 拟时序分辨扫描
# -----------------------------------------------------------------------------
# 静态 TWAS(FUSION 式):z = w'z_gwas / sqrt(w'Rw),w 为拟时序平均权重。
# 拟时序扫描:w(t) = Wmat %*% b(t)(b 为 stage-1 的 B-spline 基),在每个 t 上
#              各做一次同样的 burden 检验;跨 t 的最小 p 用 Bonferroni 校正
#              (网格点高度相关,故该校正偏保守,是对基线不利的诚实设定)。
# =============================================================================
cat("Step 2: 基线 — 静态 TWAS vs 拟时序分辨 TWAS 扫描\n")
tgrid <- seq(0.02, 0.98, length.out = NPT)
Bt <- bs_w(tgrid)                      # NPT x NBASIS
bbar <- colMeans(Bt)                   # 拟时序平均基 -> 静态权重

burden_z <- function(w, z, R) {
  den <- as.numeric(crossprod(w, R %*% w))
  if (!is.finite(den) || den <= 1e-12) return(NA_real_)
  as.numeric(crossprod(w, z)) / sqrt(den)
}

res <- data.frame(); zmat <- matrix(NA_real_, nrow(wgtlist), NPT,
                                    dimnames = list(wgtlist$ID, sprintf("t%.2f", tgrid)))
for (i in seq_len(nrow(wgtlist))) {
  g <- wgtlist$ID[i]
  wg <- wgt[wgt$ID == g, ]
  ss <- sumstat[match(wg$SNP, sumstat$SNP), ]
  W <- as.matrix(wg[, basis_cols]); z <- ss$Z
  keep <- !is.na(z)
  W <- W[keep, , drop = FALSE]; z <- z[keep]; Rk <- R[keep, keep, drop = FALSE]

  z_static <- burden_z(as.vector(W %*% bbar), z, Rk)
  zt <- vapply(seq_len(NPT), function(k) burden_z(as.vector(W %*% Bt[k, ]), z, Rk), numeric(1))
  zmat[i, ] <- zt

  p_static <- 2 * pnorm(-abs(z_static))
  p_scan   <- min(pmin(1, 2 * pnorm(-abs(zt)) * NPT), na.rm = TRUE)  # Bonferroni over grid
  res <- rbind(res, data.frame(
    ID = g, scenario = wgtlist$scenario[i], CHR = wgtlist$CHR[i], tss = wgtlist$tss[i],
    z_static = z_static, p_static = p_static,
    t_peak = tgrid[which.max(abs(zt))], z_peak = zt[which.max(abs(zt))],
    p_scan_bonf = p_scan,
    stringsAsFactors = FALSE))
}
res$fdr_static <- p.adjust(res$p_static, "BH")
res$fdr_scan   <- p.adjust(res$p_scan_bonf, "BH")
res$gain_log10p <- -log10(res$p_scan_bonf) - (-log10(res$p_static))

write.csv(res, file.path(OUTDIR, "592_baseline_gene_results.csv"), row.names = FALSE)
zdf <- data.frame(ID = rownames(zmat), zmat, check.names = FALSE)
write.csv(zdf, file.path(OUTDIR, "592_pseudotime_zscore_matrix.csv"), row.names = FALSE)

ALPHA <- 0.05
pw <- do.call(rbind, lapply(split(res, res$scenario), function(d) data.frame(
  scenario = d$scenario[1],
  static  = mean(d$fdr_static < ALPHA),
  scan    = mean(d$fdr_scan   < ALPHA),
  n = nrow(d), stringsAsFactors = FALSE)))
write.csv(pw, file.path(OUTDIR, "592_detection_rate_by_scenario.csv"), row.names = FALSE)
cat("       FDR<0.05 检出率:\n"); print(pw, row.names = FALSE)

# =============================================================================
# Step 3 · 出图(框架样式;不用条形图)
# =============================================================================
cat("Step 3: 出图\n")
pal <- pal_pub(3, "npg"); names(pal) <- c("dynamic", "null", "static")

# --- Fig1 dumbbell: 静态 -> 拟时序扫描的显著性变化 ---
d1 <- res[order(res$scenario, res$gain_log10p), ]
d1$ord <- seq_len(nrow(d1))
p1 <- ggplot(d1) +
  geom_segment(aes(x = -log10(p_static), xend = -log10(p_scan_bonf),
                   y = ord, yend = ord, colour = scenario), linewidth = 0.5, alpha = 0.75) +
  geom_point(aes(x = -log10(p_static), y = ord), colour = "grey45", size = 1.5) +
  geom_point(aes(x = -log10(p_scan_bonf), y = ord, colour = scenario), size = 2.1) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35, colour = "grey30") +
  scale_colour_manual(values = pal, name = "Simulated scenario") +
  labs(x = expression(-log[10]~italic(P)),
       y = "Gene (ordered by gain)",
       title = "Static TWAS (grey) vs pseudotime-resolved scan (coloured)",
       subtitle = "Baseline comparator; dynamic effects cancel out in a single static weight vector") +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_fig(p1, file.path(ASSETS, "fig1_static_vs_scan_dumbbell"), width = 7.6, height = 6.2)

# --- Fig2 heatmap: 基因 x 拟时序 的 TWAS z ---
long <- data.frame(
  ID = rep(rownames(zmat), times = NPT),
  t  = rep(tgrid, each = nrow(zmat)),
  z  = as.vector(zmat), stringsAsFactors = FALSE)
long$scenario <- res$scenario[match(long$ID, res$ID)]
ordv <- res$ID[order(res$scenario, res$z_peak)]
long$ID <- factor(long$ID, levels = ordv)
p2 <- ggplot(long, aes(x = t, y = ID, fill = z)) +
  geom_tile() +
  facet_grid(scenario ~ ., scales = "free_y", space = "free_y") +
  scale_fill_diverge(midpoint = 0, name = "TWAS z") +
  labs(x = "Pseudotime (scaled 0-1)", y = "Gene",
       title = "Pseudotime-resolved TWAS z-score landscape",
       subtitle = "Sign reversal along pseudotime is what a static test averages away") +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_fig(p2, file.path(ASSETS, "fig2_pseudotime_z_heatmap"), width = 7.2, height = 6.6)

# --- Fig3 轨迹: 各情景代表基因的 z(t) 曲线 ---
reps <- unlist(lapply(split(res, res$scenario), function(d)
  d$ID[order(-abs(d$z_peak))][seq_len(min(4, nrow(d)))]))
d3 <- long[long$ID %in% reps, ]
p3 <- ggplot(d3, aes(x = t, y = z, group = ID, colour = scenario)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 1.1, alpha = 0.7) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +   # 分面边界防止刻度标签相撞
  labs(x = "Pseudotime (scaled 0-1)", y = "TWAS z-score",
       title = "Gene-trait association trajectories across cell states",
       subtitle = "Representative genes per simulated scenario (top |z|)") +
  theme_pub(base_size = 11)
save_fig(p3, file.path(ASSETS, "fig3_effect_trajectories"), width = 8.4, height = 3.4)

# --- Fig4 dot plot: 检出率(不用条形图) ---
d4 <- reshape(pw[, c("scenario", "static", "scan")], direction = "long",
              varying = c("static", "scan"), v.names = "rate",
              timevar = "method", times = c("static", "scan"), idvar = "scenario")
d4$method <- factor(d4$method, levels = c("static", "scan"),
                    labels = c("Static TWAS", "Pseudotime scan"))
p4 <- ggplot(d4, aes(x = rate, y = scenario)) +
  geom_line(aes(group = scenario), colour = "grey65", linewidth = 0.7) +
  geom_point(aes(colour = method), size = 3.6) +
  scale_colour_manual(values = pal_pub(2, "lancet"), name = NULL) +
  scale_x_continuous(limits = c(-0.02, 1.02)) +
  labs(x = "Detection rate (BH-FDR < 0.05)", y = NULL,
       title = "Detection rate by simulated scenario",
       subtitle = "Slopegraph; grid-wise Bonferroni makes the scan the conservative arm") +
  theme_pub(base_size = 11)
save_fig(p4, file.path(ASSETS, "fig4_detection_rate_slopegraph"), width = 6.6, height = 3.2)

# =============================================================================
# Step 4 · 守卫式封装:真正的 TWiST 三联检验
# -----------------------------------------------------------------------------
# 下面的签名逐字取自上游文档,已核实的来源:
#   https://raw.githubusercontent.com/gqi/TWiST/HEAD/man/twist_association.Rd
#   https://raw.githubusercontent.com/gqi/TWiST/HEAD/man/twist_train_model.Rd
#   https://raw.githubusercontent.com/gqi/TWiST/HEAD/example_data/example.R
# 本模块不臆造任何参数;缺包时不降级、不伪造,只打印真实安装命令后退出该分支。
# =============================================================================
twist_real <- function() {
  miss <- c("TWiST", "plink2R", "fda", "grpreg")[
    !vapply(c("TWiST", "plink2R", "fda", "grpreg"),
            function(p) requireNamespace(p, quietly = TRUE), logical(1))]
  if (length(miss)) {
    cat("Step 4: TWiST 正式路径 — 跳过,缺少 R 包:", paste(miss, collapse = ", "), "\n")
    cat("        安装(需联网,本模块不自动安装):\n")
    cat("          install.packages(c('devtools','fda','grpreg','dplyr','data.table'))\n")
    cat("          devtools::install_github('gqi/TWiST')\n")
    cat("          devtools::install_github('gabraham/plink2R/plink2R')\n")
    cat("        预训练权重(CD4+ T / CD8+ T / B 细胞)见 repo 的 pretrained_models/。\n")
    cat("        已核实调用模板(上游 example_data/example.R):\n")
    cat("          load('pretrained_models/twist_weights_T_CD8.rda')  # -> wgtlist, weights_pred, bim_train\n")
    cat("          genos <- plink2R::read_plink('example_data/1000G.EUR.6')\n")
    cat("          ngwas <- ncase*ncontrol/(ncase+ncontrol)\n")
    cat("          res <- TWiST::twist_association(sumstat=, wgtlist=, weights_pred=,\n")
    cat("                                          bim_train=, genos=, ngwas=, opt=)\n")
    cat("          # res$out.tbl: ID CHR P0 P1 tss sigma2 p.global p.dynamic p.nonlinear degree\n")
    return(invisible(list(status = "skipped", missing = miss)))
  }
  cat("Step 4: TWiST 已安装。请按上游 example.R 提供真实 sumstat / 预训练权重 / plink 参考,\n")
  cat("        再调用 TWiST::twist_association();本模块的合成数据不构成合法输入。\n")
  invisible(list(status = "installed"))
}
if (RUN_TWIST) twist_real() else {
  cat("Step 4: 未请求 TWiST 正式路径 (--run-twist);仅运行基线。\n")
}

cat(sprintf("\n[592] 完成。结果 -> %s ;展示图 -> %s\n", OUTDIR, ASSETS))
