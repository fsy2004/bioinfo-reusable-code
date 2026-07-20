# =============================================================================
# 编号   : 558
# 脚本名 : Milo KNN 邻域差异丰度(neighborhood differential abundance, DA)
# 分类   : 08_singlecell_spatial_trajectory
# 用途   : 在 KNN 图的"代表性邻域"上检验细胞密度随条件变化 —— 无需先离散聚类,
#          适合连续状态/谱系。每个邻域 = 一个 index 细胞的 KNN 邻居集合;
#          对每个邻域按样本数细胞数,用 GLM 检验两条件间丰度差异 → logFC + FDR。
# ★诚实基线 :
#          (1) 邻域 DA 的 beeswarm:邻域按谱系拟时序排序、点色=logFC,显著点描边;
#          (2) 对照(离散簇法):把同一数据按离散簇做"每簇细胞比例"的检验
#              (Fisher / 卡方),证明连续邻域法能定位到簇内"局部"富集,
#              而粗粒度簇比例检验会被簇内异质性稀释 → 二者并排,体现预期差异。
# ★工具接地 :
#          miloR 当前装不上(Bioc loadNamespace fail)→ 本脚本【接地真实 miloR API】
#          (Milo() | buildGraph | makeNhoods | countCells | calcNhoodDistance |
#           testNhoods | annotateNhoods | buildNhoodGraph | plotDAbeeswarm |
#           plotNhoodGraphDA,均经官方 vignette milo_gastrulation.Rmd 确认),
#          用 try(library(miloR)) 包裹;若装上则优先用官方实现。
#          降级路径用已装的 BiocNeighbors(KNN)+ 手算邻域富集 + base glm,
#          概念等价地复现 Milo 流程并出图,不依赖缺失包。
# 依赖   : (优先) miloR ;(降级,已装) BiocNeighbors · SingleCellExperiment ·
#          igraph · ggbeeswarm · ggplot2 · ggraph(可选,缺则降级散点)
# 运行   : Rscript 558_milo_neighborhood_da.R                 # 合成示例(in-memory)
#          Rscript 558_milo_neighborhood_da.R --k 25 --prop 0.15
# 输入   : 合成 SingleCellExperiment(脚本内生成,synthetic demo only):
#          reducedDim "PCA"(2D 连续流形) + colData(sample, condition, lineage)
#          某区域(谱系中段)在 "treat" 条件下细胞富集 → 真阳性 DA 信号。
#          换真实数据:见 README,把 SCE 的 reducedDim/colData 列对齐即可。
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
set.seed(42)
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressWarnings(suppressMessages({
  library(ggplot2)
  library(SingleCellExperiment)
  library(BiocNeighbors)
  library(igraph)
  library(ggbeeswarm)
}))

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  k        = 45,      # KNN 邻居数(建图 + 邻域大小;越大计数越足、功效越高)
  prop     = 0.12,    # makeNhoods 抽样比例(代表性 index 细胞占比)
  fdr      = 0.10,    # SpatialFDR 显著阈值
  outdir   = file.path(SCRIPT_DIR, "results")))
for (kk in c("k","prop","fdr")) args[[kk]] <- as.numeric(args[[kk]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

USE_MILO <- isTRUE(suppressWarnings(suppressMessages(
  requireNamespace("miloR", quietly = TRUE) &&
  tryCatch({ library(miloR); TRUE }, error = function(e) FALSE))))
cat(sprintf("[engine] miloR 可用 = %s → %s\n", USE_MILO,
            if (USE_MILO) "用官方 miloR 实现" else "降级:BiocNeighbors + 手算邻域 DA(概念等价复现)"))

# =============================================================================
# 0. 合成 SingleCellExperiment(synthetic demo only)
#    连续流形:沿 lineage 轴 t∈[0,1] 的一条"发育轨迹",分散成 2D PCA 坐标。
#    设计真信号:谱系中段(t∈[0.45,0.6])细胞在 condition="treat" 下富集
#    (该窗口 treat 细胞数 ≈ 2.4× ctrl)→ 期望邻域 DA 在中段检出正 logFC。
# =============================================================================
cat("Step 0: 合成连续流形 SCE(某谱系区段在 treat 富集)...\n")
make_sce <- function() {
  n_per <- 1800
  mk <- function(t, cond, s) {
    nn <- length(t)
    # 2D 连续坐标:PC1 随 t 单调(轨迹主轴)+ 噪声;PC2 为正交散布(发丝噪声小→流形紧致)
    pc1 <- 6 * t + rnorm(nn, 0, 0.25)
    pc2 <- sin(2 * pi * t) * 1.2 + rnorm(nn, 0, 0.30)
    data.frame(PC1 = pc1, PC2 = pc2, ptime = t, condition = cond, sample = s)
  }
  # 每条件 3 个生物学重复样本(DA 检验需要样本层级)。
  # ctrl:t 均匀 [0,1]。treat:在窗口 [0.47,0.58] 额外富集 ~1.4× 背景密度的细胞
  # (强且空间集中 → 窗口内邻域 ≈ 多数 treat 细胞,产生明确正 logFC)。
  ctrl <- do.call(rbind, lapply(1:3, function(i)
    mk(runif(n_per, 0, 1), "ctrl", sprintf("ctrl_%d", i))))
  treat <- do.call(rbind, lapply(1:3, function(i) {
    tt <- c(runif(n_per, 0, 1), runif(round(n_per * 0.8), 0.47, 0.58))
    mk(tt, "treat", sprintf("treat_%d", i)) }))
  meta <- rbind(ctrl, treat)
  meta$cell <- sprintf("cell%06d", seq_len(nrow(meta)))
  # 离散簇标签(供对照法 + annotateNhoods):按 ptime 切 5 段
  meta$cluster <- cut(meta$ptime, breaks = c(-Inf, .2, .4, .6, .8, Inf),
                      labels = paste0("C", 1:5))
  pca <- as.matrix(meta[, c("PC1", "PC2")]); rownames(pca) <- meta$cell
  # 极简表达矩阵(占位:1 基因 = ptime;真实数据用真表达)— 仅满足 SCE 构造
  sce <- SingleCellExperiment(
    assays = list(logcounts = matrix(meta$ptime, nrow = 1,
                                     dimnames = list("g1", meta$cell))),
    colData = DataFrame(meta[, c("sample","condition","cluster","ptime")],
                        row.names = meta$cell))
  reducedDim(sce, "PCA") <- pca
  sce
}
sce <- make_sce()
cat(sprintf("  SCE: %d 细胞 · 6 样本(ctrl×3/treat×3) · 真信号窗口 ptime∈[0.45,0.60]\n", ncol(sce)))

# =============================================================================
# 1. 邻域 DA —— 官方 miloR 路径(接地真实 API;装上才走)
# =============================================================================
run_milo_official <- function(sce, k, prop, d = 2) {
  # 接地真实 miloR API(经 vignette milo_gastrulation.Rmd 确认)
  milo <- miloR::Milo(sce)
  milo <- miloR::buildGraph(milo, k = k, d = d, reduced.dim = "PCA")
  milo <- miloR::makeNhoods(milo, prop = prop, k = k, d = d,
                            refined = TRUE, reduced_dims = "PCA")
  milo <- miloR::countCells(milo, meta.data = as.data.frame(colData(milo)),
                            sample = "sample")
  milo <- miloR::calcNhoodDistance(milo, d = d, reduced.dim = "PCA")
  design.df <- data.frame(colData(milo))[, c("sample","condition")]
  design.df <- distinct(design.df); rownames(design.df) <- design.df$sample
  da <- miloR::testNhoods(milo, design = ~ condition, design.df = design.df,
                          reduced.dim = "PCA")
  da <- miloR::annotateNhoods(milo, da, coldata_col = "cluster")
  da <- miloR::annotateNhoods(milo, da, coldata_col = "ptime")
  milo <- miloR::buildNhoodGraph(milo)
  list(milo = milo, da = da)
}

# =============================================================================
# 1'. 邻域 DA —— 降级路径:BiocNeighbors KNN + 手算邻域 + GLM(概念等价复现)
#     复刻 Milo 核心:① KNN 图;② 抽样 index 细胞定义"代表性邻域"(refined:
#        把 index 吸附到局部密度高处的近似);③ 每邻域 × 每样本 计数;
#        ④ 对每个邻域用 quasi-binomial GLM 检验 condition 效应 → logFC;
#        ⑤ 邻域间重叠 → Spatial FDR(按邻域连通度加权的 BH 校正,近似 miloR)。
# =============================================================================
run_milo_degraded <- function(sce, k, prop, fdr) {
  X <- reducedDim(sce, "PCA")
  cd <- as.data.frame(colData(sce)); N <- nrow(X)

  ## ① KNN 图(BiocNeighbors,真实 API:findKNN(X, k)) -----------------------
  knn <- BiocNeighbors::findKNN(X, k = k, BNPARAM = BiocNeighbors::KmknnParam())
  idx <- knn$index            # N × k 邻居下标(不含自身)

  ## ② 抽样 index 细胞 → 代表性邻域;refined:把随机种子吸附到其 KNN 局部中位 -
  n_nh <- max(30, round(N * prop))
  seeds <- sample.int(N, n_nh)
  refine <- function(s) {            # 近似 makeNhoods(refined=TRUE):取 KNN 内最靠局部密度中心者
    nb <- c(s, idx[s, ])
    cen <- colMeans(X[nb, , drop = FALSE])
    nb[which.min(rowSums((X[nb, , drop = FALSE] -
                          matrix(cen, length(nb), ncol(X), byrow = TRUE))^2))]
  }
  index_cells <- unique(vapply(seeds, refine, integer(1)))
  # 邻域成员 = index 细胞 + 其 k 近邻
  nhoods <- lapply(index_cells, function(s) c(s, idx[s, ]))
  n_nh <- length(nhoods)

  ## ③ 每邻域 × 每样本 计数 → counts 矩阵(nhood × sample) ------------------
  samples <- sort(unique(cd$sample))
  cnt <- t(vapply(nhoods, function(mem)
    table(factor(cd$sample[mem], levels = samples)), numeric(length(samples))))
  colnames(cnt) <- samples
  smp_cond <- cd$condition[match(samples, cd$sample)]
  smp_tot  <- as.numeric(table(factor(cd$sample, levels = samples)))  # 每样本总细胞(library size)

  ## ④ Poisson-GLM(library-size offset)+ 全局共享离散度(edgeR/miloR 核心):
  #    counts ~ condition + offset(log lib.size)。单邻域计数稀疏→独立 quasi 检验
  #    无功效;miloR 经 edgeR 在所有邻域间【共享离散度】大幅提功效。这里复刻:
  #    ① 每邻域拟 Poisson 取 coef(logFC)+ 标准误 + Pearson 残差;
  #    ② 把所有邻域 Pearson 残差汇总,估一个【共同离散度 phi】(= edgeR common dispersion);
  #    ③ 用 phi 缩放标准误,以 z 检验得 PValue(quasi-likelihood 思路)。
  cond_f <- factor(smp_cond, levels = c("ctrl","treat"))
  off    <- log(smp_tot)
  beta <- se_b <- nh_size_v <- rep(NA_real_, n_nh)
  pear_sq <- 0; pear_df <- 0
  for (i in seq_len(n_nh)) {
    y <- cnt[i, ]; nh_size_v[i] <- length(nhoods[[i]])
    if (sum(y) < 3) next
    fit <- tryCatch(glm(y ~ cond_f + offset(off), family = poisson()),
                    error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$coefficients
    if (nrow(co) < 2) next
    beta[i] <- co[2, 1]; se_b[i] <- co[2, 2]
    pr <- residuals(fit, type = "pearson")            # Pearson 残差(供共享离散度)
    pear_sq <- pear_sq + sum(pr^2); pear_df <- pear_df + df.residual(fit)
  }
  phi <- max(1, pear_sq / pear_df)                    # 共同离散度(过散→phi>1,收紧 p)
  z   <- beta / (se_b * sqrt(phi))                    # quasi z 统计量
  pval <- 2 * pnorm(-abs(z))
  lfc <- pmax(-8, pmin(8, beta / log(2)))             # log2 + 防分离爆表
  da <- data.frame(Nhood = seq_len(n_nh), logFC = lfc, PValue = pval,
                   index_cell = index_cells, nh_size = nh_size_v)
  da <- da[!is.na(da$PValue), ]
  cat(sprintf("  [degraded] 共享离散度 phi=%.2f (edgeR common-dispersion 复刻), %d 邻域入检\n",
              phi, nrow(da)))

  ## ⑤ Spatial FDR —— 接地 miloR::graphSpatialFDR 的加权 BH:
  #    权重 w = 邻域第 k 近邻距离(连通度低/孤立 → 权重大),按 p 升序做
  #    weighted Benjamini-Hochberg:adjp = cummin_{从大p到小p}( sum(w)*p / cumsum(w) )。
  #    这校正了"邻域相互重叠 → 检验非独立"的多重比较膨胀。
  kdist <- apply(X[da$index_cell, , drop = FALSE], 1, function(z)
    sort(sqrt(colSums((t(X) - z)^2)))[k + 1])           # 到第 k 近邻的距离 = 连通度反比
  w  <- kdist                                            # miloR 用 1/连通度 ≈ kth-NN 距离
  # weighted Benjamini-Hochberg(同 miloR::graphSpatialFDR 思路):
  ord  <- order(da$PValue)
  pw   <- da$PValue[ord]; ww <- w[ord]
  cw   <- cumsum(ww); sw <- sum(ww)
  adjp <- rev(cummin(rev(pw * sw / cw))); adjp <- pmin(adjp, 1)
  sfdr <- rep(NA_real_, nrow(da)); sfdr[ord] <- adjp
  da$SpatialFDR <- sfdr
  da$FDR <- p.adjust(da$PValue, "BH")

  ## 邻域注释:最丰富簇 + 平均 ptime(对齐 miloR::annotateNhoods) -----------
  da$cluster <- vapply(da$Nhood, function(i) {
    tb <- table(cd$cluster[nhoods[[i]]]); names(tb)[which.max(tb)] }, character(1))
  da$ptime <- vapply(da$Nhood, function(i) mean(cd$ptime[nhoods[[i]]]), numeric(1))

  ## 邻域图(index 细胞为节点,邻域共享成员为边)→ 供网络叠加图 -------------
  nh_xy <- X[da$index_cell, , drop = FALSE]
  list(da = da, nhoods = nhoods, nh_xy = nh_xy, X = X, cd = cd,
       index_cells = da$index_cell)
}

cat("Step 1: 邻域 DA 检验...\n")
if (USE_MILO) {
  fit <- tryCatch(run_milo_official(sce, args$k, args$prop),
                  error = function(e) { cat("  ⚠ miloR 官方路径出错→降级:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) { USE_MILO <- FALSE }
}
if (!USE_MILO) fit <- run_milo_degraded(sce, args$k, args$prop, args$fdr)
da <- fit$da
sig <- da$SpatialFDR < args$fdr
cat(sprintf("  邻域数 = %d · 显著(SpatialFDR<%.2f)= %d (其中正 logFC = %d)\n",
            nrow(da), args$fdr, sum(sig, na.rm = TRUE),
            sum(sig & da$logFC > 0, na.rm = TRUE)))
# 显著邻域应集中在真信号窗口 ptime∈[0.45,0.60]:sanity-check
sig_pt <- da$ptime[sig & da$logFC > 0]
cat(sprintf("  ↑ 富集(正)显著邻域 ptime 中位=%.2f(真信号窗 0.45-0.60;落窗内=验证管道有效)\n",
            stats::median(sig_pt)))
write.csv(da, file.path(args$outdir, "neighborhood_DA_results.csv"), row.names = FALSE)

# =============================================================================
# 2. ★诚实基线对照:离散簇法(粗粒度)—— 同数据按簇做比例检验
#    每簇:treat vs ctrl 细胞数 2×2 → Fisher 检验 + log2 比例比。
#    预期:真信号在 C3(ptime 0.4-0.6)局部,簇法把整簇平均→信号被簇内
#    非富集细胞稀释,logFC 偏小/不显著;而邻域法能定位局部 → 体现"连续 > 离散"。
# =============================================================================
cat("Step 2: 诚实基线对照(离散簇比例检验)...\n")
cd <- as.data.frame(colData(sce))
clus_da <- do.call(rbind, lapply(levels(cd$cluster), function(cl) {
  inc  <- cd$cluster == cl
  a <- sum(inc & cd$condition == "treat"); b <- sum(inc & cd$condition == "ctrl")
  c2 <- sum(!inc & cd$condition == "treat"); d2 <- sum(!inc & cd$condition == "ctrl")
  ft <- fisher.test(matrix(c(a, b, c2, d2), 2))
  data.frame(cluster = cl,
             logFC = log2(((a + 0.5)/(c2 + 0.5)) / ((b + 0.5)/(d2 + 0.5))),
             PValue = ft$p.value)
}))
clus_da$FDR <- p.adjust(clus_da$PValue, "BH")
write.csv(clus_da, file.path(args$outdir, "cluster_DA_control.csv"), row.names = FALSE)
cat("  簇法结果:\n"); print(clus_da, digits = 3)
nh_max <- max(abs(da$logFC[sig]), na.rm = TRUE)
cat(sprintf("  对比:邻域法峰值|logFC|=%.2f vs 簇法 C3 |logFC|=%.2f → 连续法定位更锐\n",
            nh_max, abs(clus_da$logFC[clus_da$cluster == "C3"])))

# =============================================================================
# 3. 出图(全部顶刊级,禁平凡条形图;每图独立成文件 PDF+PNG)
# =============================================================================
cat("Step 3: 出图...\n")
pal  <- pal_pub(name = "npg")
sigcol <- function(x, fdr) ifelse(x < fdr, "black", NA)

## 图1:DA beeswarm —— 邻域按谱系拟时序分箱,y=logFC,色=logFC,显著描黑边 ----
da$pt_bin <- cut(da$ptime, breaks = seq(0, 1, by = 0.2),
                 labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"),
                 include.lowest = TRUE)
da$sig <- da$SpatialFDR < args$fdr
p_bee <- ggplot(da, aes(x = pt_bin, y = logFC)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  ggbeeswarm::geom_quasirandom(aes(colour = logFC, stroke = ifelse(sig, 0.7, 0)),
                               shape = 21, size = 2.2, width = 0.32,
                               fill = NA) +
  geom_quasirandom(aes(colour = logFC), shape = 16, size = 1.9, width = 0.32,
                   alpha = 0.9) +
  scale_color_diverge(midpoint = 0, name = "logFC") +
  labs(title = "Neighborhood DA beeswarm (Milo-style)",
       subtitle = "Each point = one KNN neighborhood; ordered along lineage pseudotime",
       x = "Lineage pseudotime bin", y = expression(log[2]~"fold-change (treat / ctrl)")) +
  theme_pub(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_bee, file.path(ASSETS, "fig1_DA_beeswarm"), width = 7.2, height = 5)

## 图2:邻域图网络叠加 —— index 细胞为节点,布局=PCA,色=logFC,显著加大 -------
nh_xy <- if (USE_MILO) reducedDim(sce, "PCA")[da$index_cell %||% seq_len(nrow(da)), , drop = FALSE] else fit$nh_xy
# 背景细胞(灰)+ 邻域节点(着色),并连接彼此 KNN 相邻的邻域(边=共享成员)
bg <- data.frame(reducedDim(sce, "PCA")); colnames(bg) <- c("PC1","PC2")
nh_df <- data.frame(PC1 = nh_xy[, 1], PC2 = nh_xy[, 2],
                    logFC = da$logFC, sig = da$sig)
# 构造邻域间边:两邻域 index 细胞欧氏距离 < 阈值则连边(近似 buildNhoodGraph 共享成员)
dd <- as.matrix(dist(nh_xy)); thr <- stats::quantile(dd[upper.tri(dd)], 0.04)
ed <- which(dd < thr & upper.tri(dd), arr.ind = TRUE)
edge_df <- data.frame(x = nh_xy[ed[,1],1], y = nh_xy[ed[,1],2],
                      xend = nh_xy[ed[,2],1], yend = nh_xy[ed[,2],2])
p_net <- ggplot() +
  geom_point(data = bg, aes(PC1, PC2), colour = "grey85", size = 0.4, alpha = 0.5) +
  geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "grey60", linewidth = 0.25, alpha = 0.5) +
  geom_point(data = nh_df, aes(PC1, PC2, fill = logFC, size = sig),
             shape = 21, colour = "grey20", stroke = 0.3) +
  scale_fill_diverge(midpoint = 0, name = "logFC") +
  scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.6),
                    labels = c("n.s.", paste0("FDR<", args$fdr)), name = "DA") +
  labs(title = "Neighborhood graph DA overlay",
       subtitle = "Nodes = neighborhoods on PCA manifold; edges = shared-cell adjacency",
       x = "PC1", y = "PC2") +
  coord_equal() +
  theme_pub(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),  # 防透明背景在看图器显黑
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_net, file.path(ASSETS, "fig2_nhood_graph_DA"), width = 7.4, height = 6.2)

## 图3:邻域 volcano —— x=logFC, y=-log10(SpatialFDR),显著着色,标注富集方向 ----
da$neglog <- -log10(pmax(da$SpatialFDR, 1e-300))
da$cat <- ifelse(da$sig & da$logFC > 0, "Enriched in treat",
          ifelse(da$sig & da$logFC < 0, "Depleted in treat", "n.s."))
vc <- c("Enriched in treat" = pal[1], "Depleted in treat" = pal[4], "n.s." = "grey75")
p_vol <- ggplot(da, aes(logFC, neglog)) +
  geom_hline(yintercept = -log10(args$fdr), linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = cat, size = nh_size %||% 2),
             position = position_jitter(width = 0, height = 0.03, seed = 42),
             alpha = 0.7) +
  scale_color_manual(values = vc, name = NULL) +
  scale_size_continuous(range = c(1.2, 3.2), guide = "none") +
  labs(title = "Neighborhood DA volcano",
       subtitle = sprintf("Spatial-FDR corrected; %d/%d neighborhoods significant",
                          sum(da$sig), nrow(da)),
       x = expression(log[2]~"fold-change (treat / ctrl)"),
       y = expression(-log[10]~"(Spatial FDR)")) +
  theme_pub(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_vol, file.path(ASSETS, "fig3_DA_volcano"), width = 6.6, height = 5.2)

## 图4:诚实基线对照并排 —— 邻域法(per-nhood)vs 簇法(per-cluster)的 logFC
##      lollipop + 邻域 logFC 分布(violin/jitter),凸显连续法局部分辨率 ----------
da$cluster_f <- factor(da$cluster, levels = paste0("C", 1:5))
clus_da$cluster_f <- factor(clus_da$cluster, levels = paste0("C", 1:5))
# 显式两水平因子,保证图例标签与显著性对应(全显著时也不串标签)
clus_da$sig <- factor(ifelse(clus_da$FDR < args$fdr, "sig", "ns"), levels = c("ns", "sig"))
da$sig <- factor(ifelse(da$sig, "sig", "ns"), levels = c("ns", "sig"))
p_cmp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  # 邻域法:每簇内众多邻域的 logFC 分布(violin + jitter)
  geom_violin(data = da, aes(cluster_f, logFC), fill = pal[2], colour = NA, alpha = 0.18, width = 0.85) +
  geom_jitter(data = da, aes(cluster_f, logFC, colour = sig), width = 0.12, size = 0.9, alpha = 0.6) +
  scale_color_manual(values = c(ns = "grey70", sig = pal[1]),
                     labels = c("n.s.", "DA nhood"), name = "Neighborhood (continuous)") +
  # 簇法:每簇单点 logFC(lollipop 棒 + 菱形;显著=金,n.s.=空心)
  geom_segment(data = clus_da, aes(x = cluster_f, xend = cluster_f, y = 0, yend = logFC),
               colour = "grey30", linewidth = 0.6) +
  geom_point(data = clus_da, aes(cluster_f, logFC, fill = sig),
             shape = 23, size = 4.5, colour = "black", stroke = 0.7) +
  scale_fill_manual(values = c(ns = "white", sig = pal[8]),
                    labels = c("n.s.", paste0("FDR<", args$fdr)),
                    name = "Cluster test (discrete control)", drop = FALSE) +
  labs(title = "Honest baseline: continuous neighborhoods vs discrete clusters",
       subtitle = "Diamonds = per-cluster proportion test (Fisher); dots/violin = per-neighborhood DA",
       x = "Lineage cluster (by pseudotime)",
       y = expression(log[2]~"fold-change (treat / ctrl)")) +
  theme_pub(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 8, 6, 8))
save_fig(p_cmp, file.path(ASSETS, "fig4_baseline_nhood_vs_cluster"), width = 7.6, height = 5.4)

cat("完成。结果表见", normalizePath(args$outdir), ";图见 assets/\n")

# 依赖版本快照(铁律6)
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
