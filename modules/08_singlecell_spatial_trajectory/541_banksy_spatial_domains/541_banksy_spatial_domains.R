# =============================================================================
# 编号   : 541
# 脚本名 : BANKSY 空间域识别(spatial domain segmentation)+ 诚实非空间基线对照
# 分类   : 08_singlecell_spatial_trajectory
# 用途   : 用 BANKSY(Bioconductor)做空间转录组的"空间域"识别:
#          把每个 spot 的【自身表达 + 邻域均值(harmonic m=0)+ 方位 Gabor 梯度
#          (azimuthal Gabor, harmonic m=1)】拼成 neighbor-augmented 特征,再降维聚类。
#          非深度学习、完全可解释:lambda 旋钮显式控制"邻域上下文"的权重。
# ★诚实基线 : BANKSY 的 lambda 旋钮天然给出对照——
#             lambda=0  → 退化为【普通非空间聚类】(仅用自身表达,等价 Leiden/kmeans on expression);
#             lambda=0.8 → 【空间域模式】(强邻域增强)。
#             二者对同一合成数据跑,用 ARI(vs 已知真域)量化"空间增强带来的域连贯性提升"。
#             不只报好看指标:非空间基线 ARI 会被一并打印/出图,差距即增益。
# 依赖   : Banksy · SpatialExperiment · SummarizedExperiment · S4Vectors(Bioc)
#          aricode 或 mclust(算 ARI,二选一即可)· ggplot2(+ theme_pub.R)
# 运行   : Rscript 541_banksy_spatial_domains.R                      # 合成示例(自动生成)
#          Rscript 541_banksy_spatial_domains.R --input data/你的.csv --outdir results/run1
# 输入   : 一张 long-format 空间表达 CSV(synthetic demo 自动生成,见 ① README):
#          列 = spot, x, y, domain(可选,真值,仅评估用), gene1, gene2, ... geneN
#          即每行一个 spot:坐标 + 各基因表达;有 domain 列则算 ARI,无则只出分割图。
# =============================================================================

## ---- 0. 框架 + 依赖 + 参数 -------------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(Banksy); library(SpatialExperiment); library(SummarizedExperiment)
  library(S4Vectors); library(ggplot2)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input   = file.path(DDAT, "spatial_demo.csv"),
  outdir  = file.path(SCRIPT_DIR, "results"),
  k_geom  = 18,        # 邻域 kNN 大小(m=0 用 k_geom,m=1 用 2*k_geom)
  lambda  = 0.8,       # 空间域模式权重(0=纯表达基线,0.8=域分割;BANKSY 文档推荐)
  res     = 0.8,       # Leiden 分辨率
  npcs    = 20))
for (k in c("k_geom", "lambda", "res", "npcs")) args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ARI 工具(aricode 优先,退而 mclust;均无则给 NA 不报错) ----------------------
ari_fun <- function(a, b) {
  if (requireNamespace("aricode", quietly = TRUE)) return(aricode::ARI(a, b))
  if (requireNamespace("mclust",  quietly = TRUE)) return(mclust::adjustedRandIndex(a, b))
  NA_real_
}

## ---- 1. 合成空间转录组(synthetic demo only)------------------------------
# 设计:40x40 spot 网格 = 4 个连续方块域(象限);每域 10 个 marker 基因均值抬高,
# 但单 spot 表达噪声很大(sd=3.5)→ 仅看自身表达难以分辨域(非空间基线会差);
# BANKSY 用"邻域均值 + 方位梯度"借空间上下文把噪声平掉 → 域连贯性大幅提升。
# 此设计刻意让空间增强【有用】(域内表达相似但噪声大),以诚实展示 ARI 提升,而非掩盖差距。
if (!file.exists(args$input)) {
  cat("[gen] 生成合成空间转录组 (synthetic, for demo only)...\n")
  side <- 40L
  gx <- rep(seq_len(side), times = side)
  gy <- rep(seq_len(side), each  = side)
  # 4 象限域(1=左下,2=左上,3=右下,4=右上)
  dom <- ifelse(gx <= side / 2, ifelse(gy <= side / 2, 1L, 2L),
                                 ifelse(gy <= side / 2, 3L, 4L))
  ng <- 40L
  mu <- matrix(2, nrow = ng, ncol = 4)                     # 全基因基线均值 2
  for (d in 1:4) mu[((d - 1) * 10 + 1):(d * 10), d] <- 4   # 各域 10 个 marker 抬到 4
  noise_sd <- 3.5
  expr <- vapply(seq_along(dom),
                 function(i) rnorm(ng, mu[, dom[i]], noise_sd), numeric(ng))
  expr[expr < 0] <- 0
  rownames(expr) <- sprintf("gene%02d", seq_len(ng))
  df <- data.frame(spot = sprintf("spot%04d", seq_along(dom)),
                   x = gx, y = gy, domain = dom,
                   t(expr), check.names = FALSE)
  write.csv(df, args$input, row.names = FALSE)
  cat(sprintf("  写入 %s :%d spots × %d genes,4 个方块域,噪声 sd=%.1f\n",
              basename(args$input), length(dom), ng, noise_sd))
}

## ---- 2. 读表 → 组装 SpatialExperiment -------------------------------------
cat("Step 1: 读 long-format CSV → 构建 SpatialExperiment...\n")
raw <- read.csv(args$input, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("x", "y") %in% colnames(raw)))
has_truth <- "domain" %in% colnames(raw)
spot_id <- if ("spot" %in% colnames(raw)) raw$spot else sprintf("spot%05d", seq_len(nrow(raw)))
meta_cols <- intersect(c("spot", "x", "y", "domain"), colnames(raw))
gene_cols <- setdiff(colnames(raw), meta_cols)
expr_mat  <- t(as.matrix(raw[, gene_cols, drop = FALSE]))      # 基因 × spot
colnames(expr_mat) <- spot_id; mode(expr_mat) <- "numeric"
coords <- as.matrix(raw[, c("x", "y")]); rownames(coords) <- spot_id

cd <- S4Vectors::DataFrame(row.names = spot_id)
if (has_truth) cd$domain <- factor(raw$domain)
se <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = expr_mat), spatialCoords = coords, colData = cd)
cat(sprintf("  SpatialExperiment: %d genes × %d spots%s\n",
            nrow(se), ncol(se), if (has_truth) " (含 domain 真值,将算 ARI)" else ""))

## ---- 3. 计算 BANKSY 邻域增强特征(自身 + m=0 邻域均值 + m=1 方位 Gabor)-----
cat("Step 2: computeBanksy — 邻域均值(H0)+ 方位 Gabor 梯度(H1)...\n")
se <- computeBanksy(se, assay_name = "counts", compute_agf = TRUE,
                    k_geom = c(args$k_geom, args$k_geom * 2), verbose = FALSE)
cat(sprintf("  新增 neighbour 矩阵 assays: %s\n",
            paste(setdiff(assayNames(se), "counts"), collapse = ", ")))

## ---- 4. PCA:同时算【非空间基线 lambda=0】与【空间域 lambda=0.8】两套嵌入 ----
cat("Step 3: runBanksyPCA — 同时算 lambda=0(非空间基线)与 lambda=", args$lambda, "(空间域)...\n", sep = "")
lams <- sort(unique(c(0, args$lambda)))     # 0 = 诚实非空间基线;args$lambda = 空间域
se <- runBanksyPCA(se, use_agf = TRUE, lambda = lams, npcs = args$npcs, seed = 42)
cat(sprintf("  reducedDims: %s\n", paste(reducedDimNames(se), collapse = ", ")))

## ---- 5. 聚类:两套 lambda 各跑 Leiden;并显式提取"基线 vs 空间"两列标签 -------
cat("Step 4: clusterBanksy — Leiden 聚类(基线 lambda=0 与 空间 lambda=", args$lambda, ")...\n", sep = "")
se <- clusterBanksy(se, use_agf = TRUE, lambda = lams, resolution = args$res,
                    algo = "leiden", seed = 42)
cn <- clusterNames(se)
base_name <- grep("lam0_",  cn, value = TRUE)[1]          # lambda=0 → 非空间基线
if (is.na(base_name)) base_name <- grep("lam0\\b", cn, value = TRUE)[1]
spat_name <- setdiff(cn, base_name)[1]                    # lambda>0 → 空间域
# 与真值对齐编号(便于配色一致;不改聚类内容,仅重映射标签数字)
if (has_truth) {
  for (c in c(base_name, spat_name))
    colData(se)[[c]] <- factor(connectClusters(se, verbose = FALSE)[[c]])
}
lab_base <- colData(se)[[base_name]]
lab_spat <- colData(se)[[spat_name]]
cat(sprintf("  非空间基线 clusters=%d · 空间域 clusters=%d\n",
            length(unique(lab_base)), length(unique(lab_spat))))

## ---- 6. ★诚实基线评估:ARI(域连贯性)非空间 vs 空间 ------------------------
cat("Step 5: ★诚实基线对照 — ARI(vs 已知真域)...\n")
eval_tab <- data.frame(
  method = c("Non-spatial baseline (lambda=0)",
             sprintf("BANKSY spatial (lambda=%.2g)", args$lambda)),
  lambda = c(0, args$lambda),
  n_clusters = c(length(unique(lab_base)), length(unique(lab_spat))),
  stringsAsFactors = FALSE)
if (has_truth) {
  truth <- colData(se)$domain
  eval_tab$ARI <- c(ari_fun(lab_base, truth), ari_fun(lab_spat, truth))
  cat(sprintf("  非空间基线  ARI = %.3f\n", eval_tab$ARI[1]))
  cat(sprintf("  BANKSY 空间 ARI = %.3f  →  增益 ΔARI = %+.3f\n",
              eval_tab$ARI[2], eval_tab$ARI[2] - eval_tab$ARI[1]))
} else {
  eval_tab$ARI <- NA_real_
  cat("  (输入无 domain 真值列 → 跳过 ARI;仅出分割图)\n")
}
write.csv(eval_tab, file.path(args$outdir, "baseline_vs_banksy_ARI.csv"), row.names = FALSE)

## ---- 7. UMAP(空间嵌入 + 非空间嵌入,用于特征对比图)------------------------
cat("Step 6: runBanksyUMAP — 空间 & 非空间嵌入的 UMAP...\n")
se <- tryCatch(
  runBanksyUMAP(se, use_agf = TRUE, lambda = lams, npcs = args$npcs, seed = 42),
  error = function(e) { cat("  ⚠ UMAP 失败(", conditionMessage(e), "),跳过 UMAP 图\n"); se })
umap_names <- grep("UMAP", reducedDimNames(se), value = TRUE)

## ---- 8. 顶刊级出图 ---------------------------------------------------------
cat("Step 7: 出图(空间分割 / ARI 对比 / 邻域增强 UMAP)...\n")

pdom <- if (has_truth) length(levels(droplevels(colData(se)$domain))) else 4
pal_dom <- pal_pub(max(8, pdom + 1), "npg")

spot_df <- data.frame(spatialCoords(se),
                      baseline = lab_base, spatial = lab_spat,
                      check.names = FALSE)
colnames(spot_df)[1:2] <- c("x", "y")
if (has_truth) spot_df$truth <- colData(se)$domain

# ---- 图1:空间域分割图(spot 按域上色;真值 / 非空间基线 / BANKSY 三联)---------
mk_spatial <- function(df, fillcol, title, sub) {
  ggplot(df, aes(x = x, y = y, fill = .data[[fillcol]])) +
    geom_tile(width = 1, height = 1) +
    coord_equal() +
    scale_fill_manual(values = pal_dom, name = "Domain") +
    labs(title = title, subtitle = sub, x = "Spatial X", y = "Spatial Y") +
    theme_pub(base_size = 12) +
    theme(panel.grid = element_blank())
}
panels <- list()
if (has_truth) panels$truth <- mk_spatial(spot_df, "truth", "Ground-truth domains", "Synthetic spatial layout")
panels$base <- mk_spatial(spot_df, "baseline", "Non-spatial baseline",
  if (has_truth) sprintf("Leiden on expression only (lambda=0) · ARI=%.3f", eval_tab$ARI[1]) else "Leiden on expression only (lambda=0)")
panels$spat <- mk_spatial(spot_df, "spatial", "BANKSY spatial domains",
  if (has_truth) sprintf("Neighbour-augmented (lambda=%.2g) · ARI=%.3f", args$lambda, eval_tab$ARI[2])
  else sprintf("Neighbour-augmented (lambda=%.2g)", args$lambda))
# 单独存每张
save_fig(panels$spat, file.path(ASSETS, "fig1_spatial_domains_banksy"), width = 6, height = 5.4)
save_fig(panels$base, file.path(ASSETS, "fig1b_spatial_domains_baseline"), width = 6, height = 5.4)
# 三联对比(可选合成)
trio <- compose_panels(panels[intersect(c("truth", "base", "spat"), names(panels))],
                       ncol = length(panels), tag = "A")
save_fig(trio, file.path(ASSETS, "fig2_segmentation_compare"),
         width = 5.4 * length(panels), height = 5.4)

# ---- 图3:ARI 对比 — lollipop(禁用条形图)----------------------------------
if (has_truth) {
  ari_df <- eval_tab
  ari_df$method <- factor(ari_df$method, levels = ari_df$method[order(ari_df$ARI)])
  p_ari <- ggplot(ari_df, aes(x = ARI, y = method, colour = method)) +
    geom_segment(aes(x = 0, xend = ARI, yend = method), linewidth = 1.4, colour = "grey70") +
    geom_point(size = 6) +
    geom_text(aes(label = sprintf("%.3f", ARI)), colour = "black", size = 3.6,
              fontface = "bold", hjust = -0.45) +
    scale_colour_manual(values = pal_pub(2, "lancet"), guide = "none") +
    scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Spatial augmentation improves domain coherence",
         subtitle = sprintf("Adjusted Rand Index vs ground truth · gain = %+.3f",
                            ari_df$ARI[which.max(ari_df$ARI)] - ari_df$ARI[which.min(ari_df$ARI)]),
         x = "Adjusted Rand Index (ARI)", y = NULL) +
    theme_pub(base_size = 12)
  save_fig(p_ari, file.path(ASSETS, "fig3_ari_baseline_vs_banksy"), width = 7.5, height = 3.6)
}

# ---- 图4:邻域增强特征 UMAP(空间嵌入,按 BANKSY 域上色;若有真值再叠真值形状)----
if (length(umap_names)) {
  un_spat <- grep(sprintf("lam%s", sub("0\\.", "0_", as.character(args$lambda))), umap_names, value = TRUE)
  if (!length(un_spat)) un_spat <- grep("lam0\\.", umap_names, value = TRUE)
  if (!length(un_spat)) un_spat <- setdiff(umap_names, grep("lam0_", umap_names, value = TRUE))
  un_spat <- un_spat[1]
  um <- as.data.frame(reducedDim(se, un_spat)); colnames(um)[1:2] <- c("UMAP1", "UMAP2")
  um$spatial <- lab_spat
  p_um <- ggplot(um, aes(UMAP1, UMAP2, colour = spatial)) +
    geom_point(size = 1.1, alpha = 0.85) +
    scale_colour_manual(values = pal_dom, name = "BANKSY\ndomain") +
    labs(title = "Neighbour-augmented feature UMAP",
         subtitle = sprintf("BANKSY embedding (self + mean + AGF, lambda=%.2g)", args$lambda)) +
    theme_pub(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  save_fig(p_um, file.path(ASSETS, "fig4_banksy_feature_umap"), width = 6, height = 5.2)
}

## ---- 9. 落盘 + 依赖快照 ----------------------------------------------------
write.csv(spot_df, file.path(args$outdir, "spot_domain_assignments.csv"), row.names = FALSE)
cat("完成。结果表见", normalizePath(args$outdir), ";展示图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖版本快照(铁律6)
