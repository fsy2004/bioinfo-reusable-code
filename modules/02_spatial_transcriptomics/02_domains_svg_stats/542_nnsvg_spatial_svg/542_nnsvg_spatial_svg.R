# =============================================================================
# 编号   : 542
# 脚本名 : nnSVG 空间可变基因 (Spatially Variable Genes, SVG) 识别 + 诚实 HVG 基线
# 分类   : 08_singlecell_spatial_trajectory
# 用途   : 在空间转录组 (SpatialExperiment) 上用 nnSVG (最近邻高斯过程, NNGP) 线性可
#          扩展地识别"空间可变基因"——表达随空间坐标呈梯度/斑块的基因,而非单纯高表达。
# ★诚实基线: 同时跑【非空间 HVG】(scran modelGeneVar / getTopHVGs,完全不看坐标)作对照。
#            合成数据里我们埋了两类阳性基因:① 真有空间结构(梯度/斑块) ② 仅高变无结构。
#            正确结果应是: nnSVG 把"空间结构"基因排到前面,而 HVG 会把"仅高变"基因也抬高
#            → SVG vs HVG 散点能直观看出"谁抓空间、谁抓方差",证明 SVG 抓的是空间而非高变。
# 依赖   : nnSVG · SpatialExperiment · SingleCellExperiment · scran · scater · scuttle (Bioc)
#          + ggplot2 (framework theme_pub.R)
# 运行   : Rscript 542_nnsvg_spatial_svg.R                       # 零改动跑合成示例
#          Rscript 542_nnsvg_spatial_svg.R --counts my_counts.csv --coords my_coords.csv
# 输入   : counts = 基因 x spot 计数矩阵 (行=基因, 列=spot; 首列基因名), csv
#          coords = spot 坐标表 (列: spot,x,y), csv
#          (合成示例由脚本内生成, 标注 synthetic demo only)
# 备注   : Castl(集成多 SVG 方法)/ SPARK-X 需另装,本机未装,故本模块仅用 nnSVG;
#          README 注明扩展方式。nnSVG 结果落在 rowData(spe)$LR_stat / rank / pval / padj。
# =============================================================================

## ---- 0. 载入框架主题 + 真实工具包 ------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(SpatialExperiment); library(SingleCellExperiment)
  library(nnSVG); library(scran); library(scater); library(scuttle)
  library(ggplot2)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  counts   = file.path(DDAT, "spatial_counts.csv"),
  coords   = file.path(DDAT, "spatial_coords.csv"),
  outdir   = file.path(SCRIPT_DIR, "results"),
  n_top    = 6,      # 空间表达图展示的 top SVG 数
  n_hvg    = 30,     # HVG 基线取前多少个高变基因
  n_threads = 1))
for (k in c("n_top","n_hvg","n_threads")) args[[k]] <- as.integer(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. 合成空间转录组 (synthetic demo only) -------------------------------
# 设计 (埋下"地面真值"以便检验 SVG vs HVG):
#   · 80 spot x 80 spot 网格? 用随机散点更像真实组织: n_spot 个 spot 随机分布于单位方形。
#   · gradient_*  : 沿 x 方向线性梯度 (强空间结构, 中等方差)        → nnSVG 应排最前
#   · patch_*     : 局部圆形斑块高表达 (强空间结构, 块状)            → nnSVG 应排最前
#   · hvgonly_*   : 高方差但表达与坐标【无关】(随机重排坐标标签)    → HVG 抬高, nnSVG 不该高
#   · noise_*     : 低方差泊松噪声 (既不空间也不高变)                → 两者都低
if (!(file.exists(args$counts) && file.exists(args$coords))) {
  n_spot <- 600
  cx <- runif(n_spot); cy <- runif(n_spot)
  mk <- function(lambda_vec) t(sapply(seq_len(length(lambda_vec)),
                                      function(i) rpois(1, 0)))  # placeholder, replaced below
  rows <- list(); meta <- character(0)
  add_gene <- function(name, lambda) { rows[[name]] <<- rpois(n_spot, lambda); }
  # 真空间-梯度 (6 个): 表达 = 基线 + 斜率 * x
  for (g in 1:6) { add_gene(sprintf("gradient_%02d", g), 2 + 16 * cx) }
  # 真空间-斑块 (6 个): 圆心随机, 圆内高表达, 圆外低
  for (g in 1:6) {
    ox <- runif(1, 0.25, 0.75); oy <- runif(1, 0.25, 0.75); r <- 0.22
    inside <- (cx - ox)^2 + (cy - oy)^2 < r^2
    add_gene(sprintf("patch_%02d", g), ifelse(inside, 18, 2))
  }
  # 仅高变、无空间 (8 个): 用一个有结构的 lambda 但把坐标-表达对应关系打乱 → 高方差但空间随机
  for (g in 1:8) {
    base_lambda <- 2 + 16 * cx                 # 本来有梯度
    perm <- sample(n_spot)                      # 打乱 → 空间随机
    add_gene(sprintf("hvgonly_%02d", g), base_lambda[perm])
  }
  # 噪声基因 (30 个): 低 lambda 泊松, 无结构无高变
  for (g in 1:30) add_gene(sprintf("noise_%02d", g), 4)

  counts <- do.call(rbind, rows)
  rownames(counts) <- names(rows)
  colnames(counts) <- sprintf("spot%04d", seq_len(n_spot))
  storage.mode(counts) <- "integer"
  write.csv(data.frame(gene = rownames(counts), counts, check.names = FALSE),
            args$counts, row.names = FALSE)
  write.csv(data.frame(spot = colnames(counts), x = cx, y = cy),
            args$coords, row.names = FALSE)
  cat(sprintf("[gen] 合成空间数据: %d 基因 x %d spot (6 梯度 + 6 斑块 + 8 仅高变 + 30 噪声)\n",
              nrow(counts), ncol(counts)))
}

## ---- 2. 读入 → 构建 SpatialExperiment ------------------------------------
cat("Step 1: 读 counts + coords → SpatialExperiment...\n")
cdf <- read.csv(args$counts, check.names = FALSE)
genes <- cdf[[1]]; cm <- as.matrix(cdf[, -1, drop = FALSE]); rownames(cm) <- genes
storage.mode(cm) <- "integer"
coords_df <- read.csv(args$coords)
stopifnot(all(coords_df$spot == colnames(cm)))   # spot 对齐 sanity-check
coords_mat <- as.matrix(coords_df[, c("x","y")]); rownames(coords_mat) <- coords_df$spot

spe <- SpatialExperiment(assays = list(counts = cm), spatialCoords = coords_mat)
spe$sample_id <- "demo"
rowData(spe)$gene_name <- rownames(spe)          # filter_genes 要求此列
cat(sprintf("  SPE: %d genes x %d spots\n", nrow(spe), ncol(spe)))

## ---- 3. 标准化 (logcounts) + nnSVG 基因过滤 -------------------------------
cat("Step 2: logNormCounts + filter_genes...\n")
spe <- logNormCounts(spe)
# 合成数据基因少、表达充分 → 放宽过滤阈值以保留全部信号基因 (真实数据用默认 3 / 0.5)
spe <- filter_genes(spe, filter_genes_ncounts = 2, filter_genes_pcspots = 0.2, filter_mito = FALSE)
cat(sprintf("  过滤后保留 %d 基因\n", nrow(spe)))

## ---- 4. nnSVG: 最近邻高斯过程识别空间可变基因 ------------------------------
cat("Step 3: nnSVG (NNGP, 线性可扩展)...\n")
spe <- nnSVG(spe, n_threads = args$n_threads)
svg <- as.data.frame(rowData(spe))
svg$gene <- rownames(spe)
svg <- svg[order(svg$rank), ]
write.csv(svg, file.path(args$outdir, "nnSVG_results.csv"), row.names = FALSE)
n_sig <- sum(svg$padj < 0.05, na.rm = TRUE)
cat(sprintf("  nnSVG: %d / %d 基因 padj<0.05 显著空间可变; top1 = %s (LR=%.1f)\n",
            n_sig, nrow(svg), svg$gene[1], svg$LR_stat[1]))

## ---- 5. ★诚实基线: 非空间 HVG (scran, 完全不看坐标) -----------------------
cat("Step 4: ★诚实基线 — 非空间 HVG (modelGeneVar, 忽略坐标)...\n")
dec <- modelGeneVar(spe)                          # 仅按表达方差建模, 无任何空间信息
hvg_df <- as.data.frame(dec)
hvg_df$gene <- rownames(dec)
hvg_df <- hvg_df[order(-hvg_df$bio), ]            # bio = 生物学方差成分, 越大越"高变"
hvg_df$hvg_rank <- seq_len(nrow(hvg_df))
top_hvg <- getTopHVGs(dec, n = min(args$n_hvg, nrow(dec)))
write.csv(hvg_df, file.path(args$outdir, "HVG_baseline.csv"), row.names = FALSE)

# 合并 SVG 与 HVG 排名, 标注地面真值类别
cls <- function(g) ifelse(grepl("^gradient", g), "spatial-gradient",
                   ifelse(grepl("^patch", g),    "spatial-patch",
                   ifelse(grepl("^hvgonly", g),  "HVG-only (no space)", "noise")))
cmp <- merge(
  svg[, c("gene","rank","LR_stat","padj")],
  hvg_df[, c("gene","hvg_rank","bio")], by = "gene")
cmp$class <- factor(cls(cmp$gene),
  levels = c("spatial-gradient","spatial-patch","HVG-only (no space)","noise"))
cmp$svg_sig <- cmp$padj < 0.05
write.csv(cmp, file.path(args$outdir, "SVG_vs_HVG_comparison.csv"), row.names = FALSE)

# 诚实基线实测小结: SVG top-N 与 HVG top-N 的类别构成
topN <- args$n_hvg
svg_top  <- cmp$class[order(cmp$rank)][seq_len(min(topN, nrow(cmp)))]
hvg_top  <- cmp$class[order(cmp$hvg_rank)][seq_len(min(topN, nrow(cmp)))]
cat(sprintf("  [基线对照] SVG-top%d 里 spatial=%d, HVG-only=%d, noise=%d\n", topN,
            sum(grepl("spatial", svg_top)), sum(svg_top == "HVG-only (no space)"),
            sum(svg_top == "noise")))
cat(sprintf("            HVG-top%d 里 spatial=%d, HVG-only=%d, noise=%d  ← HVG 误收 HVG-only\n", topN,
            sum(grepl("spatial", hvg_top)), sum(hvg_top == "HVG-only (no space)"),
            sum(hvg_top == "noise")))

## ===========================================================================
## 出图 (全部顶刊风格, 禁止平凡条形图; 每图独立成文件 PDF+PNG)
## ===========================================================================
coord_plot_df <- as.data.frame(spatialCoords(spe))
colnames(coord_plot_df) <- c("x","y")
logc <- as.matrix(logcounts(spe))

## ---- 图1: top SVG 空间表达图 (表达叠坐标, 多 panel facet, viridis) --------
cat("Step 5: 图1 top-SVG 空间表达 facet...\n")
top_genes <- svg$gene[seq_len(min(args$n_top, nrow(svg)))]
expr_long <- do.call(rbind, lapply(top_genes, function(g) {
  data.frame(coord_plot_df, expr = logc[g, ],
             gene = sprintf("%s (rank %d)", g, svg$rank[svg$gene == g]))
}))
expr_long$gene <- factor(expr_long$gene, levels = unique(expr_long$gene))
p1 <- ggplot(expr_long, aes(x, y, colour = expr)) +
  geom_point(size = 0.9) +
  facet_wrap(~ gene, ncol = 3) +
  scale_color_cont(option = "D", name = "logcounts") +
  coord_equal() +
  labs(title = "Top spatially variable genes (nnSVG)",
       subtitle = "Expression overlaid on spatial coordinates",
       x = "spatial x", y = "spatial y") +
  theme_pub(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
save_fig(p1, file.path(ASSETS, "fig1_top_svg_spatial_expression"), width = 8.5, height = 5.6)

## ---- 图2: SVG rank vs LR-stat lollipop (top genes) ------------------------
cat("Step 6: 图2 SVG LR-stat lollipop...\n")
lol <- svg[seq_len(min(20, nrow(svg))), ]
lol$gene <- factor(lol$gene, levels = rev(lol$gene))
lol$class <- cls(as.character(lol$gene))
p2 <- ggplot(lol, aes(x = LR_stat, y = gene, colour = class)) +
  geom_segment(aes(x = 0, xend = LR_stat, yend = gene), linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_pub("npg", name = "ground truth") +
  labs(title = "nnSVG ranking by likelihood-ratio statistic",
       subtitle = "Top 20 genes; colour = synthetic ground-truth class",
       x = "LR statistic (spatial signal)", y = NULL) +
  theme_pub(base_size = 11)
save_fig(p2, file.path(ASSETS, "fig2_svg_lrstat_lollipop"), width = 7.2, height = 6)

## ---- 图3: ★SVG vs HVG 散点 (谁抓空间? 谁只是高变?) -----------------------
cat("Step 7: 图3 ★SVG vs HVG 散点 (诚实基线核心图)...\n")
# x = HVG rank (越小越高变), y = SVG rank (越小越空间); 左下角=两者都靠前
# 关键看点: HVG-only 基因 (橙) 应落在【HVG 靠前 (x 小) 但 SVG 靠后 (y 大)】区域
cmp$neglog_padj <- -log10(pmax(cmp$padj, 1e-300))
p3 <- ggplot(cmp, aes(x = hvg_rank, y = rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = class, size = neglog_padj), alpha = 0.85) +
  scale_color_pub("npg", name = "ground truth") +
  scale_size_continuous(name = expression(-log[10]~padj~(SVG)), range = c(1.2, 5)) +
  scale_y_reverse() + scale_x_reverse() +
  labs(title = "Spatially variable vs highly variable",
       subtitle = "Lower-left = high in both | upper-left = HVG-only (high variance, no spatial structure)",
       x = "HVG rank (non-spatial, by variance)  -->  more variable",
       y = "SVG rank (nnSVG)  -->  more spatial") +
  theme_pub(base_size = 11)
save_fig(p3, file.path(ASSETS, "fig3_svg_vs_hvg_scatter"), width = 7.4, height = 6)

## ---- 图4: 各类别 SVG rank 分布 violin+raincloud (基线对照量化) ------------
cat("Step 8: 图4 类别 x SVG-rank 小提琴...\n")
p4 <- ggplot(cmp, aes(x = class, y = rank, fill = class)) +
  geom_violin(width = 0.9, alpha = 0.5, colour = NA, trim = FALSE) +
  geom_jitter(aes(colour = class), width = 0.12, size = 1.6, alpha = 0.9) +
  scale_fill_pub("npg") + scale_color_pub("npg") +
  scale_y_reverse() +
  guides(fill = "none", colour = "none") +
  labs(title = "SVG rank by ground-truth class",
       subtitle = "nnSVG ranks spatial (gradient/patch) genes high; HVG-only & noise rank low",
       x = NULL, y = "nnSVG rank  -->  more spatial (top)") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(p4, file.path(ASSETS, "fig4_svgrank_by_class_violin"), width = 7.2, height = 5.6)

cat(sprintf("\n完成。结果表见 %s ;展示图见 assets/\n", normalizePath(args$outdir)))
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖快照
