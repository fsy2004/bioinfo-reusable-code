# =============================================================================
# 编号   : 545
# 脚本名 : SPOTlight 空间 spot 细胞类型去卷积(NMF 回归)+ scatterpie 空间饼图
# 分类   : 16_spatial_communication
# 用途   : 用带标签的 scRNA 参考,对空间转录组每个 spot 做细胞类型去卷积
#          (SPOTlight = NMF 主题建模 + NNLS 回归),得到每 spot 的细胞类型比例,
#          并以 scatterpie 空间饼图 / 比例热图 直观展示组织微环境构成。
# ★诚实基线 : 合成 spot 时混合比例【已知】→ 把 SPOTlight 预测比例 vs 真实比例
#            算 RMSE / Pearson r(全局 + 逐 spot),验证去卷积是否可信;
#            并对比一个"纯丰度"朴素基线(spot 表达直接打 marker 分),
#            证明 NMF 去卷积确实优于朴素打分,避免只报好看指标。
# 依赖   : SPOTlight(NMF 去卷积) · SingleCellExperiment · scran · scuttle ·
#          scatterpie(空间饼) · ggplot2 · 框架 theme_pub.R
# 运行   : Rscript 545_spotlight_deconvolution.R                  # 合成示例,CPU 秒级
#          Rscript 545_spotlight_deconvolution.R --outdir results/run1
# 输入   : 合成 demo,脚本内生成(synthetic, for demo only):
#          ① scRNA 参考 counts(基因 × 细胞)+ 每细胞 cell_type 标签;
#          ② spatial spot counts(基因 × spot)+ 每 spot 的【真实】细胞类型比例。
#          换真实数据时:参考 SingleCellExperiment(counts + colData$cell_type),
#          spatial 为 counts 矩阵(基因行 × spot 列)+ spot 的 x/y 坐标。
# =============================================================================

## ---- 0. 载入框架(顶刊主题)+ 依赖 ----------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(SPOTlight); library(SingleCellExperiment)
  library(scran); library(scuttle); library(ggplot2)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  outdir   = file.path(SCRIPT_DIR, "results"),
  n_types  = 4,        # 细胞类型数
  n_per    = 60,       # 每类参考细胞数
  n_genes  = 200,      # 基因总数
  n_spots  = 100,      # 空间 spot 数
  cells_per_spot = 30  # 每 spot 由多少个细胞混合而成
))
for (k in c("n_types","n_per","n_genes","n_spots","cells_per_spot"))
  args[[k]] <- as.integer(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. 合成 scRNA 参考:每类细胞有专属 marker(synthetic, demo only) -----
# 思路:n_types 类细胞,基因平均分给各类作 marker;某类细胞在其 marker 上高表达,
#       其余基因低背景表达;用负二项(rnbinom)出整数 counts,贴近真实 scRNA。
cat("Step 1: 合成带 marker 的 scRNA 参考...\n")
nt <- args$n_types; ng <- args$n_genes; npc <- args$n_per
type_names <- paste0("CellType", LETTERS[seq_len(nt)])
genes <- sprintf("gene%03d", seq_len(ng))
# 每个基因指派一个"高表达归属"类型(marker),近似均分
gene_owner <- type_names[((seq_len(ng) - 1) %% nt) + 1]

bg_mu  <- 0.4    # 背景表达均值
mk_mu  <- 8.0    # marker 高表达均值
size   <- 1.5    # 负二项离散度(越小越离散)

sim_cell <- function(ctype) {
  mu <- ifelse(gene_owner == ctype, mk_mu, bg_mu)
  rnbinom(ng, mu = mu, size = size)
}
sc_counts <- do.call(cbind, lapply(rep(type_names, each = npc), sim_cell))
rownames(sc_counts) <- genes
colnames(sc_counts) <- sprintf("cell%04d", seq_len(ncol(sc_counts)))
cell_type <- factor(rep(type_names, each = npc), levels = type_names)

sce <- SingleCellExperiment(assays = list(counts = sc_counts),
                            colData = DataFrame(cell_type = cell_type))
cat(sprintf("  参考: %d 基因 × %d 细胞 (%d 类 × %d/类)\n", ng, ncol(sce), nt, npc))

## 落盘示例数据(便于复用/查看;真实数据替换此处)
write.csv(as.data.frame(as.matrix(sc_counts)),
          file.path(DDAT, "scRNA_reference_counts.csv"))
write.csv(data.frame(cell = colnames(sce), cell_type = as.character(cell_type)),
          file.path(DDAT, "scRNA_reference_labels.csv"), row.names = FALSE)

## ---- 2. 合成空间 spot:每 spot 用【已知】比例混合细胞(ground truth) -------
cat("Step 2: 合成空间 spot(已知混合比例 = 诚实基线 ground truth)...\n")
ns  <- args$n_spots; cps <- args$cells_per_spot
# 设 4 个空间"区域",各区域以某一主导细胞类型为主 + 少量其它(造空间结构,利于饼图)
# 自包含 Dirichlet 采样(gamma 归一化,无需额外依赖)
rdirichlet1 <- function(alpha) { g <- rgamma(length(alpha), shape = alpha, rate = 1); g / sum(g) }
region <- sample(seq_len(nt), ns, replace = TRUE)
true_prop <- matrix(0, nrow = ns, ncol = nt, dimnames = list(
  sprintf("spot%03d", seq_len(ns)), type_names))
for (s in seq_len(ns)) {
  alpha <- rep(1, nt); alpha[region[s]] <- 6   # Dirichlet:主导类型权重高
  true_prop[s, ] <- rdirichlet1(alpha)
}
# 按真实比例从参考细胞中抽样并把 counts 相加 → spot 表达(贴近真实 spot=细胞混合)
sp_counts <- matrix(0, nrow = ng, ncol = ns, dimnames = list(genes, rownames(true_prop)))
ref_by_type <- split(seq_len(ncol(sce)), cell_type)
for (s in seq_len(ns)) {
  n_each <- round(true_prop[s, ] * cps); n_each[n_each < 0] <- 0
  if (sum(n_each) == 0) n_each[region[s]] <- cps
  picked <- unlist(lapply(seq_len(nt), function(k)
    if (n_each[k] > 0) sample(ref_by_type[[k]], n_each[k], replace = TRUE)))
  sp_counts[, s] <- rowSums(sc_counts[, picked, drop = FALSE])
  true_prop[s, ] <- n_each / sum(n_each)   # 用实际抽到的细胞数重算真实比例
}
# 空间坐标:按区域聚成 4 个空间簇(让 scatterpie 呈现空间格局)
centers <- data.frame(cx = c(3, 9, 3, 9)[seq_len(nt)], cy = c(3, 3, 9, 9)[seq_len(nt)])
coords <- data.frame(
  x = centers$cx[region] + rnorm(ns, 0, 1.1),
  y = centers$cy[region] + rnorm(ns, 0, 1.1))
rownames(coords) <- rownames(true_prop)

spe <- SingleCellExperiment(assays = list(counts = sp_counts),
                            colData = DataFrame(x = coords$x, y = coords$y))
write.csv(cbind(coords, as.data.frame(round(true_prop, 4))),
          file.path(DDAT, "spatial_spots_truth.csv"))
cat(sprintf("  spot: %d 个 (每 spot ~%d 细胞混合, 4 空间区域)\n", ns, cps))

## ---- 3. 标记基因(scran::scoreMarkers)→ SPOTlight 输入 mgs ---------------
cat("Step 3: logNorm + scoreMarkers 取 marker gene set(mgs)...\n")
sce <- scuttle::logNormCounts(sce)
mk  <- scran::scoreMarkers(sce, groups = colData(sce)$cell_type)
mgs <- do.call(rbind, lapply(names(mk), function(g) {
  x <- mk[[g]]; x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x <- x[x$mean.AUC > 0.6, ]
  if (nrow(x) == 0) x <- head(mk[[g]][order(mk[[g]]$mean.AUC, decreasing = TRUE), ], 10)
  data.frame(gene = rownames(x), cluster = g, weight = x$mean.AUC)
}))
cat(sprintf("  marker gene set: %d 条 (%d 类)\n", nrow(mgs), length(mk)))

## ---- 4. SPOTlight 去卷积(NMF 主题建模 + NNLS 回归) ----------------------
cat("Step 4: SPOTlight NMF 去卷积...\n")
spot_res <- SPOTlight(
  x = sce, y = spe,
  groups   = as.character(colData(sce)$cell_type),
  mgs      = mgs,
  gene_id  = "gene", group_id = "cluster", weight_id = "weight")
pred <- spot_res$mat                      # spots × cell_types 预测比例(行和=1)
pred <- pred[rownames(true_prop), type_names, drop = FALSE]
write.csv(as.data.frame(round(pred, 4)),
          file.path(args$outdir, "predicted_proportions.csv"))

## ---- 5. ★诚实基线:预测 vs 真实(全局 + 逐 spot)+ 朴素基线对照 ----------
cat("Step 5: ★诚实基线评估(预测 vs 真实 RMSE / 相关)...\n")
rmse <- function(a, b) sqrt(mean((a - b)^2))
spot_rmse <- vapply(seq_len(ns), function(i) rmse(pred[i, ], true_prop[i, ]), numeric(1))
spot_cor  <- vapply(seq_len(ns), function(i) {
  v <- suppressWarnings(cor(pred[i, ], true_prop[i, ])); if (is.na(v)) 0 else v }, numeric(1))
global_rmse <- rmse(as.vector(pred), as.vector(true_prop))
global_cor  <- cor(as.vector(pred), as.vector(true_prop))

# 朴素基线:不做 NMF,直接用 spot 表达对各类 marker 求平均表达 → 归一化为"比例"
sp_log <- log1p(t(t(sp_counts) / (colSums(sp_counts) / 1e4)))   # CPM-ish logNorm
naive <- sapply(type_names, function(g) {
  mk_g <- mgs$gene[mgs$cluster == g]; colMeans(sp_log[mk_g, , drop = FALSE]) })
naive <- naive / rowSums(naive)
naive <- naive[rownames(true_prop), type_names, drop = FALSE]
naive_rmse <- rmse(as.vector(naive), as.vector(true_prop))
naive_cor  <- cor(as.vector(naive), as.vector(true_prop))

eval_tab <- data.frame(
  method = c("SPOTlight (NMF)", "Naive marker-score baseline"),
  global_RMSE = round(c(global_rmse, naive_rmse), 4),
  global_Pearson_r = round(c(global_cor, naive_cor), 4))
write.csv(eval_tab, file.path(args$outdir, "honest_baseline_eval.csv"), row.names = FALSE)
cat(sprintf("  SPOTlight : RMSE=%.3f  r=%.3f  (逐spot RMSE 中位=%.3f)\n",
            global_rmse, global_cor, median(spot_rmse)))
cat(sprintf("  Naive base: RMSE=%.3f  r=%.3f  → NMF %s 朴素基线\n",
            naive_rmse, naive_cor, ifelse(global_rmse < naive_rmse, "优于", "未优于(查参数)")))

## =============================================================================
## 出图(顶刊风格;禁止平凡条形图 → 空间饼/热图/散点)
## =============================================================================
cat("Step 6: 出图(scatterpie 空间饼 / 比例热图 / 预测vs真实散点)...\n")
pal <- pal_pub(nt, "npg"); names(pal) <- type_names

## 图1:scatterpie 空间饼图 —— 每 spot 一个饼,扇区=各细胞类型预测比例 -------
df_pie <- data.frame(coords, as.data.frame(pred), check.names = FALSE)
df_pie$x <- coords$x; df_pie$y <- coords$y
p_pie <- ggplot() +
  scatterpie::geom_scatterpie(
    data = df_pie, aes(x = x, y = y, r = 0.42),
    cols = type_names, color = NA, alpha = 0.95) +
  scale_fill_manual(values = pal, name = "Cell type") +
  coord_equal() +
  labs(title = "SPOTlight spatial deconvolution",
       subtitle = "Each pie = one spot; sectors = predicted cell-type proportions",
       x = "Spatial X", y = "Spatial Y") +
  theme_pub(base_size = 12) +
  theme(legend.position = "right")
save_fig(p_pie, file.path(ASSETS, "spatial_scatterpie"),
         width = 7.2, height = 6)

## 图2:比例热图 —— spot × cell type 预测比例(viridis 连续色)-------------
hm_long <- data.frame(
  spot = factor(rep(rownames(pred), nt), levels = rev(rownames(pred))),
  cell_type = factor(rep(type_names, each = ns), levels = type_names),
  prop = as.vector(pred))
# spot 按主导类型排序,便于看出块状结构
ord <- order(region, -apply(pred, 1, max))
hm_long$spot <- factor(rep(rownames(pred), nt),
                       levels = rownames(pred)[rev(ord)])
p_hm <- ggplot(hm_long, aes(cell_type, spot, fill = prop)) +
  geom_tile() +
  scale_fill_cont(option = "D", name = "Proportion", limits = c(0, 1)) +
  labs(title = "Predicted cell-type proportion per spot",
       subtitle = "Spots ordered by dominant cell type",
       x = "Cell type", y = "Spot") +
  theme_pub(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_hm, file.path(ASSETS, "proportion_heatmap"),
         width = 5.6, height = 7)

## 图3:★诚实基线散点 —— 预测比例 vs 真实比例(对角线=完美)----------------
df_sc <- data.frame(
  true = as.vector(true_prop), pred = as.vector(pred),
  cell_type = factor(rep(type_names, each = ns), levels = type_names))
lab <- sprintf("Global RMSE = %.3f\nPearson r = %.3f", global_rmse, global_cor)
p_sc <- ggplot(df_sc, aes(true, pred, color = cell_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_color_manual(values = pal, name = "Cell type") +
  annotate("text", x = 0.02, y = 0.97, hjust = 0, vjust = 1,
           label = lab, size = 3.6, fontface = "bold") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Honest baseline: predicted vs true proportions",
       subtitle = "Known synthetic mixtures; points on dashed line = perfect",
       x = "True proportion", y = "Predicted proportion") +
  theme_pub(base_size = 12)
save_fig(p_sc, file.path(ASSETS, "pred_vs_true_scatter"),
         width = 6.2, height = 6)

## 图4:方法对比 —— SPOTlight vs 朴素基线 逐 spot RMSE(raincloud/violin+jitter)
rc <- data.frame(
  method = factor(rep(c("SPOTlight\n(NMF)", "Naive\nbaseline"), each = ns),
                  levels = c("Naive\nbaseline", "SPOTlight\n(NMF)")),
  rmse = c(spot_rmse,
           vapply(seq_len(ns), function(i) rmse(naive[i, ], true_prop[i, ]), numeric(1))))
p_rc <- ggplot(rc, aes(method, rmse, fill = method, color = method)) +
  geom_violin(width = 0.85, alpha = 0.25, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.16, alpha = 0.9, outlier.shape = NA, color = "grey20") +
  geom_jitter(width = 0.07, size = 1.1, alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35"), guide = "none") +
  scale_color_manual(values = c("#4DBBD5", "#E64B35"), guide = "none") +
  labs(title = "Per-spot deconvolution error",
       subtitle = "Lower = better; NMF deconvolution vs naive marker scoring",
       x = NULL, y = "Per-spot RMSE (pred vs true)") +
  theme_pub(base_size = 12)
save_fig(p_rc, file.path(ASSETS, "method_rmse_violin"),
         width = 5.2, height = 5.4)

cat("完成。结果表见", normalizePath(args$outdir), ";展示图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖快照(铁律6)
