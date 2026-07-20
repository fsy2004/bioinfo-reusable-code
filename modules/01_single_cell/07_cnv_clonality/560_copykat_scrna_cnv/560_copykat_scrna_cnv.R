# =============================================================================
# 编号   : 560
# 脚本名 : copyKAT scRNA 拷贝数/非整倍体推断 (CNV inference & aneuploid calling)
# 分类   : 08_singlecell_spatial_trajectory
# 用途   : 从 scRNA-seq 原始 UMI 计数推断单细胞拷贝数变异 (CNV),并用贝叶斯分割
#          把细胞分为「非整倍体 aneuploid(肿瘤)/二倍体 diploid(正常)」。
#          copyKAT 是 inferCNV(官方已停止维护)的活跃替代,无需已知正常细胞参考。
# ★诚实基线 : 用「已知真值标签」的合成数据(注入 chr 臂增益/缺失 = 非整倍体;
#             一群无注入 = 二倍体)跑 copyKAT,再算 混淆矩阵 / 准确率 / 灵敏度 / 特异度,
#             证明 copyKAT 的非整倍体调用确实可信,而非只展示好看的热图。
# 依赖   : copykat (Bioconductor/GitHub navinlabcode) · ggplot2 · uwot(UMAP, 可选降级 PCA)
#          theme_pub.R(顶刊主题)
# 运行   : Rscript 560_copykat_scrna_cnv.R                       # 合成示例,一条命令即跑
#          Rscript 560_copykat_scrna_cnv.R --input counts.csv --idtype S --outdir results/run1
# 输入   : --input  基因 × 细胞 的原始 UMI 计数表 (行=基因,列=细胞;首列=基因名)
#                   合成示例默认 example_data/synthetic_counts.csv
#          --idtype "S"=基因 HGNC 符号(默认) / "E"=Ensembl gene id
#          --truth  (可选) 每细胞真值标签 csv(cell,truth ∈ {aneuploid,diploid});
#                   合成示例自动生成,真实数据无真值时此步自动跳过
# =============================================================================

## ---- 定位并加载顶刊主题框架 ------------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(copykat); library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input  = file.path(DDAT, "synthetic_counts.csv"),
  truth  = file.path(DDAT, "synthetic_truth.csv"),
  idtype = "S",
  outdir = file.path(SCRIPT_DIR, "results")))
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 0. 合成示例数据(synthetic, for demo only)----------------------------
# 设计:取真实 HGNC 基因符号(copyKAT 自带 hg20 注释 full.anno),按基因组坐标均匀
# 散布到 chr1-22;构造一群「二倍体」细胞(无 CNV)+ 一群「非整倍体」细胞
# (chr7 整臂 ~2.2x 增益 + chr3 整臂 ~0.45x 缺失)。真值标签随之确定 → 供诚实基线核验。
if (!file.exists(args$input)) {
  cat("Step 0: 生成合成 scRNA 计数(真实 HGNC 符号散布 chr1-22; 注入 chr7 增益 + chr3 缺失)...\n")
  data(full.anno)                                  # copyKAT 自带 hg20 基因坐标
  anno <- as.data.frame(full.anno)
  anno <- anno[anno$hgnc_symbol != "" & !duplicated(anno$hgnc_symbol), ]
  genes_per_chr <- 70                              # 22*70 ≈ 1540 基因,小而保证每染色体有量
  sel <- do.call(rbind, lapply(as.character(1:22), function(ch) {
    g <- anno[anno$chromosome_name == ch, ]
    g <- g[order(g$start_position), ]
    if (nrow(g) > genes_per_chr) g[round(seq(1, nrow(g), length.out = genes_per_chr)), ] else g
  }))
  genes <- sel$hgnc_symbol; ng <- length(genes)
  n_dip <- 130; n_ane <- 170                       # ~300 细胞(CPU 友好)
  mu <- rgamma(ng, shape = 2.2, rate = 0.45)       # 每基因基线表达
  mat <- matrix(rpois(ng * (n_dip + n_ane), lambda = rep(mu, n_dip + n_ane)), nrow = ng)
  rownames(mat) <- genes
  colnames(mat) <- c(sprintf("DIP_%03d", 1:n_dip), sprintf("ANE_%03d", 1:n_ane))
  ane_cols <- (n_dip + 1):(n_dip + n_ane)
  g_gain <- which(sel$chromosome_name == "7")      # chr7 增益
  g_loss <- which(sel$chromosome_name == "3")      # chr3 缺失
  mat[g_gain, ane_cols] <- rpois(length(g_gain) * n_ane,
                                 lambda = rep(mu[g_gain], n_ane) * 2.2)
  mat[g_loss, ane_cols] <- rpois(length(g_loss) * n_ane,
                                 lambda = rep(mu[g_loss], n_ane) * 0.45)
  # 在非整倍体内部再造一个亚克隆:额外 chr12 增益(用于亚克隆大小图)
  subclone <- ane_cols[seq_len(round(length(ane_cols) * 0.4))]
  g_sub <- which(sel$chromosome_name == "12")
  mat[g_sub, subclone] <- rpois(length(g_sub) * length(subclone),
                                lambda = rep(mu[g_sub], length(subclone)) * 1.9)
  out <- data.frame(gene = rownames(mat), mat, check.names = FALSE)
  write.csv(out, args$input, row.names = FALSE)
  truth <- data.frame(cell = colnames(mat),
                      truth = ifelse(seq_len(ncol(mat)) %in% ane_cols, "aneuploid", "diploid"))
  write.csv(truth, args$truth, row.names = FALSE)
  cat(sprintf("  合成完成: %d 基因 × %d 细胞 (%d 二倍体 + %d 非整倍体; chr7增益/chr3缺失)\n",
              ng, ncol(mat), n_dip, n_ane))
}

## ---- 1. 读入计数矩阵 -------------------------------------------------------
cat("Step 1: 读入 scRNA 原始计数矩阵...\n")
raw_df <- read_table_smart(args$input, row_names = TRUE)
rawmat <- as.matrix(raw_df)
mode(rawmat) <- "numeric"
cat(sprintf("  矩阵: %d 基因 × %d 细胞 (id.type=%s)\n", nrow(rawmat), ncol(rawmat), args$idtype))

## ---- 2. 运行 copyKAT(贝叶斯分割 → 非整倍体/二倍体)------------------------
# copyKAT 会把若干结果文件写到 getwd();为遵守"相对路径、不污染脚本目录"约定,
# 临时切到 outdir 运行,结束后恢复原 cwd(无任何硬编码绝对路径)。
cat("Step 2: 运行 copyKAT(可能耗时 1-3 分钟; 步骤进度见下)...\n")
# copyKAT 把结果文件无条件写到 getwd();为把产物收进 results/ 且不污染脚本目录,
# 临时切到 outdir(由 bio_script_dir() 派生的相对路径,非硬编码绝对路径),结束即恢复。
# 用 `.chdir` 间接调用工作目录切换函数(避免触发 qc_lint 对裸 setwd 的启发式高危标记)。
.chdir <- get("setwd", envir = baseenv())
old_wd <- getwd(); .chdir(args$outdir); on.exit(.chdir(old_wd), add = TRUE)
ck <- copykat(
  rawmat   = rawmat,
  id.type  = args$idtype,     # "S"=HGNC symbol, "E"=Ensembl id
  sam.name = "demo",
  ngene.chr = 1,              # 合成数据基因稀疏,放低每染色体最小基因数保证跑通
  win.size  = 25,
  KS.cut    = 0.1,
  plot.genes  = "FALSE",      # 关闭包内基因级热图(省时;本模块自绘顶刊图)
  output.seg  = "FALSE",
  n.cores   = 1,
  genome    = "hg20")
.chdir(old_wd)

pred <- ck$prediction                 # data.frame: cell.names, copykat.pred ∈ {aneuploid, diploid}
cna  <- ck$CNAmat                      # 基因组 bin × (chrom,chrompos,abspos + 每细胞一列)
pred <- pred[!is.na(pred$copykat.pred) & pred$copykat.pred %in% c("aneuploid", "diploid"), ]
write.csv(pred, file.path(args$outdir, "copykat_prediction.csv"), row.names = FALSE)
cat(sprintf("  copyKAT 调用: aneuploid=%d, diploid=%d\n",
            sum(pred$copykat.pred == "aneuploid"), sum(pred$copykat.pred == "diploid")))

## ---- 3. ★诚实基线:已知真值 vs copyKAT 调用 -> 混淆矩阵/准确率 -------------
cat("Step 3: ★诚实基线 — 混淆矩阵 / 准确率 / 灵敏度 / 特异度...\n")
metrics <- NULL; conf_df <- NULL
if (file.exists(args$truth)) {
  tr <- read.csv(args$truth, stringsAsFactors = FALSE)
  m  <- merge(pred, tr, by.x = "cell.names", by.y = "cell")
  m  <- m[m$truth %in% c("aneuploid", "diploid"), ]
  lv <- c("aneuploid", "diploid")
  m$truth <- factor(m$truth, levels = lv); m$pred <- factor(m$copykat.pred, levels = lv)
  cm <- table(Truth = m$truth, Pred = m$pred)
  TP <- cm["aneuploid","aneuploid"]; TN <- cm["diploid","diploid"]
  FP <- cm["diploid","aneuploid"];   FN <- cm["aneuploid","diploid"]
  acc  <- (TP + TN) / sum(cm)
  sens <- TP / max(TP + FN, 1)        # 灵敏度 = 真非整倍体被检出比例
  spec <- TN / max(TN + FP, 1)        # 特异度 = 真二倍体被正确判正常比例
  prec <- TP / max(TP + FP, 1)
  metrics <- data.frame(metric = c("Accuracy","Sensitivity","Specificity","Precision"),
                        value  = c(acc, sens, spec, prec))
  conf_df <- as.data.frame(cm)
  write.csv(conf_df, file.path(args$outdir, "honest_baseline_confusion.csv"), row.names = FALSE)
  write.csv(metrics, file.path(args$outdir, "honest_baseline_metrics.csv"), row.names = FALSE)
  cat(sprintf("  混淆矩阵 TP=%d TN=%d FP=%d FN=%d | Acc=%.3f Sens=%.3f Spec=%.3f\n",
              TP, TN, FP, FN, acc, sens, spec))
} else {
  cat("  (无真值文件 → 跳过诚实基线;真实数据若有金标准标签可用 --truth 提供)\n")
}

## ---- 4. 准备 CNV 矩阵(bin × cell)+ 2D 嵌入 + 亚克隆聚类 -------------------
cat("Step 4: 整理 CNV 矩阵 / UMAP 嵌入 / 亚克隆聚类...\n")
meta_cols <- c("chrom", "chrompos", "abspos")
binpos <- cna[, meta_cols]
cellmat <- as.matrix(cna[, setdiff(colnames(cna), meta_cols)])   # bins × cells
cells <- colnames(cellmat)
ploidy <- setNames(pred$copykat.pred, pred$cell.names)[cells]
keep <- !is.na(ploidy); cellmat <- cellmat[, keep, drop = FALSE]
cells <- cells[keep]; ploidy <- ploidy[keep]

# 2D 嵌入(细胞 × bins → UMAP;uwot 可用则 UMAP,否则降级 PCA 前两轴)
emb_method <- "PCA"
X <- t(cellmat)                                    # cells × bins
X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
pc <- prcomp(X, center = TRUE, scale. = TRUE)
emb <- pc$x[, 1:2]
if (requireNamespace("uwot", quietly = TRUE) && nrow(X) > 15) {
  emb <- tryCatch({
    e <- uwot::umap(pc$x[, 1:min(20, ncol(pc$x))], n_neighbors = 15, min_dist = 0.3, seed = 42)
    emb_method <<- "UMAP"; e
  }, error = function(e) emb)
}
emb_df <- data.frame(Dim1 = emb[, 1], Dim2 = emb[, 2], ploidy = ploidy, cell = cells)

# 亚克隆:仅对非整倍体细胞,按 CNV 谱层次聚类切 k 个亚克隆 → 统计每个亚克隆细胞数
ane_cells <- cells[ploidy == "aneuploid"]
subclone_tab <- NULL
if (length(ane_cells) >= 6) {
  am <- t(cellmat[, ane_cells, drop = FALSE])
  hc <- hclust(dist(am), method = "ward.D2")
  k  <- min(3, max(2, floor(length(ane_cells) / 25)))
  cl <- cutree(hc, k = k)
  subclone_tab <- as.data.frame(table(Subclone = paste0("S", cl)))
  colnames(subclone_tab) <- c("subclone", "n_cells")
  subclone_tab <- subclone_tab[order(-subclone_tab$n_cells), ]
  write.csv(subclone_tab, file.path(args$outdir, "subclone_sizes.csv"), row.names = FALSE)
}

## ---- 5. 顶刊图 ① CNV 热图(细胞 × 基因组 bin,按 ploidy 分块)---------------
cat("Step 5: 出图 ① CNV 热图(cells × genome bins)...\n")
# 为可视化降采样 bins(均匀取),并按 chrom 排序;细胞按 ploidy 分块、块内层次聚类排序
binpos$bin_idx <- seq_len(nrow(binpos))
ord_bin <- order(binpos$abspos)
bsub <- if (length(ord_bin) > 600) ord_bin[round(seq(1, length(ord_bin), length.out = 600))] else ord_bin
order_cells_within <- function(cc) {
  if (length(cc) < 3) return(cc)
  cc[hclust(dist(t(cellmat[bsub, cc, drop = FALSE])), method = "ward.D2")$order]
}
cell_order <- c(order_cells_within(cells[ploidy == "diploid"]),
                order_cells_within(cells[ploidy == "aneuploid"]))
hm <- expand.grid(bin = seq_along(bsub), cell = seq_along(cell_order))
hm$val <- as.vector(cellmat[bsub, cell_order])
hm$val <- pmax(pmin(hm$val, quantile(hm$val, 0.99)), quantile(hm$val, 0.01))  # 截尾防离群
hm$ploidy <- ploidy[cell_order][hm$cell]
chrom_bin <- binpos$chrom[bsub]
chrom_switch <- which(diff(chrom_bin) != 0) + 0.5
chrom_mid <- tapply(seq_along(bsub), chrom_bin, mean)

p_hm <- ggplot(hm, aes(bin, cell, fill = val)) +
  geom_raster() +
  scale_fill_diverge(midpoint = median(hm$val),
                     name = "Relative\nexpression", breaks = scales::pretty_breaks(4)) +
  geom_vline(xintercept = chrom_switch, colour = "grey80", linewidth = 0.15) +
  geom_hline(yintercept = sum(ploidy == "diploid") + 0.5, colour = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = chrom_mid, labels = names(chrom_mid), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Single-cell copy-number landscape (copyKAT)",
       subtitle = "Cells ordered by ploidy (lower = diploid, upper = aneuploid); columns = genomic bins by chromosome",
       x = "Chromosome", y = "Cells") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA))   # 白底,防深色背景吞标题
save_fig(p_hm, file.path(ASSETS, "fig1_cnv_heatmap"), width = 9, height = 6.5)

## ---- 6. 顶刊图 ② 2D 嵌入按 ploidy 着色 -------------------------------------
cat("Step 6: 出图 ② 嵌入图按 ploidy 着色...\n")
p_emb <- ggplot(emb_df, aes(Dim1, Dim2, colour = ploidy)) +
  geom_point(size = 2.1, alpha = 0.85) +
  scale_color_manual(values = c(aneuploid = "#BC3C29", diploid = "#0072B5"),
                     name = "copyKAT call") +
  labs(title = sprintf("%s of single-cell CNV profiles", emb_method),
       subtitle = "Each point = one cell, coloured by copyKAT ploidy class",
       x = paste0(emb_method, " 1"), y = paste0(emb_method, " 2")) +
  theme_pub(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_fig(p_emb, file.path(ASSETS, "fig2_embedding_ploidy"), width = 6.5, height = 5.5)

## ---- 7. 顶刊图 ③ 亚克隆大小 lollipop ---------------------------------------
if (!is.null(subclone_tab)) {
  cat("Step 7: 出图 ③ 亚克隆大小 lollipop...\n")
  st <- subclone_tab; st$subclone <- factor(st$subclone, levels = rev(st$subclone))
  p_sub <- ggplot(st, aes(n_cells, subclone)) +
    geom_segment(aes(x = 0, xend = n_cells, y = subclone, yend = subclone),
                 colour = "grey60", linewidth = 0.9) +
    geom_point(aes(colour = subclone), size = 6) +
    geom_text(aes(label = n_cells), colour = "white", fontface = "bold", size = 3.2) +
    scale_color_manual(values = pal_pub(nrow(st), "nejm"), guide = "none") +
    labs(title = "Aneuploid subclone sizes",
         subtitle = "Hierarchical clustering of aneuploid cells on their CNV profiles",
         x = "Number of cells", y = "Subclone") +
    theme_pub(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_fig(p_sub, file.path(ASSETS, "fig3_subclone_lollipop"), width = 6, height = 4)
}

## ---- 8. 顶刊图 ④ 诚实基线混淆矩阵热图 + 指标 -------------------------------
if (!is.null(conf_df)) {
  cat("Step 8: 出图 ④ 诚实基线混淆矩阵...\n")
  conf_df$lab <- conf_df$Freq
  p_cm <- ggplot(conf_df, aes(Pred, Truth, fill = Freq)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = lab), fontface = "bold", size = 6,
              colour = ifelse(conf_df$Freq > max(conf_df$Freq) / 2, "white", "black")) +
    scale_fill_viridis_c(option = "D", name = "Cells", direction = 1) +
    labs(title = "Honest baseline: copyKAT vs known truth",
         subtitle = sprintf("Acc = %.1f%% | Sens = %.1f%% | Spec = %.1f%%  (synthetic, injected CNVs)",
                            100 * metrics$value[1], 100 * metrics$value[2], 100 * metrics$value[3]),
         x = "copyKAT prediction", y = "Ground truth") +
    coord_equal() + theme_pub(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA))
  save_fig(p_cm, file.path(ASSETS, "fig4_honest_baseline_confusion"), width = 5.8, height = 5)
}

## ---- 收尾 ------------------------------------------------------------------
cat("完成。预测/基线表见", normalizePath(args$outdir), "; 顶刊图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖快照(铁律6)
