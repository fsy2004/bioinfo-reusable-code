# =============================================================================
# 编号   : 559
# 脚本名 : 单细胞多样本多条件 pseudobulk 差异状态(DS) — muscat 金标准
# 分类   : 03_transcriptomics_deg
# 用途   : 多样本×多条件 scRNA 的差异状态分析:把同一(细胞类型 × 样本)的细胞
#          聚合成 pseudobulk profile,再在样本水平上做 edgeR / limma-voom / DESeq2
#          差异检验。这是 2026 公认金标准,避免"把单个细胞当独立样本"造成的
#          伪重复(pseudo-replication)与 I 类错误膨胀。
#          (muscat: Crowell et al. Nat Commun 2020, "muscat detects ... DS")
#
# ★诚实基线(第7类铁律 · 防伪重复):脚本内置对照——同一份数据,
#          ① pseudobulk(样本级聚合 + edgeR)=正确做法;
#          ② cell-level DE(把每个细胞当独立样本 + edgeR)=错误做法。
#          合成数据在"零真实条件效应"的基因上注入了样本(供体)随机效应。
#          预期结果:cell-level 把大量这类 null 基因误判为显著(假阳性爆炸),
#          pseudobulk 把 I 类错误拉回名义水平 → 直接量化展示伪重复的危害。
#
# 依赖   : muscat · SingleCellExperiment · scater · edgeR · limma · (DESeq2 可选)
#          · ggplot2 (+ ggrepel 可选);全部已安装
# 运行   : Rscript 559_muscat_pseudobulk_ds.R                  # 零改动跑合成示例
#          Rscript 559_muscat_pseudobulk_ds.R --input my_sce.rds --outdir results/run1
# 输入   : --input 一个 .rds(SingleCellExperiment),colData 须含三列:
#          cluster_id(细胞类型) / sample_id(样本) / group_id(条件,2 水平)。
#          assays 须含 counts。缺省 --input 时脚本内生成合成 SCE(synthetic, demo only)。
# =============================================================================

## ---- 0. 载入框架(theme_pub.R)+ 依赖 + 参数 -------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(muscat); library(SingleCellExperiment); library(scater)
  library(edgeR);  library(ggplot2)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input  = file.path(DDAT, "demo_sce.rds"),   # 缺省合成;真实数据传 SCE .rds
  outdir = file.path(SCRIPT_DIR, "results"),
  method = "edgeR",                            # pseudobulk 引擎: edgeR / limma-voom / DESeq2
  min_cells = 10,                              # pbDS 每(簇,样本)最少细胞数
  fdr = 0.05))
args$min_cells <- as.numeric(args$min_cells); args$fdr <- as.numeric(args$fdr)
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
PAL <- pal_pub(name = "npg")

## ---- 1. 合成 SCE(synthetic, demo only)或读入真实 SCE ----------------------
# 设计(为演示伪重复而刻意构造):
#   · 6 样本 = 2 条件(ctrl/stim)× 3 生物学重复;2 细胞类型 cluster1/cluster2。
#   · 真实 DS 基因(de_genes):仅在 cluster1 的 stim 中均值×fold,是真信号。
#   · ★null_donor 基因:条件间无系统差异,但带有"样本(供体)随机效应"——
#     即每个样本整体上调/下调。cell-level DE 会把供体波动误当条件效应 → 假阳性;
#     pseudobulk 在 6 个样本水平检验,正确判其不显著。
make_demo_sce <- function(path) {
  ng <- 400; nk <- 2; reps <- 3
  grp_lv <- c("ctrl", "stim")
  samples <- do.call(rbind, lapply(grp_lv, function(g)
    data.frame(grp = g, sid = paste0(g, 1:reps), stringsAsFactors = FALSE)))
  cells_per <- 70                         # 每(样本×簇)细胞数
  cell_sid <- cell_grp <- cell_kid <- character(0)
  for (i in seq_len(nrow(samples))) for (k in seq_len(nk)) {
    cell_sid <- c(cell_sid, rep(samples$sid[i], cells_per))
    cell_grp <- c(cell_grp, rep(samples$grp[i], cells_per))
    cell_kid <- c(cell_kid, rep(paste0("cluster", k), cells_per))
  }
  nc <- length(cell_sid)
  base_mu <- rgamma(ng, shape = 2, rate = 0.4) + 0.5
  de_genes   <- 1:40                       # 真 DS(cluster1, stim↑)
  null_donor <- 41:120                     # 无条件效应但有供体随机效应(伪重复陷阱)
  # 每个样本一份供体随机乘子(仅作用于 null_donor 基因)
  donor_fac <- matrix(exp(rnorm(length(null_donor) * nrow(samples), 0, 0.6)),
                      nrow = length(null_donor),
                      dimnames = list(NULL, samples$sid))
  counts <- matrix(0L, ng, nc)
  for (j in seq_len(nc)) {
    mu <- base_mu
    if (cell_kid[j] == "cluster1" && cell_grp[j] == "stim")
      mu[de_genes] <- mu[de_genes] * 3.0                       # 真信号
    mu[null_donor] <- mu[null_donor] * donor_fac[, cell_sid[j]]# 供体随机效应(null)
    counts[, j] <- rnbinom(ng, mu = mu, size = 2)
  }
  rownames(counts) <- paste0("gene", seq_len(ng))
  colnames(counts) <- paste0("cell", seq_len(nc))
  cd <- DataFrame(cluster_id = cell_kid, sample_id = cell_sid, group_id = cell_grp)
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = cd)
  metadata(sce)$truth <- list(de_genes = rownames(counts)[de_genes],
                              null_donor = rownames(counts)[null_donor])
  saveRDS(sce, path)
  cat(sprintf("[gen] 合成 SCE: %d genes × %d cells | 6 样本(2 条件×3 重复)× 2 细胞类型 | 真 DS=%d, null+供体效应=%d\n",
              ng, nc, length(de_genes), length(null_donor)))
  sce
}

cat("Step 1: 载入 SingleCellExperiment...\n")
if (file.exists(args$input)) {
  sce <- readRDS(args$input)
  cat(sprintf("  读入 %s | %d genes × %d cells\n", basename(args$input), nrow(sce), ncol(sce)))
} else {
  sce <- make_demo_sce(args$input)
}
stopifnot(all(c("cluster_id","sample_id","group_id") %in% colnames(colData(sce))))
if (!"logcounts" %in% assayNames(sce)) sce <- logNormCounts(sce)

## ---- 2. muscat 标准 pseudobulk DS 流程 ------------------------------------
cat("Step 2: prepSCE → aggregateData(by cluster×sample) → pbDS...\n")
sce <- prepSCE(sce, kid = "cluster_id", sid = "sample_id", gid = "group_id", drop = TRUE)
ei  <- metadata(sce)$experiment_info
# 样本级聚合:同一(cluster_id, sample_id)的细胞 counts 求和 = pseudobulk profile
pb  <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
cat(sprintf("  pseudobulk: %d 个细胞类型 × %d 个样本/类型 (assay=cluster)\n",
            length(assayNames(pb)), ncol(pb)))
# DS 检验:每个细胞类型内,样本水平 edgeR/limma-voom/DESeq2,对 group_id 做对比
ds <- pbDS(pb, method = args$method, min_cells = args$min_cells, verbose = FALSE)
res_pb <- resDS(sce, ds, bind = "row")            # tidy: gene/cluster_id/logFC/p_val/p_adj.loc/...
res_pb$padj <- res_pb$p_adj.loc                   # 簇内 BH 校正(局部)
write.csv(res_pb, file.path(args$outdir, "pseudobulk_DS_results.csv"), row.names = FALSE)
contrast_nm <- names(ds$table)[1]
cat(sprintf("  pbDS(%s) 完成 | 对比=%s | 显著 DS 基因(padj<%.2g): %d\n",
            args$method, contrast_nm, args$fdr, sum(res_pb$padj < args$fdr, na.rm = TRUE)))

## ---- 3. ★诚实基线:把每个细胞当独立样本的 cell-level DE(错误做法)---------
# 对每个细胞类型,直接在细胞水平用 edgeR QLF 检验 group_id(忽略样本结构)。
# 这正是常被诟病的伪重复:有效样本量被虚增到细胞数,p 值被严重低估。
cat("Step 3: ★诚实基线 — cell-level DE(伪重复对照)...\n")
cell_level_de <- function(sce, kid) {
  sub <- sce[, sce$cluster_id == kid]
  cnt <- as.matrix(counts(sub))
  grp <- factor(sub$group_id)
  keep <- rowSums(cnt >= 1) >= 10
  cnt <- cnt[keep, , drop = FALSE]
  y <- DGEList(cnt, group = grp); y <- calcNormFactors(y)
  des <- model.matrix(~ grp)
  y <- estimateDisp(y, des); fit <- glmQLFit(y, des)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tt), cluster_id = kid,
             logFC = tt$logFC, p_val = tt$PValue, padj = tt$FDR,
             stringsAsFactors = FALSE)
}
res_cell <- do.call(rbind, lapply(levels(sce$cluster_id), function(k) cell_level_de(sce, k)))
write.csv(res_cell, file.path(args$outdir, "celllevel_DE_results.csv"), row.names = FALSE)

## ---- 3b. 量化对照:在 truth 已知的合成数据上算 I 类错误 / 召回 -------------
truth <- metadata(sce)$truth
cmp_tbl <- NULL
if (!is.null(truth)) {
  # 只看含伪重复陷阱的 cluster1(真 DS + null+供体效应 都在 cluster1)
  pick <- function(df, clu) df[df$cluster_id == clu, ]
  c1_pb   <- pick(res_pb,   "cluster1"); c1_cell <- pick(res_cell, "cluster1")
  fpr <- function(df) {            # null+供体基因里被判显著的比例 = 实测 I 类错误
    n <- df[df$gene %in% truth$null_donor, ]
    mean(n$padj < args$fdr, na.rm = TRUE) }
  tpr <- function(df) {            # 真 DS 基因被检出的比例 = 召回
    d <- df[df$gene %in% truth$de_genes, ]
    mean(d$padj < args$fdr, na.rm = TRUE) }
  cmp_tbl <- data.frame(
    method = c("Pseudobulk (muscat, correct)", "Cell-level DE (pseudo-replication)"),
    FalsePositiveRate_nullDonor = c(fpr(c1_pb), fpr(c1_cell)),
    Recall_trueDS               = c(tpr(c1_pb), tpr(c1_cell)),
    n_sig_total = c(sum(c1_pb$padj   < args$fdr, na.rm = TRUE),
                    sum(c1_cell$padj < args$fdr, na.rm = TRUE)))
  write.csv(cmp_tbl, file.path(args$outdir, "honest_baseline_comparison.csv"), row.names = FALSE)
  cat("  —— 诚实基线实测(cluster1, FDR<", args$fdr, ") ——\n", sep = "")
  print(cmp_tbl, row.names = FALSE)
  cat(sprintf("  解读: 供体随机效应基因本无真实条件差异;cell-level 假阳性率 %.0f%% (伪重复),pseudobulk %.0f%%。\n",
              100 * cmp_tbl$FalsePositiveRate_nullDonor[2], 100 * cmp_tbl$FalsePositiveRate_nullDonor[1]))
}

## ---- 4. 顶刊级图(全部非平凡条形;每图独立成文件)---------------------------
cat("Step 4: 出图(pseudobulk MDS / 火山 / DS 热图 / top-DS lollipop / 基线对照)...\n")

# 4.1 pseudobulk MDS(按样本)—— 样本是否按条件分离(质控核心图)
mds_df <- NULL
try({
  y <- DGEList(assay(pb, "cluster1"))
  y <- calcNormFactors(y)
  lcpm <- edgeR::cpm(y, log = TRUE)
  mds  <- limma::plotMDS(lcpm, plot = FALSE)
  sid  <- colnames(assay(pb, "cluster1"))
  grp  <- ei$group_id[match(sid, ei$sample_id)]
  mds_df <- data.frame(MDS1 = mds$x, MDS2 = mds$y, sample = sid, group = grp)
}, silent = TRUE)
if (!is.null(mds_df)) {
  p_mds <- ggplot(mds_df, aes(MDS1, MDS2, color = group, label = sample)) +
    geom_point(size = 4, alpha = 0.9) +
    { if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_text_repel(size = 3, show.legend = FALSE)
      else geom_text(vjust = -1, size = 3, show.legend = FALSE) } +
    scale_color_manual(values = PAL) +
    labs(title = "Pseudobulk MDS (cluster1)", subtitle = "Each point = one sample's aggregated profile",
         x = "MDS dim 1", y = "MDS dim 2", color = "Condition") +
    theme_pub(base_size = 12)
  save_fig(p_mds, file.path(ASSETS, "fig1_pseudobulk_mds"), width = 6, height = 5)
}

# 4.2 pathway-guided 火山(每细胞类型 facet;来自 pseudobulk DS 结果)
volc <- res_pb
volc$sig <- ifelse(volc$padj < args$fdr & abs(volc$logFC) > 1,
                   ifelse(volc$logFC > 0, "Up", "Down"), "n.s.")
volc$neglog10p <- -log10(pmax(volc$p_val, 1e-300))
top_lab <- do.call(rbind, lapply(split(volc, volc$cluster_id), function(d) {
  d <- d[d$sig != "n.s.", ]; if (!nrow(d)) return(NULL)
  d[order(d$p_val), ][seq_len(min(8, nrow(d))), ] }))
p_volc <- ggplot(volc, aes(logFC, neglog10p, color = sig)) +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(args$fdr), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  { if (!is.null(top_lab) && requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 2.6,
                               max.overlaps = 20, show.legend = FALSE) } +
  facet_wrap(~ cluster_id, scales = "free_y") +
  scale_color_manual(values = c(Up = "#E64B35", Down = "#3C5488", "n.s." = "grey75")) +
  labs(title = "Pseudobulk DS volcano per cell type",
       subtitle = sprintf("muscat %s · contrast = %s", args$method, contrast_nm),
       x = expression(log[2]~fold~change), y = expression(-log[10]~italic(p)), color = NULL) +
  theme_pub(base_size = 12)
save_fig(p_volc, file.path(ASSETS, "fig2_volcano_per_celltype"), width = 8, height = 4.6)

# 4.3 行 z-score DS 热图(top DS 基因 × 样本,按条件注释)
sig_genes <- res_pb[res_pb$cluster_id == "cluster1" & res_pb$padj < args$fdr, ]
sig_genes <- head(sig_genes[order(sig_genes$p_val), "gene"], 25)
if (length(sig_genes) >= 2) {
  lcpm1 <- edgeR::cpm(DGEList(assay(pb, "cluster1")), log = TRUE)[sig_genes, , drop = FALSE]
  z <- t(scale(t(lcpm1)))                                  # 行 z-score
  sid <- colnames(z); grp <- ei$group_id[match(sid, ei$sample_id)]
  ord <- order(grp, sid)
  hm <- data.frame(gene = rep(rownames(z), ncol(z)),
                   sample = rep(sid, each = nrow(z)),
                   z = as.vector(z))
  hm$sample <- factor(hm$sample, levels = sid[ord])
  hm$gene   <- factor(hm$gene, levels = rownames(z)[order(rowMeans(z[, grp == "stim", drop = FALSE]))])
  ann <- data.frame(sample = factor(sid[ord], levels = sid[ord]),
                    group = grp[ord], y = 0)
  p_hm <- ggplot(hm, aes(sample, gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_diverge(midpoint = 0, name = "Row z") +
    labs(title = "Top DS genes (cluster1) — pseudobulk z-score",
         subtitle = "Columns = samples, grouped by condition", x = NULL, y = NULL) +
    theme_pub(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p_hm, file.path(ASSETS, "fig3_ds_heatmap"), width = 6.5, height = 7)
}

# 4.4 top-DS lollipop(按 |logFC| 排序的显著基因;替代条形图)
lol <- res_pb[res_pb$padj < args$fdr, ]
lol <- lol[order(-abs(lol$logFC)), ]
lol <- head(lol, 20)
if (nrow(lol) >= 2) {
  lol$gene <- factor(lol$gene, levels = rev(lol$gene[order(lol$logFC)]))
  lol$dir  <- ifelse(lol$logFC > 0, "Up", "Down")
  p_lol <- ggplot(lol, aes(logFC, gene, color = dir)) +
    geom_segment(aes(x = 0, xend = logFC, y = gene, yend = gene),
                 color = "grey70", linewidth = 0.5) +
    geom_point(aes(size = -log10(p_val))) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#3C5488"), name = NULL) +
    scale_size_continuous(name = expression(-log[10]~italic(p)), range = c(2, 6)) +
    labs(title = "Top differential-state genes (pseudobulk)",
         subtitle = sprintf("|log2FC|-ranked · FDR<%.2g", args$fdr),
         x = expression(log[2]~fold~change), y = NULL) +
    theme_pub(base_size = 12)
  save_fig(p_lol, file.path(ASSETS, "fig4_top_ds_lollipop"), width = 6.5, height = 6)
}

# 4.5 ★诚实基线对照图:实测假阳性率 dumbbell + 显著基因 p 值分布
if (!is.null(cmp_tbl)) {
  # (a) dumbbell: 两法在 null(供体效应)基因上的实测假阳性率 vs 真 DS 召回
  # 注:列名用 PB / Cell(无连字符),否则 reshape 生成 value.Cell-level 致 aes() 取不到列。
  wide <- data.frame(
    metric = c("False-positive rate\n(null donor-effect genes)", "Recall\n(true DS genes)"),
    PB     = c(cmp_tbl$FalsePositiveRate_nullDonor[1], cmp_tbl$Recall_trueDS[1]),   # pseudobulk
    Cell   = c(cmp_tbl$FalsePositiveRate_nullDonor[2], cmp_tbl$Recall_trueDS[2]),   # cell-level
    stringsAsFactors = FALSE)
  p_db <- ggplot(wide, aes(y = metric)) +
    geom_segment(aes(x = PB, xend = Cell, yend = metric),
                 color = "grey65", linewidth = 1.1) +
    geom_point(aes(x = PB, color = "Pseudobulk"), size = 5) +
    geom_point(aes(x = Cell, color = "Cell-level"), size = 5) +
    geom_vline(xintercept = args$fdr, linetype = "dashed", color = "#E64B35", linewidth = 0.4) +
    annotate("text", x = args$fdr, y = 0.6, label = sprintf("nominal FDR=%.2g", args$fdr),
             color = "#E64B35", size = 3, hjust = -0.05) +
    scale_color_manual(values = c(Pseudobulk = "#00A087", "Cell-level" = "#DC0000"), name = NULL) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = "Honest baseline: pseudo-replication inflates false positives",
         subtitle = "Cell-level DE calls donor-driven null genes 'significant'; pseudobulk does not",
         x = "Proportion of genes called significant", y = NULL) +
    theme_pub(base_size = 11)
  save_fig(p_db, file.path(ASSETS, "fig5_honest_baseline_dumbbell"), width = 7.5, height = 4)

  # (b) raincloud/violin: null 供体基因的 -log10(p) 分布,两法对照(越靠右假信号越强)
  rc <- rbind(
    data.frame(method = "Pseudobulk",
               neglogp = -log10(pmax(res_pb[res_pb$cluster_id == "cluster1" &
                          res_pb$gene %in% truth$null_donor, "p_val"], 1e-300))),
    data.frame(method = "Cell-level",
               neglogp = -log10(pmax(res_cell[res_cell$cluster_id == "cluster1" &
                          res_cell$gene %in% truth$null_donor, "p_val"], 1e-300))))
  p_rc <- ggplot(rc, aes(method, neglogp, fill = method, color = method)) +
    geom_violin(width = 0.85, alpha = 0.35, color = NA, trim = FALSE) +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.6) +
    geom_boxplot(width = 0.12, alpha = 0.7, outlier.shape = NA, color = "black") +
    geom_hline(yintercept = -log10(args$fdr), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = c(Pseudobulk = "#00A087", "Cell-level" = "#DC0000"), guide = "none") +
    scale_color_manual(values = c(Pseudobulk = "#00A087", "Cell-level" = "#DC0000"), guide = "none") +
    labs(title = "Null genes (donor-effect only): significance distribution",
         subtitle = "Cell-level DE pushes true-null genes far past the threshold",
         x = NULL, y = expression(-log[10]~italic(p)~"(null donor-effect genes)")) +
    theme_pub(base_size = 12)
  save_fig(p_rc, file.path(ASSETS, "fig6_null_pvalue_distribution"), width = 6, height = 5)
}

cat("\n完成。结果表见", normalizePath(args$outdir), "; 展示图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖快照(铁律6)
