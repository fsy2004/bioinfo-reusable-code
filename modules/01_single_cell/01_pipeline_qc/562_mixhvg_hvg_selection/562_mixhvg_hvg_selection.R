# =============================================================================
# 562 · mixHVG — 混合多种高变基因(HVG)选择方法
# -----------------------------------------------------------------------------
# 上游   : Zhao R, Lu J, Li Y, Zhou W, Zhao N, Ji H.
#          "A systematic evaluation of highly variable gene selection methods for
#          single-cell RNA-sequencing." Genome Biology 2025;26(1):424.
#          PMID 41382205 · doi:10.1186/s13059-025-03887-x  (已用 NCBI E-utilities 核实)
#          R 包 mixhvg (CRAN 1.0.1) · https://github.com/RuzhangZhao/mixhvg
# -----------------------------------------------------------------------------
# 本模块做什么:
#   ① 基线 (baseline):Seurat vst (= mixhvg 的 "seuratv3"),单一方法,即 Seurat
#      FindVariableFeatures 的默认做法。任何"混合更好"的说法都必须打得过它。
#   ② 逐方法打分:按 mixhvg 源码里各方法真实调用的底层函数复算 feature score。
#   ③ 混合 (comb_rank):把各方法分数各自 rank,逐基因取 max,再取 top-n。
#   ④ 评测:对合成数据的 ground-truth marker 算 recall,并用所选基因跑
#      PCA+kmeans 与真实细胞类型比 ARI(下游可用性,而非只看基因名重合)。
#   ⑤ 若本机装了 mixhvg,额外调用官方 FindVariableFeaturesMix() 并报告与本模块
#      本地实现的 Jaccard 一致性(守卫式:没装就跳过,绝不假装跑过)。
#
# API 出处(逐字读自上游源码,未臆造):
#   签名  https://github.com/RuzhangZhao/mixhvg/blob/HEAD/man/FindVariableFeaturesMix.Rd
#   实现  https://raw.githubusercontent.com/RuzhangZhao/mixhvg/HEAD/R/FindVariableFeaturesMix.R
#   FindVariableFeaturesMix(object, method.names = c("scran","scran_pos","seuratv1"),
#       nfeatures = 2000, loess.span = 0.3, clip.max = "auto", num.bin = 20,
#       binning.method = "equal_width", extra.rank = NULL, verbose = FALSE)
#   注:上游 README 页写的默认组合是 c("scran","seuratv1","mv_PFlogPF","scran_pos"),
#       与 .Rd / 源码里的 c("scran","scran_pos","seuratv1") 不一致 —— 本模块以源码为准,
#       并在 README 中如实记录该出入。
#
# 图中文字英文,代码注释中文。禁止条形图(用 lollipop / slopegraph / 散点 / heatmap)。
# =============================================================================

suppressWarnings(suppressMessages({
  library(Matrix)
}))

# ---- 定位脚本目录 + 载入统一出图框架 ----------------------------------------
.this_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
})
.fw <- normalizePath(file.path(.this_dir, "..", "..", "..", "_framework", "theme_pub.R"),
                     mustWork = FALSE)
if (!file.exists(.fw)) stop("找不到框架文件 theme_pub.R:", .fw)
source(.fw)

set.seed(562)   # 固定随机种子(kmeans / PCA 可复现)

# ---- 参数区(默认指向 example_data/,支持 --key value 覆盖)------------------
args <- bio_args(list(
  counts   = file.path(.this_dir, "example_data", "counts.csv"),
  meta     = file.path(.this_dir, "example_data", "cell_metadata.csv"),
  truth    = file.path(.this_dir, "example_data", "ground_truth_hvg.csv"),
  outdir   = file.path(.this_dir, "results"),
  assets   = file.path(.this_dir, "assets"),
  nfeatures = "100",                                  # 示例数据 700 基因,取 100
  methods  = "scran,scran_pos,seuratv1",              # 混合组合(源码默认)
  baseline = "seuratv3"                               # 基线 = Seurat vst
))
NFEAT <- as.integer(args$nfeatures)
# 上游源码 FindVariableFeaturesMix.R:341-343 的三个别名映射,原样照搬
.alias <- function(v) {
  v[v == "disp_nc"]  <- "seuratv1"
  v[v == "logmv_ct"] <- "seuratv3"
  v[v == "mv_lognc"] <- "scran"
  unique(v)
}
MIX_METHODS <- .alias(trimws(strsplit(args$methods, ",")[[1]]))
BASELINE <- .alias(args$baseline)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$assets, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Step 1 · 读入表达矩阵
# =============================================================================
cat("Step 1 · 读入 counts 矩阵\n")
cnt_df <- read_table_smart(args$counts, row_names = TRUE)
counts <- as(as.matrix(cnt_df), "dgCMatrix")
storage.mode(counts@x) <- "double"
cat(sprintf("       genes = %d, cells = %d\n", nrow(counts), ncol(counts)))
if (nrow(counts) <= NFEAT)
  stop(sprintf("nfeatures(%d) 必须小于基因数(%d)", NFEAT, nrow(counts)))

has_meta  <- file.exists(args$meta)
has_truth <- file.exists(args$truth)
labels <- NULL; truth_genes <- NULL; truth_tier <- NULL
if (has_meta) {
  m <- read_table_smart(args$meta)
  labels <- setNames(as.character(m[[2]]), as.character(m[[1]]))[colnames(counts)]
  cat(sprintf("       cell types: %s\n", paste(names(table(labels)), collapse = ", ")))
}
if (has_truth) {
  t_df <- read_table_smart(args$truth)
  truth_genes <- as.character(t_df[[1]])
  truth_tier  <- if (ncol(t_df) >= 2) setNames(as.character(t_df[[2]]), truth_genes) else NULL
  truth_genes <- intersect(truth_genes, rownames(counts))
  cat(sprintf("       ground-truth HVG: %d\n", length(truth_genes)))
}

# =============================================================================
# Step 2 · 预处理出 mixhvg 各方法所需的四种矩阵
#   (变换方式逐行照抄上游源码 FindVariableFeaturesMix.R,不自行发明)
# =============================================================================
cat("Step 2 · 构造 counts / normalized / lognormalized / PFlog1pPF 四套矩阵\n")
suppressWarnings(suppressMessages(library(Seurat)))

lognormalizedcounts <- Seurat::NormalizeData(counts, verbose = FALSE)
normalizedcounts <- lognormalizedcounts
normalizedcounts@x <- exp(normalizedcounts@x) - 1     # 源码即如此反 log 得到 nc
# PFlog1pPF:proportional fitting → log1p → proportional fitting(源码写法)
PFlog1pPF <- Seurat::NormalizeData(counts, scale.factor = mean(Matrix::colSums(counts)),
                                   verbose = FALSE)
PFlog1pPF <- Seurat::NormalizeData(PFlog1pPF, scale.factor = mean(Matrix::colSums(PFlog1pPF)),
                                   normalization.method = "RC", verbose = FALSE)

HAS_SCRAN <- requireNamespace("scran", quietly = TRUE) &&
             requireNamespace("SingleCellExperiment", quietly = TRUE) &&
             requireNamespace("SummarizedExperiment", quietly = TRUE) &&
             requireNamespace("scuttle", quietly = TRUE)
if (!HAS_SCRAN)
  cat("       [guard] scran/SingleCellExperiment/scuttle 缺失 → scran 系方法将跳过\n")

# =============================================================================
# Step 3 · 逐方法 feature score
#   每个 case 的底层调用与参数,均取自上游源码 FindFeatureVal()
# =============================================================================
feature_val <- function(method) {
  vst <- function(mat) Seurat::FindVariableFeatures(
    mat, loess.span = 0.3, clip.max = "auto", num.bin = 20,
    binning.method = "equal_width", verbose = FALSE)$vst.variance.standardized
  disp <- function(mat) {
    v <- Seurat::FindVariableFeatures(mat, selection.method = "disp", num.bin = 20,
                                      binning.method = "equal_width",
                                      verbose = FALSE)$mvp.dispersion
    v[is.na(v)] <- 0; v[v < 0] <- 0; v
  }
  # scran 系:取 bio 分量,负值/NA 截 0(上游 dec.keep 逻辑)
  #   mode = "asis"    → 把传入矩阵直接当 logcounts(上游 mv_ct/mv_nc/mv_PFlogPF 分支:
  #                      sce@assays@data$logcounts <- sce@assays@data$counts)
  #   mode = "poisson" → modelGeneVarByPoisson,直接吃 counts(上游 scran_pos 分支)
  #   mode = "lognorm" → scuttle::logNormCounts 后建模(上游 "scran" 分支在 counts
  #                      非空时走的就是这条,不是 Seurat 的 NormalizeData)
  mv <- function(mat, mode = c("asis", "poisson", "lognorm")) {
    mode <- match.arg(mode)
    if (!HAS_SCRAN) return(NULL)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
    dec <- switch(mode,
      "poisson" = scran::modelGeneVarByPoisson(sce),
      "lognorm" = scran::modelGeneVar(scuttle::logNormCounts(sce)),
      "asis"    = {
        SummarizedExperiment::assay(sce, "logcounts") <- mat
        scran::modelGeneVar(sce)
      })
    v <- dec$bio; v[is.na(v) | v <= 0] <- 0; v
  }
  switch(method,
    "seuratv3"      = vst(counts),
    "logmv_nc"      = vst(normalizedcounts),
    "logmv_lognc"   = vst(lognormalizedcounts),
    "logmv_PFlogPF" = vst(PFlog1pPF),
    "seuratv1"      = disp(lognormalizedcounts),
    "disp_ct"       = disp(log1p(counts)),
    "disp_lognc"    = disp(log1p(lognormalizedcounts)),
    "disp_PFlogPF"  = disp(log1p(PFlog1pPF)),
    "mv_ct"         = mv(counts,          "asis"),
    "mv_nc"         = mv(normalizedcounts,"asis"),
    "mv_PFlogPF"    = mv(PFlog1pPF,       "asis"),
    # 上游 "scran" 分支:counts 非空 → SingleCellExperiment + scuttle::logNormCounts
    # + modelGeneVar。本模块传的就是 counts,故走 logNormCounts(而非 Seurat NormalizeData)
    "scran"         = mv(counts,          "lognorm"),
    "scran_pos"     = mv(counts,          "poisson"),
    "mean_max_ct"     = Matrix::rowMeans(counts),
    "mean_max_nc"     = Matrix::rowMeans(normalizedcounts),
    "mean_max_lognc"  = Matrix::rowMeans(lognormalizedcounts),
    "mean_max_PFlogPF"= Matrix::rowMeans(PFlog1pPF),
    stop("未知方法: ", method)
  )
}

ALL_METHODS <- unique(c(BASELINE, MIX_METHODS,
                        "logmv_PFlogPF", "disp_PFlogPF", "mean_max_lognc"))
cat("Step 3 · 逐方法打分:", paste(ALL_METHODS, collapse = ", "), "\n")
scores <- list()
for (mth in ALL_METHODS) {
  v <- tryCatch(feature_val(mth), error = function(e) {
    cat(sprintf("       [skip] %s 失败: %s\n", mth, conditionMessage(e))); NULL })
  if (is.null(v)) next
  names(v) <- rownames(counts)
  scores[[mth]] <- v
  cat(sprintf("       %-16s ok\n", mth))
}
if (!BASELINE %in% names(scores)) stop("基线方法未能计算,无法比较")

# ---- comb_rank:上游混合逻辑(rank 各方法,逐基因取 max)---------------------
comb_rank <- function(lst) {
  R <- vapply(lst, function(v) rank(v, ties.method = "min"), numeric(length(lst[[1]])))
  setNames(apply(R, 1, max), names(lst[[1]]))
}
top_n <- function(v, n) names(sort(v, decreasing = TRUE))[seq_len(n)]

mix_ok <- intersect(MIX_METHODS, names(scores))
if (length(mix_ok) < 2)
  cat("       [guard] 可用方法不足 2 个,混合退化为单方法\n")
mix_score <- if (length(mix_ok) >= 2) comb_rank(scores[mix_ok]) else scores[[mix_ok[1]]]

sel <- lapply(scores, top_n, n = NFEAT)
sel[["MIX"]] <- top_n(mix_score, NFEAT)
cat(sprintf("Step 4 · 各方法各选 top-%d;混合组合 = %s\n", NFEAT,
            paste(mix_ok, collapse = " + ")))

# =============================================================================
# Step 5 · 评测:ground-truth recall + 下游聚类 ARI
# =============================================================================
cat("Step 5 · 评测(recall + 下游 kmeans ARI)\n")
logmat <- as.matrix(lognormalizedcounts)

# 平均轮廓宽度(以真实细胞类型为簇);越高说明所选基因保留的类型结构越干净
mean_silhouette <- function(pcs, lab) {
  D <- as.matrix(dist(pcs)); lab <- as.character(lab); u <- unique(lab)
  if (length(u) < 2) return(NA_real_)
  s <- vapply(seq_along(lab), function(i) {
    same <- which(lab == lab[i] & seq_along(lab) != i)
    if (!length(same)) return(0)
    a <- mean(D[i, same])
    b <- min(vapply(setdiff(u, lab[i]), function(g) mean(D[i, lab == g]), numeric(1)))
    (b - a) / max(a, b)
  }, numeric(1))
  mean(s)
}

eval_one <- function(genes) {
  out <- list()
  if (!is.null(truth_genes) && length(truth_genes)) {
    out$recall <- length(intersect(genes, truth_genes)) / length(truth_genes)
    if (!is.null(truth_tier)) {
      for (tt in unique(truth_tier)) {
        g_t <- names(truth_tier)[truth_tier == tt]
        g_t <- intersect(g_t, rownames(counts))
        out[[paste0("recall_", tt)]] <- length(intersect(genes, g_t)) / max(1, length(g_t))
      }
    }
  }
  if (!is.null(labels)) {
    X <- t(scale(t(logmat[genes, , drop = FALSE])))     # 基因内 z-score
    X[!is.finite(X)] <- 0
    pcs <- prcomp(t(X), center = FALSE, scale. = FALSE)$x[, seq_len(min(15, length(genes) - 1)), drop = FALSE]
    k <- length(unique(labels))
    km <- kmeans(pcs, centers = k, nstart = 25, iter.max = 100)
    out$ARI <- if (requireNamespace("mclust", quietly = TRUE))
      mclust::adjustedRandIndex(km$cluster, labels) else NA_real_
    # ARI 在容易的数据上会饱和(全 1),分辨不出方法优劣;补一个连续量:
    # 用真实标签算平均轮廓宽度 = 所选基因把已知细胞类型分开的程度(手写,不加依赖)
    out$silhouette <- mean_silhouette(pcs, labels)
  }
  out
}

res <- do.call(rbind, lapply(names(sel), function(m) {
  e <- eval_one(sel[[m]])
  # 无 meta 且无 truth 时 e 是空 list,as.data.frame(空) 是 0 行,
  # 与 method 列的 1 行拼不起来 → 只留 method 列(选基因清单照常落盘)
  if (!length(e)) return(data.frame(method = m, check.names = FALSE))
  data.frame(method = m, as.data.frame(e), check.names = FALSE)
}))
res$is_mix <- res$method == "MIX"
# 只给 counts(无 meta / 无 truth)时 recall / silhouette 两列都不存在,
# 此时不排序 —— 否则 order() 拿到 0 长度向量会把整张表悄悄清空。
.sort_keys <- Filter(function(k) k %in% names(res), c("recall", "silhouette"))
if (length(.sort_keys))
  res <- res[do.call(order, lapply(.sort_keys, function(k) -res[[k]])), ]
print(res, row.names = FALSE)
write.csv(res, file.path(args$outdir, "562_method_metrics.csv"), row.names = FALSE)

# 选中基因清单落盘
sel_df <- do.call(rbind, lapply(names(sel), function(m)
  data.frame(method = m, rank = seq_along(sel[[m]]), gene = sel[[m]])))
write.csv(sel_df, file.path(args$outdir, "562_selected_features.csv"), row.names = FALSE)

# =============================================================================
# Step 6 · 守卫式调用官方 mixhvg(未安装则跳过,不伪造结果)
# =============================================================================
cat("Step 6 · 官方 mixhvg 包核对\n")
upstream <- list(status = "skipped")
if (requireNamespace("mixhvg", quietly = TRUE)) {
  up <- tryCatch({
    mixhvg::FindVariableFeaturesMix(counts, method.names = MIX_METHODS,
                                    nfeatures = NFEAT, verbose = FALSE)
  }, error = function(e) { cat("       官方调用报错:", conditionMessage(e), "\n"); NULL })
  if (!is.null(up)) {
    j <- length(intersect(up, sel$MIX)) / length(union(up, sel$MIX))
    upstream <- list(status = "ok", jaccard_vs_local = round(j, 3),
                     version = as.character(utils::packageVersion("mixhvg")))
    cat(sprintf("       官方 vs 本地实现 Jaccard = %.3f\n", j))
    writeLines(up, file.path(args$outdir, "562_upstream_mixhvg_features.txt"))
  }
} else {
  upstream$reason <- "mixhvg 未安装;install.packages('mixhvg')"
  cat("       [guard] mixhvg 未安装 → 使用本模块按上游源码复现的 comb_rank 混合逻辑\n")
  cat("       安装后可复跑本步做一致性核对:install.packages('mixhvg')\n")
}

# =============================================================================
# Step 7 · 出图(lollipop / slopegraph / 散点 / heatmap;不用条形图)
# =============================================================================
cat("Step 7 · 出图\n")
theme_set(theme_pub(base_size = 11))
PLOT <- function(p, stem, w, h) {
  save_fig(p, file.path(args$outdir, stem), width = w, height = h)  # 矢量 PDF + PNG
  # assets/ 只放 README 预览用的 PNG(库规范:assets 提交,results 不提交)
  invisible(file.copy(file.path(args$outdir, paste0(stem, ".png")),
                      file.path(args$assets, paste0(stem, ".png")), overwrite = TRUE))
  cat(sprintf("       fig -> %s.png\n", stem))
}

# --- 图1 lollipop:各方法 ground-truth recall,MIX 高亮 -----------------------
if ("recall" %in% names(res)) {
  d1 <- res[, c("method", "recall", "is_mix")]
  d1$method <- factor(d1$method, levels = d1$method[order(d1$recall)])
  p1 <- ggplot(d1, aes(x = recall, y = method)) +
    geom_segment(aes(x = 0, xend = recall, yend = method), colour = "grey75", linewidth = 0.6) +
    geom_point(aes(colour = is_mix), size = 4) +
    scale_colour_manual(values = c(`FALSE` = "#4DBBD5", `TRUE` = "#E64B35"),
                        labels = c("single method", "mixture"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = sprintf("Recall of ground-truth HVGs (top %d)", NFEAT), y = NULL,
         title = "Single-method vs mixture HVG recall")
  PLOT(p1, "562_fig1_recall_lollipop", 6.6, 4.6)
}

# --- 图2 slopegraph:每档 marker (HI/MID/LO) 上 baseline → MIX 的变化 ---------
tier_cols <- grep("^recall_", names(res), value = TRUE)
if (length(tier_cols) && BASELINE %in% res$method) {
  d2 <- do.call(rbind, lapply(tier_cols, function(cc) data.frame(
    tier = sub("^recall_", "", cc),
    baseline = res[[cc]][res$method == BASELINE],
    mixture  = res[[cc]][res$method == "MIX"])))
  d2l <- rbind(data.frame(tier = d2$tier, stage = "baseline\n(Seurat vst)", value = d2$baseline),
               data.frame(tier = d2$tier, stage = "mixture\n(comb_rank)",  value = d2$mixture))
  d2l$stage <- factor(d2l$stage, levels = c("baseline\n(Seurat vst)", "mixture\n(comb_rank)"))
  p2 <- ggplot(d2l, aes(x = stage, y = value, group = tier, colour = tier)) +
    geom_line(linewidth = 1.1) + geom_point(size = 3.4) +
    ggrepel::geom_text_repel(data = subset(d2l, stage == levels(d2l$stage)[2]),
                             aes(label = tier), nudge_x = 0.12, size = 3.4,
                             direction = "y", segment.colour = NA, show.legend = FALSE) +
    scale_colour_manual(values = pal_pub(length(unique(d2l$tier)), "npg"), guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Recall within expression tier",
         title = "Where the mixture gains: recall by marker expression tier",
         subtitle = "HI / MID / LO = high / medium / low mean expression markers")
  PLOT(p2, "562_fig2_tier_slopegraph", 5.4, 5.0)
}

# --- 图3 散点:mean-variance 平面,标出 baseline / MIX 各自独有的选择 ---------
gmean <- Matrix::rowMeans(lognormalizedcounts)
gvar  <- apply(as.matrix(lognormalizedcounts), 1, var)
only_mix  <- setdiff(sel$MIX, sel[[BASELINE]])
only_base <- setdiff(sel[[BASELINE]], sel$MIX)
d3 <- data.frame(gene = rownames(counts), mean = gmean, var = gvar)
d3$grp <- "not selected"
d3$grp[d3$gene %in% intersect(sel$MIX, sel[[BASELINE]])] <- "both"
d3$grp[d3$gene %in% only_base] <- "baseline only"
d3$grp[d3$gene %in% only_mix]  <- "mixture only"
d3$grp <- factor(d3$grp, levels = c("not selected", "both", "baseline only", "mixture only"))
d3 <- d3[order(d3$grp), ]
p3 <- ggplot(d3, aes(x = mean, y = var, colour = grp, size = grp, alpha = grp)) +
  geom_point() +
  scale_colour_manual(values = c("not selected" = "grey82", "both" = "#3C5488",
                                 "baseline only" = "#4DBBD5", "mixture only" = "#E64B35"),
                      name = NULL) +
  scale_size_manual(values = c(0.7, 1.3, 2.1, 2.1), guide = "none") +
  scale_alpha_manual(values = c(0.5, 0.7, 0.95, 0.95), guide = "none") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Mean expression (lognorm, log10)", y = "Variance (log10)",
       title = "Which genes the mixture adds over the vst baseline")
PLOT(p3, "562_fig3_mean_variance_scatter", 6.6, 4.8)

# --- 图4 heatmap:方法两两 Jaccard 重合度 ------------------------------------
mm <- names(sel)
J <- outer(seq_along(mm), seq_along(mm), Vectorize(function(i, j)
  length(intersect(sel[[mm[i]]], sel[[mm[j]]])) / length(union(sel[[mm[i]]], sel[[mm[j]]]))))
dimnames(J) <- list(mm, mm)
d4 <- as.data.frame(as.table(J)); names(d4) <- c("m1", "m2", "jaccard")
ord <- hclust(as.dist(1 - J))$order
d4$m1 <- factor(d4$m1, levels = mm[ord]); d4$m2 <- factor(d4$m2, levels = mm[ord])
p4 <- ggplot(d4, aes(m1, m2, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  # viridis 低值端为深紫、高值端为亮黄 → 深底配白字、亮底配深字
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.7,
            colour = ifelse(d4$jaccard > 0.62, "grey10", "white")) +
  scale_fill_cont(name = "Jaccard", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "Method agreement on selected features",
       subtitle = "Low off-diagonal overlap is why mixing helps") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PLOT(p4, "562_fig4_method_jaccard_heatmap", 6.4, 5.6)

# --- 图5 散点:recall 与下游结构保留 (silhouette) 的联合表现 -----------------
if (all(c("recall", "silhouette") %in% names(res))) {
  p5 <- ggplot(res, aes(x = recall, y = silhouette)) +
    geom_point(aes(colour = is_mix), size = 4) +
    ggrepel::geom_text_repel(aes(label = method), size = 3.2, max.overlaps = 20,
                             box.padding = 0.45, segment.colour = "grey65") +
    scale_colour_manual(values = c(`FALSE` = "#4DBBD5", `TRUE` = "#E64B35"),
                        labels = c("single method", "mixture"), name = NULL) +
    labs(x = "Recall of ground-truth HVGs", y = "Mean silhouette (true cell types)",
         title = "Feature recall vs downstream structure preservation",
         subtitle = "Upper right = recovers true HVGs and keeps cell types separable")
  PLOT(p5, "562_fig5_recall_vs_silhouette", 6.2, 4.8)
}

# =============================================================================
# Step 8 · 汇总
# =============================================================================
summ <- list(
  n_genes = nrow(counts), n_cells = ncol(counts), nfeatures = NFEAT,
  baseline = BASELINE, mixture_methods = mix_ok,
  baseline_recall = if ("recall" %in% names(res)) res$recall[res$method == BASELINE] else NA,
  mixture_recall  = if ("recall" %in% names(res)) res$recall[res$method == "MIX"] else NA,
  baseline_ARI = res$ARI[res$method == BASELINE], mixture_ARI = res$ARI[res$method == "MIX"],
  baseline_silhouette = res$silhouette[res$method == BASELINE],
  mixture_silhouette  = res$silhouette[res$method == "MIX"],
  upstream_mixhvg = upstream)
writeLines(paste(names(unlist(summ)), unlist(summ), sep = ": "),
           file.path(args$outdir, "562_summary.txt"))
# 依赖版本快照(可复现:结果对不上时先查这里的包版本)
capture.output(sessionInfo(), file = file.path(args$outdir, "562_sessionInfo.txt"))
cat("\n[562] done ->", args$outdir, "\n")
# 缺 meta / truth 时对应指标压根没算,打印成 "n/a" 而不是留空或伪造一个数
.fmt <- function(x) if (length(x) != 1 || !is.finite(x)) "n/a" else sprintf("%.3f", x)
cat(sprintf("[562] baseline(%s) recall=%s ARI=%s sil=%s\n", BASELINE,
            .fmt(summ$baseline_recall), .fmt(summ$baseline_ARI), .fmt(summ$baseline_silhouette)))
cat(sprintf("[562] mixture (%s) recall=%s ARI=%s sil=%s\n", paste(mix_ok, collapse = "+"),
            .fmt(summ$mixture_recall), .fmt(summ$mixture_ARI), .fmt(summ$mixture_silhouette)))
