# =============================================================================
# 566 · Phi-Space (PhiSpace) —— 连续表型注释 / soft cell-type annotation
# -----------------------------------------------------------------------------
# 用途   : 把 query 细胞投影到由参考集细胞类型张成的「表型空间」,每个细胞对每个
#          参考类型给出一个**连续得分**,而不是一个硬标签。中间态/过渡态细胞因此
#          不再被强行按到某一个类别上。
# 上游   : Mao J, Deng Y, Lê Cao KA. Φ-Space: continuous phenotyping of single-cell
#          multi-omics data. Genome Biology 2025;26(1):323.
#          doi:10.1186/s13059-025-03755-8 · PMID 41029411   (已用 NCBI esummary 核实)
#          仓库 https://github.com/jiadongm/PhiSpace · 文档 https://jiadongm.github.io/PhiSpace/
# 结构   : 本脚本始终跑一条**本机可跑的朴素基线**(质心相关 + PCA 回归软注释),
#          PhiSpace 本体走守卫式封装:装了才调,没装就打印真实安装命令并跳过。
#          基线存在的意义:任何"连续注释更好"的说法都必须有对照,不能孤零零地报。
# 依赖   : 必需 ggplot2(框架 theme_pub.R);可选 PhiSpace + SingleCellExperiment
# 运行   : Rscript 566_phispace_soft_annotation.R
#          Rscript 566_phispace_soft_annotation.R --ref_expr a.csv --query_expr b.csv
#          Rscript 566_phispace_soft_annotation.R --run_phispace   # 需先装 PhiSpace
# =============================================================================

suppressWarnings(suppressMessages({
  set.seed(566)
}))

# ---- 定位脚本目录 + 载入框架绘图样式 ----------------------------------------
.this_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", a[m[1]]))))
  getwd()
}
HERE <- .this_dir()
source(file.path(HERE, "..", "..", "..", "_framework", "theme_pub.R"))

# ---- 参数区(全部支持 --key value 覆盖)------------------------------------
opt <- bio_args(list(
  ref_expr    = file.path(HERE, "example_data", "reference_counts.csv"),
  ref_meta    = file.path(HERE, "example_data", "reference_metadata.csv"),
  query_expr  = file.path(HERE, "example_data", "query_counts.csv"),
  query_meta  = file.path(HERE, "example_data", "query_metadata.csv"),
  label_col   = "cell_type",   # 参考集 metadata 里的类型列名
  ncomp       = 10,            # PCA/PLS 成分数
  outdir      = file.path(HERE, "results"),
  assets      = file.path(HERE, "assets"),
  run_phispace = FALSE
))
opt$ncomp <- as.integer(opt$ncomp)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$assets, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Step 1  读数据 + 归一化
# =============================================================================
cat("Step 1  读取表达矩阵与 metadata\n")

read_expr <- function(path) {
  df <- read_table_smart(path)
  g  <- as.character(df[[1]])
  m  <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- make.unique(g)
  m                                  # genes x cells
}

ref_x <- read_expr(opt$ref_expr)
qry_x <- read_expr(opt$query_expr)
ref_m <- read_table_smart(opt$ref_meta)
qry_m <- if (file.exists(opt$query_meta)) read_table_smart(opt$query_meta) else NULL

# 只保留两边共有的基因(对应上游的 KeepCommonGenes 思路)
common <- intersect(rownames(ref_x), rownames(qry_x))
if (length(common) < 10) stop("参考集与 query 的共有基因不足 10 个,检查基因命名是否一致")
ref_x <- ref_x[common, , drop = FALSE]
qry_x <- qry_x[common, , drop = FALSE]
cat(sprintf("        共有基因 %d · 参考细胞 %d · query 细胞 %d\n",
            length(common), ncol(ref_x), ncol(qry_x)))

# 参考集标签对齐
stopifnot(opt$label_col %in% colnames(ref_m))
rownames(ref_m) <- as.character(ref_m[[1]])
ref_lab <- factor(ref_m[colnames(ref_x), opt$label_col])
types   <- levels(ref_lab)
cat(sprintf("        参考类型 %d: %s\n", length(types), paste(types, collapse = ", ")))

# CPM + log1p(两侧同一套流程,避免参考/query 处理不一致造成的假差异)
norm_log <- function(m) {
  cs <- colSums(m); cs[cs == 0] <- 1
  log1p(sweep(m, 2, cs, "/") * 1e4)
}
ref_l <- norm_log(ref_x)
qry_l <- norm_log(qry_x)

# =============================================================================
# Step 2  基线 A —— 参考类型质心的 Spearman 相关(最朴素的软打分)
# =============================================================================
cat("Step 2  基线 A: 质心相关软打分 (centroid Spearman correlation)\n")

centroids <- sapply(types, function(k) rowMeans(ref_l[, ref_lab == k, drop = FALSE]))
score_cor <- cor(qry_l, centroids, method = "spearman")   # cells x types
colnames(score_cor) <- types

# =============================================================================
# Step 3  基线 B —— PCA 回归软注释(对应上游 regMethod="PCA" 的朴素同构版)
#   dummy-code 参考标签 Y → 参考集 PCA → OLS 回归 → 把 query 投到同一组载荷上预测。
#   这是我们自己写的朴素对照,**不是** PhiSpace 的实现,不要当成复现 PhiSpace。
# =============================================================================
cat(sprintf("Step 3  基线 B: PCA 回归软注释 (ncomp=%d)\n", opt$ncomp))

Y <- model.matrix(~ 0 + ref_lab); colnames(Y) <- types      # dummy 编码
gene_mu <- rowMeans(ref_l)
Xr <- t(ref_l - gene_mu)                                    # cells x genes,基因中心化
Xq <- t(qry_l - gene_mu)                                    # query 用**参考集**均值中心化

nc <- min(opt$ncomp, ncol(Xr) - 1, nrow(Xr) - 1)
sv <- svd(Xr, nu = nc, nv = nc)
Tr <- sv$u %*% diag(sv$d[seq_len(nc)], nc, nc)              # 参考集得分
Yc <- scale(Y, center = TRUE, scale = FALSE)
B  <- qr.solve(crossprod(Tr) + diag(1e-8, nc), crossprod(Tr, Yc))
Tq <- Xq %*% sv$v                                           # query 投影到同一载荷
score_pca <- sweep(Tq %*% B, 2, attr(Yc, "scaled:center"), "+")
colnames(score_pca) <- types

# 归一化成「组成比例」以便和真实比例比较:截断到非负后按行归一
to_comp <- function(s) {
  s <- pmax(s, 0); rs <- rowSums(s); rs[rs == 0] <- 1; s / rs
}
comp_cor <- to_comp(score_cor)
comp_pca <- to_comp(score_pca)

# 硬标签(argmax)—— 用来展示硬标签把过渡态抹掉了什么
hard_lab <- types[max.col(score_pca, ties.method = "first")]

methods <- list(`Centroid corr` = comp_cor, `PCA regression` = comp_pca)

# =============================================================================
# Step 4  PhiSpace 本体(守卫式:未安装则跳过,绝不静默降级)
# =============================================================================
cat("Step 4  PhiSpace 本体\n")

run_phispace <- function() {
  # 上游 pkg/DESCRIPTION 写死 Depends: R (>= 4.5.0);本机 R 4.4.3 装不上,先讲清楚
  if (getRversion() < "4.5.0") {
    return(list(status = "skipped", reason = sprintf(
      "上游 DESCRIPTION 要求 R >= 4.5.0,本机为 R %s,PhiSpace 在此环境无法安装。",
      getRversion())))
  }
  need <- c("PhiSpace", "SingleCellExperiment", "S4Vectors")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    return(list(status = "skipped", reason = sprintf(
      "未安装: %s。安装: BiocManager::install('jiadongm/PhiSpace/pkg')", paste(miss, collapse = ", "))))
  }
  # 以下调用的函数名/参数名/默认值均已对照**克隆下来的上游源码**逐个核对(源码为准,非文档):
  #   PhiSpace()    pkg/R/PhiSpaceR.R:55   —— reference/query/phenotypes/refAssay/regMethod/ncomp/updateRef 均为真实形参
  #   RankTransf()  pkg/R/RankTransf.R:11  —— 形参是 assayname(小写 n)、targetAssay 默认 'rank'
  #   两者均在 pkg/NAMESPACE 里 export
  # ★注意:上游自带的 PhiSpace_Guide_for_VibeCoding.md 把参数写成 assayName(大写 N),与源码不符。
  #   这里用**位置参数**传 "counts",不受该笔误影响。
  # 教程用 counts 起手 → RankTransf() 得到 "rank" assay → PhiSpace(..., refAssay="rank")。
  ref_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ref_x),
    colData = S4Vectors::DataFrame(cell_type = as.character(ref_lab)))
  qry_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = qry_x))
  ref_sce <- PhiSpace::RankTransf(ref_sce, "counts")
  qry_sce <- PhiSpace::RankTransf(qry_sce, "counts")
  # ncomp:上游 pkg/R/PhiSpaceR_1ref.R:213 `if(is.null(ncomp)) ncomp <- ncol(YY)`,
  # 即默认 = 表型总数。本模块的 --ncomp 是给 PCA 基线调的(默认 10),直接塞给
  # PhiSpace 会远超 3 个表型的默认值,语义不是一回事 → 只在 <= 表型数时才传,
  # 否则传 NULL 让上游用自己的默认,不替上游做主。
  phi_ncomp <- if (opt$ncomp <= length(types)) opt$ncomp else NULL
  res <- PhiSpace::PhiSpace(
    reference  = ref_sce,
    query      = qry_sce,
    phenotypes = "cell_type",
    refAssay   = "rank",
    regMethod  = "PLS",
    ncomp      = phi_ncomp,
    updateRef  = FALSE)
  sc <- as.matrix(SingleCellExperiment::reducedDim(res, "PhiSpace"))
  list(status = "ok", scores = sc,
       phispace_version = as.character(utils::packageVersion("PhiSpace")))
}

phi <- if (isTRUE(opt$run_phispace) || identical(opt$run_phispace, "TRUE")) {
  tryCatch(run_phispace(), error = function(e)
    list(status = "error", reason = conditionMessage(e)))
} else list(status = "not requested", reason = "加 --run_phispace 开启(需先安装 PhiSpace)")

cat(sprintf("        status: %s\n", phi$status))
if (!is.null(phi$reason)) cat(sprintf("        %s\n", phi$reason))
if (identical(phi$status, "ok")) {
  sc <- phi$scores
  # 只有列名能和参考类型对上时才纳入比较,避免张冠李戴
  if (all(types %in% colnames(sc))) {
    methods[["PhiSpace (PLS)"]] <- to_comp(sc[, types, drop = FALSE])
  } else {
    cat(sprintf("        PhiSpace 返回列名 = %s,与参考类型不完全一致,已单独落盘不参与并列比较\n",
                paste(colnames(sc), collapse = ",")))
  }
  write.csv(sc, file.path(opt$outdir, "phispace_scores.csv"))
}

# =============================================================================
# Step 5  评估:连续得分能否还原真实混合比例
# =============================================================================
cat("Step 5  评估与落盘\n")

truth_cols <- paste0("true_", types)
has_truth  <- !is.null(qry_m) && all(truth_cols %in% colnames(qry_m))
eval_df <- NULL

if (has_truth) {
  rownames(qry_m) <- as.character(qry_m[[1]])
  qm <- qry_m[colnames(qry_x), , drop = FALSE]
  truth <- as.matrix(qm[, truth_cols, drop = FALSE]); colnames(truth) <- types

  eval_df <- do.call(rbind, lapply(names(methods), function(nm) {
    s <- methods[[nm]]
    do.call(rbind, lapply(types, function(k) data.frame(
      method = nm, cell_type = k,
      pearson  = cor(s[, k], truth[, k]),
      spearman = cor(s[, k], truth[, k], method = "spearman"),
      rmse     = sqrt(mean((s[, k] - truth[, k])^2)))))
  }))
  write.csv(eval_df, file.path(opt$outdir, "score_vs_truth_metrics.csv"), row.names = FALSE)
  print(eval_df, row.names = FALSE, digits = 3)

  # 硬标签在过渡态上必然二选一 —— 量化这件事
  is_mid <- qm$true_group == "Intermediate_AB"
  cat(sprintf("        过渡态细胞 %d 个,硬标签分配: %s\n", sum(is_mid),
              paste(sprintf("%s=%d", names(table(hard_lab[is_mid])),
                            as.integer(table(hard_lab[is_mid]))), collapse = " ")))
  cat(sprintf("        纯细胞硬标签准确率: %.3f\n",
              mean(hard_lab[!is_mid] == sub("^Pure_", "Type", qm$true_group[!is_mid]))))
} else {
  cat("        query_meta 无 true_* 列,跳过定量评估(真实数据的常态)\n")
}

for (nm in names(methods)) {
  f <- file.path(opt$outdir, sprintf("soft_scores_%s.csv",
                                     gsub("[^A-Za-z0-9]+", "_", tolower(nm))))
  write.csv(data.frame(cell = colnames(qry_x), methods[[nm]], check.names = FALSE),
            f, row.names = FALSE)
}
write.csv(data.frame(cell = colnames(qry_x), hard_label = hard_lab),
          file.path(opt$outdir, "hard_labels_argmax.csv"), row.names = FALSE)

# =============================================================================
# Step 6  出图(框架样式;无条形图)
# =============================================================================
cat("Step 6  出图\n")

ord <- if (has_truth) order(qm$true_group != "Intermediate_AB", qm$true_TypeB) else seq_len(ncol(qry_x))
cells_ord <- colnames(qry_x)[ord]

# --- Fig 1  软得分热图:细胞按真实混合比例排序,过渡态呈现平滑梯度 ---
hm <- do.call(rbind, lapply(types, function(k) data.frame(
  cell = factor(cells_ord, levels = cells_ord), cell_type = k,
  score = comp_pca[ord, k])))
p1 <- ggplot(hm, aes(cell, cell_type, fill = score)) +
  geom_tile() +
  scale_fill_cont(name = "Soft score") +
  labs(title = "Continuous phenotype scores across query cells",
       subtitle = "Cells ordered by true TypeB fraction; transitional cells form a gradient, not a block",
       x = "Query cells (ordered)", y = NULL) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p1, file.path(opt$assets, "fig1_soft_score_heatmap"), width = 8.4, height = 3.2)

if (has_truth) {
  # --- Fig 2  散点:真实混合比例 vs 预测得分,按方法分面 ---
  sc_df <- do.call(rbind, lapply(names(methods), function(nm) data.frame(
    method = nm, truth = truth[, "TypeB"], score = methods[[nm]][, "TypeB"],
    group = qm$true_group)))
  p2 <- ggplot(sc_df, aes(truth, score, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(size = 1.5, alpha = 0.75) +
    facet_wrap(~ method) +
    scale_color_pub("npg") +
    labs(colour = NULL, title = "Recovery of the true TypeB fraction",
         subtitle = "Dashed line = identity; a soft annotation should track the whole continuum",
         x = "True TypeB fraction", y = "Predicted TypeB score") +
    theme_pub(base_size = 11, border = TRUE, legend = "bottom")
  save_fig(p2, file.path(opt$assets, "fig2_truth_vs_score_scatter"), width = 8.0, height = 4.0)

  # --- Fig 3  小提琴 + 抖动点:过渡态坐落在两个纯群之间 ---
  vd <- data.frame(group = factor(qm$true_group,
                                  levels = c("Pure_A", "Intermediate_AB", "Pure_B", "Pure_C")),
                   score = comp_pca[, "TypeB"])
  p3 <- ggplot(vd, aes(group, score, fill = group)) +
    geom_violin(width = 0.85, alpha = 0.5, colour = NA, trim = FALSE) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.55, colour = "grey20") +
    stat_summary(fun = median, geom = "point", size = 2.4, shape = 21,
                 fill = "white", colour = "black") +
    scale_fill_pub("npg", guide = "none") +
    labs(title = "What a hard label would erase",
         subtitle = "TypeB soft score by true population (PCA-regression baseline)",
         x = NULL, y = "TypeB soft score") +
    theme_pub(base_size = 11)
  save_fig(p3, file.path(opt$assets, "fig3_score_by_population_violin"), width = 5.6, height = 4.2)

  # --- Fig 4  棒棒糖图:各方法 × 各类型 与真值的相关(明确禁用条形图) ---
  ed <- eval_df[order(eval_df$pearson), ]
  ed$lab <- factor(sprintf("%s · %s", ed$method, ed$cell_type),
                   levels = sprintf("%s · %s", ed$method, ed$cell_type))
  p4 <- ggplot(ed, aes(pearson, lab, colour = method)) +
    geom_segment(aes(x = 0, xend = pearson, yend = lab), linewidth = 0.7) +
    geom_point(size = 3.4) +
    scale_color_pub("lancet") +
    scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(colour = NULL, title = "Agreement with the known composition",
         subtitle = "Pearson r between predicted soft score and true cell-type fraction",
         x = "Pearson r", y = NULL) +
    theme_pub(base_size = 11, legend = "bottom")
  save_fig(p4, file.path(opt$assets, "fig4_method_agreement_lollipop"), width = 6.2, height = 4.4)
}

cat(sprintf("完成 · 结果 -> %s · 图 -> %s\n", opt$outdir, opt$assets))
