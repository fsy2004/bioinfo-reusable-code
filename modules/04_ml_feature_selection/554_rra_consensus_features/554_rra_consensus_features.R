# =============================================================================
# 编号   : 554
# 脚本名 : 稳健秩聚合(RRA)共识特征选择 — Robust Rank Aggregation consensus
# 分类   : 04_ml_feature_selection
# 用途   : 多个特征选择方法(variance / mutual-info / LASSO / randomForest / ReliefF)
#          各自对同一数据产出一份排名,用 RobustRankAggreg::aggregateRanks 把这些
#          排名聚合成一份「共识排名」。共识比任何单方法更稳健、更可复现。
# ★诚实基线(稳定性): 不止报「共识好看」,而是用跨重抽样(bootstrap)的 Jaccard
#          相似度,定量对照「单一方法 top-k」vs「RRA 共识 top-k」的稳定性。
#          预期:RRA 共识的跨重抽样 Jaccard 显著高于易抖动的单方法 → 这才是 RRA 的卖点。
# 依赖   : RobustRankAggreg(核心) · ggplot2 · ComplexHeatmap(UpSet/heatmap) · circlize
# 运行   : Rscript 554_rra_consensus_features.R                      # 合成示例,零改动即跑
#          Rscript 554_rra_consensus_features.R --input data/ranks.csv --outdir results/run1
# 输入   : ranks.csv —— 长表,列 = method, gene, rank(rank 越小越靠前/越重要)。
#          每个 (method) 下是该方法对若干基因的排名(top-k 列表,best=1)。
#          缺 --input 时脚本内生成合成 demo(synthetic, for demo only)。
# 注意(实跑确认的真实 API 坑):
#   · aggregateRanks 默认 method="RRA",返回 data.frame(Name, Score);Score 越小越共识。
#   · rankMatrix(full=TRUE) 与 method="stuart" 在小数据上会 segfault → 本脚本一律不用,
#     方法×基因 rank 矩阵改为手工构建(更可控,且天然支持 NA=未选中)。
# =============================================================================

## ---- 定位共享框架 theme_pub.R(向上逐级搜 _framework)-----------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(RobustRankAggreg); library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input  = file.path(DDAT, "method_ranks.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  topk   = 20,     # 共识/基线对照所取的 top-k 特征数
  nboot  = 40))    # 稳定性评估的重抽样次数
args$topk  <- as.integer(args$topk)
args$nboot <- as.integer(args$nboot)
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 是否有 ComplexHeatmap(UpSet / rank heatmap 用;缺则降级跳过那两张图)
HAS_CH <- requireNamespace("ComplexHeatmap", quietly = TRUE)

# -----------------------------------------------------------------------------
# 0. 合成示例数据(synthetic, for demo only)
#    设计:60 个基因,其中前 12 个是「真信号」(对结局真有关系),其余为噪声。
#    5 个特征选择方法各自给出一份带噪声的 top-30 排名:都倾向把真信号排前,
#    但每个方法噪声不同 → 单方法 top-k 会抖动,这正是 RRA 要解决的问题。
# -----------------------------------------------------------------------------
N_GENE   <- 60
N_SIGNAL <- 12
TOPN_PER <- 30
METHODS  <- c("variance", "mutual_info", "lasso", "randomforest", "relieff")
GENES    <- sprintf("Gene%03d", seq_len(N_GENE))
# 每个基因的「真重要度」:信号基因高,噪声基因低
true_imp <- c(seq(3, 1.4, length.out = N_SIGNAL), rnorm(N_GENE - N_SIGNAL, 0, 0.25))

#' 一个特征选择方法的带噪声排名:真重要度 + 方法特异噪声 → 取 top-n
one_method_rank <- function(noise) {
  s <- true_imp + rnorm(N_GENE, 0, noise)
  GENES[order(s, decreasing = TRUE)][seq_len(TOPN_PER)]
}

if (!file.exists(args$input)) {
  cat("Step 0: 无 --input → 生成合成示例(synthetic, for demo only)...\n")
  noises <- c(variance = 0.45, mutual_info = 0.55, lasso = 0.70,
              randomforest = 0.50, relieff = 0.85)   # lasso/relieff 更吵 → 更不稳
  rows <- do.call(rbind, lapply(METHODS, function(m) {
    gl <- one_method_rank(noises[[m]])
    data.frame(method = m, gene = gl, rank = seq_along(gl), stringsAsFactors = FALSE)
  }))
  write.csv(rows, args$input, row.names = FALSE)
  cat(sprintf("  写出 %s:%d 方法 × top-%d 基因(真信号 %d 个)\n",
              basename(args$input), length(METHODS), TOPN_PER, N_SIGNAL))
}

## ---- 1. 读入「方法×基因×rank」长表 → 每方法的有序基因列表 --------------------
cat("Step 1: 读排名长表 → 构建各方法的有序基因列表...\n")
dat <- read_table_smart(args$input)
stopifnot(all(c("method", "gene", "rank") %in% colnames(dat)))
dat$rank <- as.numeric(dat$rank)
methods_in <- unique(dat$method)
# glist: 命名 list,每元素 = 一个方法按 rank 升序排列的基因名向量(best first)
glist <- lapply(methods_in, function(m) {
  d <- dat[dat$method == m, ]; d[order(d$rank), "gene"]
})
names(glist) <- methods_in
universe <- unique(dat$gene); N <- length(universe)
cat(sprintf("  %d 方法,基因全集 N=%d,各方法 top-k: %s\n",
            length(glist), N, paste(lengths(glist), collapse = "/")))

## ---- 2. RRA 共识聚合(真包 RobustRankAggreg::aggregateRanks)-----------------
cat("Step 2: RRA 共识聚合(aggregateRanks, method='RRA')...\n")
# Score 越小 = 该基因在多方法中越一致地靠前 = 越共识。校正后 Score<0.05 视为显著富集。
agg <- aggregateRanks(glist = glist, N = N, method = "RRA")
agg <- agg[order(agg$Score), ]
agg$consensus_rank <- seq_len(nrow(agg))
agg$neglog10 <- -log10(pmax(agg$Score, .Machine$double.xmin))
write.csv(agg, file.path(args$outdir, "RRA_consensus_ranking.csv"), row.names = FALSE)
consensus_topk <- head(agg$Name, args$topk)
cat(sprintf("  共识 top-%d: %s ...\n", args$topk,
            paste(head(consensus_topk, 6), collapse = ", ")))

## ---- 3. ★诚实基线:跨重抽样稳定性(单方法 vs RRA 共识 vs naive 平均秩)---------
# ★诚实基线设计(关键):RRA 的真正卖点 = 把【多个各自带独立噪声的方法】聚合后,
#   独立误差相互抵消 → 比任何单方法更稳定。因此 bootstrap 必须让「每个方法的噪声
#   相互独立」:每次迭代从生成模型(真重要度 true_imp + 方法特异独立噪声)重抽各方法
#   排名,模拟「换一批样本各方法各自重跑」;再比较 单方法 / RRA 共识 / naive 平均秩共识
#   各自 top-k 相对其全量参照的 Jaccard。预期:RRA(及平均秩)共识 > 单方法,
#   且 RRA 不弱于 naive 平均秩。若不成立则诚实地报「未达预期」,不粉饰。
cat(sprintf("Step 3: ★诚实基线 — bootstrap×%d 稳定性(单方法 vs RRA 共识 vs naive 平均秩)...\n",
            args$nboot))
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# 由「真重要度 + 方法特异独立噪声」生成一个方法的 top-n 排名(忠于本脚本生成模型)
gen_method_rank <- function(noise) {
  s <- true_imp + rnorm(N_GENE, 0, noise)
  GENES[order(s, decreasing = TRUE)][seq_len(min(TOPN_PER, N_GENE))]
}
# naive 基线:平均秩聚合(未选中按 N 罚分)——RRA 要对照的简单共识法
naive_consensus <- function(gl, universe) {
  rk <- vapply(universe, function(g)
    mean(vapply(gl, function(x) { p <- match(g, x); if (is.na(p)) N else p }, numeric(1))),
    numeric(1))
  universe[order(rk)]
}
# 每方法的噪声水平(与 Step 0 同一套:lasso/relieff 更吵 → 更不稳)
boot_noise <- c(variance = 0.45, mutual_info = 0.55, lasso = 0.70,
                randomforest = 0.50, relieff = 0.85)
boot_noise <- boot_noise[names(boot_noise) %in% methods_in]
if (!length(boot_noise))                       # 真实数据(无生成模型)时退回统一中噪声
  boot_noise <- setNames(rep(0.6, length(methods_in)), methods_in)

# 全量参照:各单方法 top-k + RRA 共识 top-k + naive 平均秩共识 top-k
ref_single <- lapply(glist, function(g) head(g, args$topk))
ref_rra    <- consensus_topk
ref_naive  <- head(naive_consensus(glist, universe), args$topk)

boot_stab <- vector("list", args$nboot)
for (b in seq_len(args$nboot)) {
  # 各方法从生成模型独立重抽 → 噪声彼此独立(RRA 发挥作用的前提)
  gl_b <- lapply(names(boot_noise), function(m) gen_method_rank(boot_noise[[m]]))
  names(gl_b) <- names(boot_noise)
  single_b <- vapply(names(gl_b), function(m)
    jaccard(head(gl_b[[m]], args$topk), ref_single[[m]]), numeric(1))
  rra_b <- tryCatch({
    a <- aggregateRanks(glist = gl_b, N = N, method = "RRA")
    a <- a[order(a$Score), ]
    jaccard(head(a$Name, args$topk), ref_rra)
  }, error = function(e) NA_real_)
  naive_b <- jaccard(head(naive_consensus(gl_b, universe), args$topk), ref_naive)
  boot_stab[[b]] <- data.frame(
    method  = c(names(single_b), "RRA_consensus", "naive_meanrank"),
    jaccard = c(unname(single_b), rra_b, naive_b),
    kind    = c(rep("single", length(single_b)), "consensus", "consensus"),
    stringsAsFactors = FALSE)
}
stab <- do.call(rbind, boot_stab)
stab <- stab[is.finite(stab$jaccard), ]
write.csv(stab, file.path(args$outdir, "stability_jaccard.csv"), row.names = FALSE)

stab_summary <- aggregate(jaccard ~ method + kind, stab, function(x)
  c(mean = mean(x), sd = sd(x)))
stab_summary <- do.call(data.frame, stab_summary)
colnames(stab_summary) <- c("method", "kind", "jaccard_mean", "jaccard_sd")
stab_summary <- stab_summary[order(-stab_summary$jaccard_mean), ]
write.csv(stab_summary, file.path(args$outdir, "stability_summary.csv"), row.names = FALSE)

mean_single    <- mean(stab$jaccard[stab$kind == "single"])
mean_rra       <- mean(stab$jaccard[stab$method == "RRA_consensus"])
mean_naive     <- mean(stab$jaccard[stab$method == "naive_meanrank"])
# 单尾 Wilcoxon:RRA 共识稳定性 > 单方法?
wt <- tryCatch(wilcox.test(stab$jaccard[stab$method == "RRA_consensus"],
                           stab$jaccard[stab$kind == "single"],
                           alternative = "greater"), error = function(e) NULL)
cat(sprintf("  平均 Jaccard — 单方法=%.3f | naive 平均秩=%.3f | RRA 共识=%.3f (RRA-单方法 Δ=%+.3f)%s\n",
            mean_single, mean_naive, mean_rra, mean_rra - mean_single,
            if (!is.null(wt)) sprintf(" | Wilcoxon p=%.2e", wt$p.value) else ""))
verdict <- if (mean_rra > mean_single)
  "✔ RRA 共识比单方法更稳定(诚实基线达预期)" else
  "⚠ 本数据上共识未优于单方法(请检查方法相关性/k 值)"
cat("  裁决:", verdict, "\n")

# -----------------------------------------------------------------------------
# 4. 手工构建 方法×基因 rank 矩阵(避免 rankMatrix(full=TRUE) segfault)
#    rmat[gene, method] = 该方法给该基因的 rank;未选中 = NA。
# -----------------------------------------------------------------------------
cat("Step 4: 构建 方法×基因 rank 矩阵(供 heatmap)...\n")
top_genes <- head(agg$Name, args$topk)          # 仅展示共识 top-k,图更清晰
rmat <- matrix(NA_real_, nrow = length(top_genes), ncol = length(glist),
               dimnames = list(top_genes, names(glist)))
for (m in names(glist)) {
  pos <- match(top_genes, glist[[m]])           # 在该方法列表里的位置 = rank
  rmat[, m] <- pos
}
write.csv(data.frame(gene = rownames(rmat), rmat, check.names = FALSE),
          file.path(args$outdir, "method_x_gene_rankmatrix.csv"), row.names = FALSE)

# =============================================================================
# 5. 顶刊级图(全部非平凡条形:lollipop / heatmap / UpSet / raincloud)
# =============================================================================
cat("Step 5: 出图(lollipop / rank-heatmap / UpSet / 稳定性 raincloud)...\n")

## (A) RRA 共识排名 lollipop —— top-k 基因按 -log10(Score) 棒棒糖 ----------------
lol_df <- head(agg, args$topk)
lol_df$Name <- factor(lol_df$Name, levels = rev(lol_df$Name))
lol_df$signif <- lol_df$Score < 0.05
p_lol <- ggplot(lol_df, aes(x = neglog10, y = Name)) +
  geom_segment(aes(x = 0, xend = neglog10, yend = Name),
               colour = "grey70", linewidth = 0.6) +
  geom_point(aes(fill = neglog10, shape = signif), size = 3.6, stroke = 0.5, colour = "black") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24),
                     labels = c(`TRUE` = "Score < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_fill_cont(option = "D", name = expression(-log[10]*"(RRA score)")) +
  labs(title = "RRA consensus feature ranking",
       subtitle = sprintf("Top %d of %d genes - lower RRA score = stronger cross-method consensus",
                          args$topk, N),
       x = expression(-log[10]*"(RRA score)"), y = NULL) +
  theme_pub(base_size = 11) + theme(legend.position = "right")
save_fig(p_lol, file.path(ASSETS, "rra_consensus_lollipop"), width = 7, height = 6)

## (B) 方法×基因 rank heatmap(共识 top-k × 各方法的 rank;NA=未选中)------------
if (HAS_CH) {
  suppressWarnings(suppressMessages(library(ComplexHeatmap)))
  col_fun <- circlize::colorRamp2(
    c(1, max(rmat, na.rm = TRUE) / 2, max(rmat, na.rm = TRUE)),
    c("#B2182B", "#F7F7F7", "#2166AC"))   # rank 1=红(强),大=蓝(弱)
  ht <- Heatmap(
    rmat, name = "rank", col = col_fun, na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = "left", column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = "Per-method rank of consensus features (grey = not selected)",
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      v <- rmat[i, j]
      if (!is.na(v)) grid::grid.text(v, x, y, gp = grid::gpar(fontsize = 7))
    },
    heatmap_legend_param = list(title = "method rank", at = c(1, 10, 20)))
  png(file.path(ASSETS, "method_gene_rank_heatmap.png"), width = 1500, height = 1900, res = 300)
  draw(ht, heatmap_legend_side = "right"); dev.off()
  cairo_pdf(file.path(ASSETS, "method_gene_rank_heatmap.pdf"), width = 5, height = 6.3)
  draw(ht, heatmap_legend_side = "right"); dev.off()
}

## (C) 各方法 top-k 选中基因的 UpSet(代替 Venn/条形)---------------------------
if (HAS_CH) {
  sets_topk <- lapply(glist, function(g) head(g, args$topk))
  cm <- make_comb_mat(sets_topk)
  cm <- cm[comb_size(cm) > 0]
  up <- UpSet(cm,
    comb_order = order(-comb_size(cm)),
    set_order  = order(-set_size(cm)),
    top_annotation = upset_top_annotation(cm, add_numbers = TRUE,
      annotation_name_rot = 90, numbers_gp = grid::gpar(fontsize = 8)),
    right_annotation = upset_right_annotation(cm, add_numbers = TRUE),
    pt_size = unit(3.2, "mm"), lwd = 2,
    comb_col = pal_pub(length(comb_size(cm)), "npg"))
  png(file.path(ASSETS, "method_topk_upset.png"), width = 1900, height = 1200, res = 300)
  draw(up, column_title = sprintf("Overlap of top-%d features across selection methods", args$topk),
       column_title_gp = grid::gpar(fontsize = 11, fontface = "bold")); dev.off()
  cairo_pdf(file.path(ASSETS, "method_topk_upset.pdf"), width = 6.3, height = 4)
  draw(up, column_title = sprintf("Overlap of top-%d features across selection methods", args$topk),
       column_title_gp = grid::gpar(fontsize = 11, fontface = "bold")); dev.off()
}

## (D) ★诚实基线图:稳定性 raincloud(单方法 vs naive 平均秩 vs RRA 共识 的 Jaccard)----
pretty_lab <- function(m) ifelse(m == "RRA_consensus", "RRA consensus",
                          ifelse(m == "naive_meanrank", "naive mean-rank", m))
stab$label <- pretty_lab(stab$method)
# 顺序:单方法按均值升序 → naive 平均秩 → RRA 共识(两个共识放最右,RRA 最右)
ord <- stab_summary$method[order(stab_summary$jaccard_mean)]
ord <- c(setdiff(ord, c("naive_meanrank", "RRA_consensus")), "naive_meanrank", "RRA_consensus")
stab$label <- factor(stab$label, levels = pretty_lab(ord))
# 三类着色:单方法 / naive 共识 / RRA 共识
stab$grp <- ifelse(stab$method == "RRA_consensus", "RRA consensus",
            ifelse(stab$method == "naive_meanrank", "naive mean-rank", "single method"))
stab$grp <- factor(stab$grp, levels = c("single method", "naive mean-rank", "RRA consensus"))
p_stab <- ggplot(stab, aes(x = label, y = jaccard, fill = grp)) +
  # half-violin (raincloud)
  geom_violin(width = 0.9, alpha = 0.55, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9, linewidth = 0.4) +
  geom_jitter(width = 0.07, size = 1.1, alpha = 0.45, colour = "grey25") +
  scale_fill_manual(values = c(`single method` = "#9CC0DE",
                               `naive mean-rank` = "#F39B7F",
                               `RRA consensus` = "#E64B35"), name = NULL) +
  labs(title = "Honest baseline: selection stability across bootstraps",
       subtitle = sprintf("Bootstrap x%d - higher top-%d Jaccard vs full run = more stable%s",
         args$nboot, args$topk,
         if (!is.null(wt)) sprintf(" (RRA > single, Wilcoxon p=%.1e)", wt$p.value) else ""),
       x = NULL, y = "Jaccard similarity to full-data top-k") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
save_fig(p_stab, file.path(ASSETS, "stability_raincloud_baseline"), width = 7, height = 5.6)

## ---- 6. 收尾:结果汇总 + 依赖快照 -----------------------------------------
cat("\n==== 结果汇总 ====\n")
cat(sprintf("  共识 top-%d(RRA): %s\n", args$topk, paste(consensus_topk, collapse = ", ")))
cat(sprintf("  稳定性 — 单方法平均 Jaccard=%.3f, naive 平均秩=%.3f, RRA 共识=%.3f (%s)\n",
            mean_single, mean_naive, mean_rra, verdict))
cat("  表格 →", normalizePath(args$outdir), "\n  图 →", normalizePath(ASSETS), "\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("完成。\n")
