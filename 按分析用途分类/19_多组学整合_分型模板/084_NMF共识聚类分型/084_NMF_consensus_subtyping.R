# =============================================================================
# 编号       : 084
# 脚本名     : NMF + ConsensusClusterPlus 分子分型 (turnkey + 顶刊图)
# 分类       : 19_多组学整合_分型模板
# 用途       : 对特征×样本矩阵(表达/免疫评分/通路评分)做无监督分子分型,输出
#              NMF 秩选择曲线、共识矩阵热图、分型注释特征热图。
# 方法/包    : NMF + ConsensusClusterPlus + ComplexHeatmap;主题 theme_pub.R
# 结果图     : NMF_rank_survey;Consensus_matrix(k=best);Subtype_heatmap
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 084_NMF_consensus_subtyping.R
# 运行(自己): Rscript 084_NMF_consensus_subtyping.R --input data/feature_matrix.csv --kmax 6
# 可选参数 : --kmin 2 --kmax 6 --nrun 10
# 输入规格 : CSV,首列特征名,其余列=样本(值非负;表达/评分矩阵)。
# 整理日期 : 2026-06-23(turnkey 重构;补 theme_pub 顶刊图)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(ggplot2); library(ComplexHeatmap); library(circlize) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "feature_matrix.csv"),
                      outdir = file.path(SCRIPT_DIR, "results"), kmin = "2", kmax = "6", nrun = "10", k = "0"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)
ranks <- as.integer(args$kmin):as.integer(args$kmax)

cat("Step 1/4: 读取矩阵...\n")
x <- as.matrix(read_table_smart(args$input, row_names = TRUE)); storage.mode(x) <- "double"
x <- x[apply(x, 1, var) > 0, , drop = FALSE]; x <- x - min(x, na.rm = TRUE)
cat("  ", nrow(x), "特征 x", ncol(x), "样本 · 测试 k =", paste(range(ranks), collapse = "-"), "\n")

# ---- Step 2. NMF 秩选择 ----
cat("Step 2/4: NMF 秩选择(可能稍慢)...\n")
suppressWarnings(suppressMessages(library(NMF)))
surv <- nmf(x, rank = ranks, nrun = as.integer(args$nrun), seed = 1, .options = "v0")
coph <- summary(surv)$cophenetic
best_k <- if (as.integer(args$k) > 0) as.integer(args$k) else ranks[which.max(coph)]  # --k 手动指定;否则 cophenetic 最高
cat("  最佳秩 k =", best_k, "(cophenetic", sprintf("%.3f", max(coph)), ")\n")
sdf <- data.frame(k = ranks, cophenetic = coph)
p_rank <- ggplot(sdf, aes(k, cophenetic)) + geom_line(colour = pal_pub(1, "npg"), linewidth = .8) +
  geom_point(size = 2.5, colour = pal_pub(1, "npg")) +
  geom_point(data = sdf[sdf$k == best_k, ], size = 4, shape = 21, fill = "#E64B35", colour = "black") +
  scale_x_continuous(breaks = ranks) +
  labs(title = "NMF rank selection", x = "Rank (k)", y = "Cophenetic correlation") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_rank, file.path(ASSETS, "NMF_rank_survey"), 6, 5); save_fig(p_rank, file.path(args$outdir, "NMF_rank_survey"), 6, 5)
best <- nmf(x, rank = best_k, nrun = as.integer(args$nrun) * 2, seed = 1)
sub <- predict(best)
write.csv(data.frame(sample = colnames(x), subtype = paste0("C", as.integer(sub))), file.path(args$outdir, "sample_subtypes.csv"), row.names = FALSE)

# ---- Step 3. 共识矩阵热图(ConsensusClusterPlus)----
cat("Step 3/4: 共识聚类...\n")
ccp <- ConsensusClusterPlus::ConsensusClusterPlus(x, maxK = max(ranks), reps = 100, pItem = .8, pFeature = .8,
        clusterAlg = "hc", distance = "pearson", seed = 1, plot = NULL, title = tempfile())
cons <- ccp[[best_k]]$consensusMatrix; ord <- ccp[[best_k]]$consensusTree$order
cm <- cons[ord, ord]; cl <- ccp[[best_k]]$consensusClass[ord]
ha <- HeatmapAnnotation(Cluster = paste0("C", cl), col = list(Cluster = setNames(pal_pub(best_k, "npg"), paste0("C", 1:best_k))))
ht <- Heatmap(cm, name = "consensus", col = colorRamp2(c(0, 1), c("white", "#3C5488")),
              cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
              top_annotation = ha, column_title = paste0("Consensus matrix (k=", best_k, ")"),
              column_title_gp = gpar(fontsize = 12, fontface = "bold"))
for (dest in c(file.path(ASSETS, "Consensus_matrix"), file.path(args$outdir, "Consensus_matrix"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 6.5, height = 6); draw(ht); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 6.5, height = 6, units = "in", res = 300); draw(ht); dev.off()
}

# ---- Step 4. 分型注释特征热图 ----
cat("Step 4/4: 分型特征热图...\n")
so <- order(as.integer(sub)); xm <- t(scale(t(x[, so])))
ha2 <- HeatmapAnnotation(Subtype = paste0("C", as.integer(sub)[so]), col = list(Subtype = setNames(pal_pub(best_k, "npg"), paste0("C", 1:best_k))))
ht2 <- Heatmap(xm, name = "Z-score", col = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")),
               top_annotation = ha2, cluster_columns = FALSE, show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 7), column_title = "Feature expression by subtype")
for (dest in c(file.path(ASSETS, "Subtype_heatmap"), file.path(args$outdir, "Subtype_heatmap"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 8, height = 6); draw(ht2); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 8, height = 6, units = "in", res = 300); draw(ht2); dev.off()
}
cat("完成。分型图/表见", normalizePath(args$outdir), "\n")
