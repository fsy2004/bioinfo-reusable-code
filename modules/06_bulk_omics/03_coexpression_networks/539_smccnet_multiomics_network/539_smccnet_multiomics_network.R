# =============================================================================
# 编号   : 539
# 脚本名 : 表型驱动多组学稀疏典则相关网络 (SmCCNet · SmCCA + LASSO)
# 分类   : 11_wgcna
# 用途   : 用 SmCCNet 在「两组学(mRNA + miRNA)+ 连续表型」上重建**性状特异**的
#          跨组学子网络:稀疏多重典则相关分析(SmCCA)对 omics-omics 与 omics-trait
#          的成对典则相关加权求和 + LASSO 稀疏化 + 子采样稳健化 → 相似度矩阵 Abar →
#          层次聚类取模块 → trait-driven 跨组学网络(节点=特征,边=共选相似度)。
# ★诚实基线 : 内置「无监督相关网络」对照(getRobustWeightsMulti(NoTrait=TRUE),
#            即不看表型、纯 omics-omics 典则相关)。合成数据故意埋了两块结构:
#            ① trait 块(真表型驱动的跨组学特征)② confounder 块(更强但与表型无关的
#            批次/混杂轴)。无监督网络会被强的 confounder 块「带偏」,而表型驱动子网
#            把 hub 聚焦到 trait 块 —— 报告两者 top-hub 对 gold/confounder 的命中数,
#            用数字证明「表型驱动 = 更聚焦」,而非只画好看的图。
# 依赖   : SmCCNet(2.0.7) · igraph · ggraph · ggplot2 · ggrepel · reshape2
# 运行   : Rscript 539_smccnet_multiomics_network.R                       # 合成示例,CPU 秒级
#          Rscript 539_smccnet_multiomics_network.R --mrna m.csv --mirna mi.csv --pheno y.csv
# 输入   : mrna/mirna = 行=样本、列=特征 的组学矩阵 CSV(首列样本ID);
#          pheno = 行=样本、含 1 列连续表型 的 CSV(首列样本ID)。三表样本须对齐。
# 备注   : 合成数据 synthetic, for demo only。表型驱动用 CCcoef 上调 omics-trait 典则
#          相关权重(对应 fastAutoSmCCNet 的 BetweenShrinkage 思想)。
# =============================================================================

## ---- 定位共享顶刊主题库 _framework/theme_pub.R ------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(SmCCNet); library(igraph); library(ggraph); library(ggplot2)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  mrna   = file.path(DDAT, "mrna.csv"),
  mirna  = file.path(DDAT, "mirna.csv"),
  pheno  = file.path(DDAT, "pheno.csv"),
  outdir = file.path(SCRIPT_DIR, "results"),
  lambda      = 0.5,   # 每组学 LASSO 惩罚 (0-1)
  subsamp     = 200,   # 子采样次数 (越大越稳,越慢)
  s_frac      = 0.7,   # 每次子采样的特征比例
  trait_weight = 5,    # omics-trait 典则相关的相对权重 (>1 = 更聚焦表型)
  top_edges   = 40))   # 网络图保留的最强边数
for (k in c("lambda","subsamp","s_frac","trait_weight","top_edges"))
  args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 0. 合成示例数据(synthetic, for demo only) ----------------------------
# 两组学 mRNA(30) + miRNA(24) + 连续表型 Trait。埋两块跨组学结构:
#   trait 块  : mRNA_1..4 + miR_1..3  受真表型 Y 驱动 (gold standard)
#   conf  块  : mRNA_15..20 + miR_15..20 受 trait-无关潜轴 Z 驱动 (更强的混杂/批次)
# 目的:让「无监督网络被强 confounder 带偏」「表型驱动网络聚焦 trait 块」可被量化。
if (!(file.exists(args$mrna) && file.exists(args$mirna) && file.exists(args$pheno))) {
  cat("Step 0: 生成合成两组学 + 连续表型 (synthetic demo only)...\n")
  n <- 150; p1 <- 30; p2 <- 24
  X1 <- matrix(rnorm(n * p1), n, p1); colnames(X1) <- sprintf("mRNA_%02d", seq_len(p1))
  X2 <- matrix(rnorm(n * p2), n, p2); colnames(X2) <- sprintf("miR_%02d",  seq_len(p2))
  Y  <- rnorm(n)            # 真表型
  Z  <- rnorm(n)            # 与表型无关的潜在混杂轴(批次/人口结构)
  tg1 <- 1:4; tg2 <- 1:3; cf1 <- 15:20; cf2 <- 15:20
  for (j in tg1) X1[, j] <- X1[, j] + 0.9 * Y     # trait 块 mRNA
  for (j in tg2) X2[, j] <- X2[, j] - 0.9 * Y     # trait 块 miR(负向)
  for (j in cf1) X1[, j] <- X1[, j] + 1.0 * Z     # confounder 块(更强)
  for (j in cf2) X2[, j] <- X2[, j] + 1.0 * Z
  sid <- sprintf("S%03d", seq_len(n))
  write.csv(data.frame(SampleID = sid, X1, check.names = FALSE), args$mrna,  row.names = FALSE)
  write.csv(data.frame(SampleID = sid, X2, check.names = FALSE), args$mirna, row.names = FALSE)
  write.csv(data.frame(SampleID = sid, Trait = Y),                args$pheno, row.names = FALSE)
}

## ---- 1. 读入并对齐三张表 ---------------------------------------------------
cat("Step 1: 读 mRNA / miRNA / 表型 并对齐样本...\n")
read_omics <- function(p) { d <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- d[[1]]; m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- rn; m }
X1 <- read_omics(args$mrna); X2 <- read_omics(args$mirna)
ph <- read.csv(args$pheno, check.names = FALSE, stringsAsFactors = FALSE)
rownames(ph) <- ph[[1]]
common <- Reduce(intersect, list(rownames(X1), rownames(X2), rownames(ph)))
X1 <- X1[common, , drop = FALSE]; X2 <- X2[common, , drop = FALSE]
Y  <- as.matrix(ph[common, setdiff(colnames(ph), colnames(ph)[1])[1], drop = FALSE])
colnames(Y) <- "Trait"
# 标准化(SmCCA 对尺度敏感)
X1 <- scale(X1); X2 <- scale(X2); Y <- scale(Y)
Xlist <- list(mRNA = X1, miRNA = X2)
feat  <- c(colnames(X1), colnames(X2))
omics_of <- c(rep("mRNA", ncol(X1)), rep("miRNA", ncol(X2))); names(omics_of) <- feat
cat(sprintf("  n=%d 样本 · mRNA %d 特征 · miRNA %d 特征 · 表型 '%s'\n",
            length(common), ncol(X1), ncol(X2), colnames(Y)))

# gold / confounder 标签(仅合成示例可知,用于诚实基线量化;真实数据自动为空)
gold <- intersect(c(sprintf("mRNA_%02d", 1:4), sprintf("miR_%02d", 1:3)), feat)
conf <- intersect(c(sprintf("mRNA_%02d", 15:20), sprintf("miR_%02d", 15:20)), feat)
has_truth <- length(gold) > 0 && length(conf) > 0

## ---- 2. SmCCA 稳健典则权重:表型驱动 vs 无监督基线 -------------------------
# CCcoef 顺序 = combn(T+1,2),T=2 组学+表型 → 列序: (mRNA-miRNA),(mRNA-trait),(miRNA-trait)
cc_trait <- c(1, args$trait_weight, args$trait_weight)   # 上调 omics-trait → 更聚焦表型
run_smcca <- function(no_trait, cccoef) {
  Ws <- getRobustWeightsMulti(Xlist, Trait = Y, Lambda = c(args$lambda, args$lambda),
          s = c(args$s_frac, args$s_frac), NoTrait = no_trait,
          SubsamplingNum = args$subsamp, CCcoef = cccoef)
  getAbar(Ws, FeatureLabel = feat)              # p×p 相似度(共选频率加权)矩阵
}
cat("Step 2: SmCCA 子采样(表型驱动)...\n")
Abar_t <- run_smcca(no_trait = FALSE, cccoef = cc_trait)
cat("Step 2': SmCCA 子采样(★诚实基线:无监督, NoTrait=TRUE)...\n")
Abar_u <- run_smcca(no_trait = TRUE,  cccoef = NULL)

## ---- 3. 诚实基线量化:hub 聚焦度对比 ---------------------------------------
cat("Step 3: 诚实基线对比(表型驱动 vs 无监督)...\n")
deg_t <- rowSums(Abar_t); deg_u <- rowSums(Abar_u)   # 节点连接强度 = 网络 hub 度
K <- min(length(gold) + length(conf), length(feat)); if (K < 6) K <- min(13, length(feat))
top_t <- names(sort(deg_t, decreasing = TRUE))[seq_len(K)]
top_u <- names(sort(deg_u, decreasing = TRUE))[seq_len(K)]
contrast <- data.frame(
  network = c("Trait-driven SmCCA", "Unsupervised (baseline)"),
  gold_in_top = c(sum(top_t %in% gold), sum(top_u %in% gold)),
  conf_in_top = c(sum(top_t %in% conf), sum(top_u %in% conf)),
  gold_total  = length(gold), conf_total = length(conf), topK = K)
if (has_truth) {
  contrast$focus_ratio <- with(contrast, gold_in_top / pmax(conf_in_top, 1))
  write.csv(contrast, file.path(args$outdir, "honest_baseline_contrast.csv"), row.names = FALSE)
  cat(sprintf("  [诚实基线] top-%d hub 中真表型特征/混杂特征命中:\n", K))
  cat(sprintf("    表型驱动 SmCCA : gold %d/%d · confounder %d/%d\n",
              contrast$gold_in_top[1], length(gold), contrast$conf_in_top[1], length(conf)))
  cat(sprintf("    无监督  基线   : gold %d/%d · confounder %d/%d\n",
              contrast$gold_in_top[2], length(gold), contrast$conf_in_top[2], length(conf)))
  cat("    → 表型驱动子网把 hub 聚焦到真表型特征;无监督网络被更强的混杂块带偏。\n")
} else cat("  (真实数据无 gold/conf 标注,跳过量化对比;仅对比两网络拓扑)\n")

## ---- 4. 表型驱动模块 + 节点表型相关 ----------------------------------------
cat("Step 4: 取跨组学模块 + 节点-表型相关...\n")
mods <- getOmicsModules(Abar_t, CutHeight = 1 - 0.1^10, PlotTree = FALSE)
mods <- mods[order(-vapply(mods, length, 1L))]            # 大模块在前
mod_of <- setNames(rep(NA_integer_, length(feat)), feat)
for (i in seq_along(mods)) mod_of[feat[mods[[i]]]] <- i
# 每个特征与表型的 |相关| (用于着色/排序;不改变网络拓扑)
Xall <- cbind(X1, X2)
trait_cor <- apply(Xall, 2, function(v) suppressWarnings(cor(v, Y[, 1])))
node_tab <- data.frame(feature = feat, omics = omics_of[feat], module = mod_of[feat],
  degree_trait = deg_t[feat], degree_unsup = deg_u[feat],
  trait_cor = trait_cor[feat], abs_trait_cor = abs(trait_cor[feat]),
  is_gold = feat %in% gold, is_conf = feat %in% conf, row.names = NULL)
write.csv(node_tab, file.path(args$outdir, "node_table.csv"), row.names = FALSE)

## ---- 5. 构网络(只保留跨组学边,top_edges 最强)----------------------------
build_net <- function(Abar) {
  A <- Abar; diag(A) <- 0
  em <- which(upper.tri(A) & A > 0, arr.ind = TRUE)
  if (!nrow(em)) return(NULL)
  ed <- data.frame(from = feat[em[, 1]], to = feat[em[, 2]], weight = A[em])
  ed$cross <- omics_of[ed$from] != omics_of[ed$to]    # 跨组学边标记
  ed[order(-ed$weight), ]
}
edges_t <- build_net(Abar_t)
edges_cross <- edges_t[edges_t$cross, ]
edges_top <- head(edges_t, args$top_edges)
write.csv(edges_t, file.path(args$outdir, "trait_network_edges.csv"), row.names = FALSE)

## ===========================================================================
## 顶刊级图(全部 lollipop / dumbbell / 网络 / 热图;禁止平凡条形图)
## ===========================================================================
PAL <- pal_pub(name = "npg")
om_cols <- c(mRNA = PAL[1], miRNA = PAL[2])

# ---- 图1:trait-specific 跨组学子网络(ggraph)----------------------------
cat("Step 5: 图1 表型驱动跨组学子网络...\n")
keep_feat <- unique(c(edges_top$from, edges_top$to))
g <- igraph::graph_from_data_frame(edges_top[, c("from", "to", "weight", "cross")],
       vertices = node_tab[node_tab$feature %in% keep_feat,
         c("feature", "omics", "trait_cor", "abs_trait_cor", "degree_trait")], directed = FALSE)
set.seed(42)
p_net <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight, colour = cross)) +
  scale_edge_width(range = c(0.2, 2.2), name = "Co-selection") +
  scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
  scale_edge_colour_manual(values = c(`TRUE` = PAL[8], `FALSE` = "grey75"),
                           labels = c(`TRUE` = "Cross-omics", `FALSE` = "Within-omics"),
                           name = "Edge type") +
  geom_node_point(aes(fill = omics, size = degree_trait), shape = 21, colour = "grey20", stroke = 0.4) +
  scale_fill_manual(values = om_cols, name = "Omics") +
  scale_size_continuous(range = c(2.5, 9), name = "Hub degree") +
  ggrepel::geom_text_repel(aes(x = x, y = y, label = name), size = 2.6,
                           max.overlaps = 20, segment.size = 0.2) +
  labs(title = "Trait-driven cross-omics subnetwork (SmCCA)",
       subtitle = sprintf("Top %d co-selection edges · node size = hub degree", nrow(edges_top))) +
  theme_pub(base_size = 11) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
save_fig(p_net, file.path(ASSETS, "fig1_trait_subnetwork"), width = 8, height = 6.5)

# ---- 图2:SmCCA 模块邻接热图(最大模块的 Abar 子块)-----------------------
cat("Step 5: 图2 模块邻接热图...\n")
big_mod <- mods[[1]]; mf <- feat[big_mod]
# 按组学 + 表型相关排序,便于看跨组学块结构
ord <- order(omics_of[mf], -abs(trait_cor[mf]))
mf <- mf[ord]
sub <- Abar_t[mf, mf]
hm <- reshape2::melt(sub, varnames = c("row", "col"), value.name = "sim")
hm$row <- factor(hm$row, levels = mf); hm$col <- factor(hm$col, levels = rev(mf))
ax_cols <- om_cols[omics_of[mf]]
p_hm <- ggplot(hm, aes(col, row, fill = sim)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "D", name = "Co-selection\nsimilarity") +
  labs(title = "SmCCA module adjacency (largest trait-driven module)",
       subtitle = "Rows/cols ordered by omics layer then |trait correlation|",
       x = NULL, y = NULL) +
  coord_equal() + theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = ax_cols, size = 6),
        axis.text.y = element_text(colour = rev(ax_cols), size = 6),
        axis.line = element_blank(), axis.ticks = element_line(linewidth = 0.2),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.key       = element_rect(fill = "white", colour = NA))
save_fig(p_hm, file.path(ASSETS, "fig2_module_adjacency_heatmap"),
         width = 7.5, height = 7)

# ---- 图3:edge-weight lollipop(最强跨组学边)------------------------------
cat("Step 5: 图3 跨组学边权 lollipop...\n")
el <- head(edges_cross, 20)
el$pair <- paste(el$from, "—", el$to)
el$pair <- factor(el$pair, levels = rev(el$pair))
p_lol <- ggplot(el, aes(weight, pair)) +
  geom_segment(aes(x = 0, xend = weight, y = pair, yend = pair), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(colour = weight), size = 3.6) +
  scale_colour_viridis_c(option = "C", name = "Co-selection") +
  labs(title = "Top cross-omics edges by SmCCA co-selection",
       subtitle = "mRNA–miRNA pairs with strongest trait-driven similarity",
       x = "Co-selection similarity (Abar)", y = NULL) +
  theme_pub(base_size = 10)
save_fig(p_lol, file.path(ASSETS, "fig3_edge_weight_lollipop"), width = 7.5, height = 6)

# ---- 图4:诚实基线 dumbbell(表型驱动 vs 无监督的 hub 度)------------------
cat("Step 5: 图4 诚实基线 hub 度 dumbbell...\n")
sel <- unique(c(head(names(sort(deg_t, decreasing = TRUE)), 12),
                head(names(sort(deg_u, decreasing = TRUE)), 12)))
db <- data.frame(feature = sel,
  trait = deg_t[sel] / max(deg_t), unsup = deg_u[sel] / max(deg_u),
  cls = ifelse(sel %in% gold, "trait-gold", ifelse(sel %in% conf, "confounder", "other")))
db <- db[order(db$trait), ]; db$feature <- factor(db$feature, levels = db$feature)
cls_cols <- c(`trait-gold` = PAL[3], confounder = PAL[8], other = "grey55")
p_db <- ggplot(db) +
  geom_segment(aes(x = unsup, xend = trait, y = feature, yend = feature), colour = "grey75", linewidth = 0.7) +
  geom_point(aes(x = unsup, y = feature, shape = "Unsupervised (baseline)"), colour = "grey45", size = 2.8) +
  geom_point(aes(x = trait, y = feature, colour = cls, shape = "Trait-driven"), size = 3.2) +
  scale_colour_manual(values = cls_cols, name = "Feature class") +
  scale_shape_manual(values = c(`Unsupervised (baseline)` = 1, `Trait-driven` = 16), name = "Network") +
  labs(title = "Honest baseline: hub focusing under trait guidance",
       subtitle = "Trait-driven SmCCA lifts true trait features; baseline favors the confounder block",
       x = "Normalized hub degree", y = NULL) +
  theme_pub(base_size = 10)
save_fig(p_db, file.path(ASSETS, "fig4_honest_baseline_dumbbell"), width = 8, height = 6.5)

## ---- 6. 收尾 ---------------------------------------------------------------
cat("\n=== 完成 ===\n")
cat("结果表:", normalizePath(args$outdir), "\n")
cat("展示图:", normalizePath(ASSETS), "\n")
if (has_truth)
  cat(sprintf("★诚实基线实测: 表型驱动 top-%d hub 命中真表型 %d/%d(混杂 %d);无监督命中真表型 %d/%d(混杂 %d)。\n",
      K, contrast$gold_in_top[1], length(gold), contrast$conf_in_top[1],
      contrast$gold_in_top[2], length(gold), contrast$conf_in_top[2]))
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
