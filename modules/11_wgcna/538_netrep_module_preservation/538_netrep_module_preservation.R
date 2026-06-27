# =============================================================================
# 编号   : 538
# 脚本名 : NetRep 跨数据集模块保守性 / 可重复性置换检验 (turnkey + 顶刊图)
# 分类   : 11_wgcna
# 用途   : 给定 discovery + test 两套表达数据 + 一套模块标签(来自 054 WGCNA),
#          用 NetRep::modulePreservation 做 **无分布假设的置换检验**,评估每个共表达
#          模块在独立队列中是否保守 / 可重复。输出 7 个保守性统计的无偏置换 p 值,
#          并据置换 null 计算 WGCNA 式 Zsummary(Z=2/10 阈值线)。
#          比 WGCNA::modulePreservation 快很多,且不假设正态。
# ★诚实基线: NetRep 的【置换 null 分布】本身即基线 —— 每个统计的观测值都与「把
#            test 队列节点标签随机打乱(null='overlap')」得到的 null 分布对照,
#            得到无偏 p 与 Z(observed vs permuted)。模块若只是「碰巧像」,其观测
#            统计会落在 null 分布内 (p≫0.05, Z<2)。脚本合成数据中【刻意混入不保守
#            模块】,基线必须把它们判为不保守,才证明管道有判别力(非只报好看指标)。
#            外部真实队列验证则配对库内 054 WGCNA 模块产出。
# 依赖   : NetRep (核心) · ggplot2 · (theme_pub.R 顶刊主题)
# 运行   : Rscript 538_netrep_module_preservation.R                       # 合成示例
#          Rscript 538_netrep_module_preservation.R \
#              --disc_expr disc.csv --test_expr test.csv --modules mod.csv --nperm 10000
# 输入   : disc_expr / test_expr = 表达矩阵 CSV(首列基因名,其余列=样本;两队列基因须可对齐)
#          modules = 两列 CSV(gene,module);module 为整数标签(0=背景/未分配,不检验)
# 整理日期: 2026-06-27(turnkey 新建;接地 NetRep 1.2.x 真实 API,实跑验证)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(NetRep); library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  disc_expr = file.path(DDAT, "discovery_expr.csv"),
  test_expr = file.path(DDAT, "test_expr.csv"),
  modules   = file.path(DDAT, "module_labels.csv"),
  outdir    = file.path(SCRIPT_DIR, "results"),
  nperm     = 10000,         # 置换次数;p 精度上限 ≈ 1/(nperm+1)
  nthreads  = 1))            # 单核即可,合成规模秒级
args$nperm    <- as.integer(args$nperm)
args$nthreads <- as.integer(args$nthreads)
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Step 0. 合成示例数据 (synthetic, for demo only) — 内生成两队列 + 模块标签
#   设计:5 个模块,每模块 30 基因 + 60 背景基因 = 210 基因。
#     · 模块 M1,M2,M3 = **保守**:discovery 与 test 共享同一潜在因子 → 两队列里
#       模块内基因强相关 → 保守性统计应显著 (p<0.05, Zsummary>10)。
#     · 模块 M4,M5 = **不保守**:只有 discovery 队列里相关(共享潜因子),test 队列
#       里改用各基因独立噪声(无共享因子)→ 模块结构在 test 中瓦解 → 诚实基线必须
#       把它们判为 not preserved (p≫0.05, Zsummary<2)。这是管道判别力的阳性对照。
# =============================================================================
# loadings = 每模块的基因载荷向量(list);对【保守】模块两队列须共用同一 loadings,
# 这样不仅模块密度(avg.weight/coherence/avg.cor)保守,连「节点连接度的相对排序」
# (cor.cor/cor.degree/cor.contrib 这三个 topology 统计)也保守 → 7 统计齐显著、Z 高 + p 小,
# 才是干净的阳性对照(否则 loadings 各队列独立则 topology 不保守、max-p 失真)。
make_cohort <- function(n_samp, mod_sizes, n_bg, preserved, loadings, seed) {
  set.seed(seed)
  blocks <- list()
  for (k in seq_along(mod_sizes)) {
    sz <- mod_sizes[k]
    if (preserved[k]) {
      # 共享潜在因子 f(队列内自有,样本不同) + 共享 loadings(跨队列同一向量)→
      # 队列内强共表达 + 跨队列同一连接拓扑 → 全 7 统计保守
      f   <- rnorm(n_samp)
      mat <- outer(f, loadings[[k]]) + matrix(rnorm(n_samp * sz, 0, 0.5), n_samp, sz)
    } else {
      # 无共享因子 → 纯独立噪声 → 模块结构在该队列瓦解
      mat <- matrix(rnorm(n_samp * sz, 0, 1), n_samp, sz)
    }
    blocks[[k]] <- mat
  }
  bg <- matrix(rnorm(n_samp * n_bg, 0, 1), n_samp, n_bg)  # 背景基因:始终独立噪声
  do.call(cbind, c(blocks, list(bg)))
}

if (!(file.exists(args$disc_expr) && file.exists(args$test_expr) && file.exists(args$modules))) {
  cat("Step 0: 生成合成示例数据 (synthetic demo only)...\n")
  mod_sizes <- c(20L, 30L, 45L, 25L, 40L); n_bg <- 60L   # 不同模块大小,便于散点图区分
  n_genes <- sum(mod_sizes) + n_bg
  gene_id <- sprintf("Gene_%03d", seq_len(n_genes))
  labels  <- c(rep(seq_along(mod_sizes), times = mod_sizes), rep(0L, n_bg))   # 0 = 背景
  # 模块 1,2,3 保守;模块 4,5 不保守(test 队列里结构瓦解)
  pres_disc <- c(TRUE, TRUE, TRUE, TRUE, TRUE)   # discovery 里 5 个模块都有结构
  pres_test <- c(TRUE, TRUE, TRUE, FALSE, FALSE) # test 里仅前 3 个保留结构
  # 各模块 loadings 一次性生成,两队列共用 → 保守模块连拓扑都一致(干净阳性对照)
  set.seed(7); loadings <- lapply(mod_sizes, function(sz) runif(sz, 0.4, 1.4))
  disc <- make_cohort(120, mod_sizes, n_bg, pres_disc, loadings, seed = 42)
  test <- make_cohort(110, mod_sizes, n_bg, pres_test, loadings, seed = 99)
  colnames(disc) <- colnames(test) <- gene_id
  rownames(disc) <- sprintf("Disc_%03d", seq_len(nrow(disc)))
  rownames(test) <- sprintf("Test_%03d", seq_len(nrow(test)))
  # 写出:表达矩阵 = 基因 x 样本(首列基因名);标签 = gene,module 两列
  write.csv(cbind(gene = gene_id, as.data.frame(t(disc))), args$disc_expr, row.names = FALSE)
  write.csv(cbind(gene = gene_id, as.data.frame(t(test))), args$test_expr, row.names = FALSE)
  write.csv(data.frame(gene = gene_id, module = labels), args$modules, row.names = FALSE)
  cat(sprintf("  合成: %d 基因 x (disc %d / test %d 样本); 模块大小 %s; 模块 1-3=保守, 4-5=不保守(阳性对照)\n",
              n_genes, nrow(disc), nrow(test), paste(mod_sizes, collapse = "/")))
}

# =============================================================================
# Step 1. 读入两队列表达 + 模块标签,对齐基因
# =============================================================================
cat("Step 1: 读取 discovery / test 表达 + 模块标签...\n")
read_expr <- function(path) {
  df <- read_table_smart(path, row_names = TRUE)   # 行=基因,列=样本
  as.matrix(df)
}
disc_g <- read_expr(args$disc_expr)
test_g <- read_expr(args$test_expr)
mod_df <- read_table_smart(args$modules, row_names = FALSE)
stopifnot(all(c("gene", "module") %in% colnames(mod_df)))
module_labels <- setNames(as.integer(mod_df$module), mod_df$gene)

common <- Reduce(intersect, list(rownames(disc_g), rownames(test_g), names(module_labels)))
disc_g <- disc_g[common, , drop = FALSE]
test_g <- test_g[common, , drop = FALSE]
module_labels <- module_labels[common]
cat(sprintf("  对齐后基因=%d; 模块(非0): %s\n",
            length(common), paste(sort(setdiff(unique(module_labels), 0)), collapse = ",")))

# NetRep 要求 data = 样本 x 基因;correlation/network = 基因 x 基因
disc_data <- t(disc_g); test_data <- t(test_g)
cat("Step 2: 构建相关 + 邻接(network)矩阵 (Pearson, |r| 作权重)...\n")
disc_corr <- cor(disc_data); test_corr <- cor(test_data)
disc_net  <- abs(disc_corr); test_net  <- abs(test_corr)   # 无符号邻接 = |相关|

# =============================================================================
# Step 3. ★核心:NetRep 置换检验 (诚实基线 = 观测 vs 置换 null)
#   network/data/correlation 三套均为 named list(每队列一项)。
#   null='overlap' = 在两队列共有节点中随机打乱 test 标签生成 null。
# =============================================================================
cat(sprintf("Step 3: NetRep::modulePreservation  (nPerm=%d, null=overlap)...\n", args$nperm))
data_list <- list(discovery = disc_data, test = test_data)
corr_list <- list(discovery = disc_corr, test = test_corr)
net_list  <- list(discovery = disc_net,  test = test_net)

pres <- modulePreservation(
  network = net_list, data = data_list, correlation = corr_list,
  moduleAssignments = module_labels,
  discovery = "discovery", test = "test",
  nPerm = args$nperm, null = "overlap", alternative = "greater",
  nThreads = args$nthreads, simplify = TRUE, verbose = FALSE)

stat_names <- colnames(pres$observed)
mods <- rownames(pres$observed)
cat(sprintf("  完成: %d 模块 x %d 保守性统计\n", length(mods), length(stat_names)))

# ---- 置换 Z-score:每统计 (obs - mean_null)/sd_null;Zsummary = 7 统计均值 ----
# 这是 WGCNA Zsummary 的 NetRep 等价:基于真实置换 null,不假设正态生成 null。
nu <- pres$nulls   # 维度 = 模块 x 统计 x nPerm
zmat <- matrix(NA_real_, length(mods), length(stat_names),
               dimnames = list(mods, stat_names))
for (m in seq_along(mods)) for (s in seq_along(stat_names)) {
  v <- nu[m, s, ]; sdv <- sd(v)
  zmat[m, s] <- if (is.finite(sdv) && sdv > 0) (pres$observed[m, s] - mean(v)) / sdv else NA_real_
}
zsummary <- rowMeans(zmat, na.rm = TRUE)
mod_size <- pres$nVarsPresent[mods]
# 模块级 p:7 统计取最大 p(最保守判定,WGCNA modulePreservation 惯例)
max_p <- apply(pres$p.values, 1, max)

summary_df <- data.frame(
  module   = mods,
  size     = as.integer(mod_size[mods]),
  Zsummary = round(zsummary[mods], 3),
  maxP     = signif(max_p[mods], 3),
  preserved = ifelse(zsummary[mods] >= 10, "strong (Z>=10)",
               ifelse(zsummary[mods] >= 2, "moderate (2<=Z<10)", "not preserved (Z<2)")),
  row.names = NULL)
write.csv(summary_df, file.path(args$outdir, "preservation_summary.csv"), row.names = FALSE)
write.csv(cbind(module = mods, round(pres$observed, 4)), file.path(args$outdir, "observed_statistics.csv"), row.names = FALSE)
write.csv(cbind(module = mods, signif(pres$p.values, 4)), file.path(args$outdir, "permutation_pvalues.csv"), row.names = FALSE)
write.csv(cbind(module = mods, round(zmat, 3)),           file.path(args$outdir, "permutation_zscores.csv"), row.names = FALSE)

cat("  --- 模块保守性小结 (诚实基线判定) ---\n")
print(summary_df)
cat(sprintf("  ★基线判别力: 模块 %s 判为 not-preserved (阳性对照应命中合成的不保守模块 4,5)\n",
            paste(summary_df$module[summary_df$Zsummary < 2], collapse = ",")))

# =============================================================================
# Step 4. 顶刊级图(禁平凡条形图)
# =============================================================================
cat("Step 4: 出图 (assets/, PDF+PNG)...\n")
PAL <- pal_pub(name = "npg")
band_cols <- c("strong (Z>=10)" = PAL[3], "moderate (2<=Z<10)" = PAL[5], "not preserved (Z<2)" = PAL[8])

# ---- 图1. Zsummary 散点 vs 模块大小(Z=2 / Z=10 阈值线)----------------------
plt1 <- summary_df
plt1$lab <- paste0("M", plt1$module)
p1 <- ggplot(plt1, aes(size, Zsummary)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2,
           fill = PAL[8], alpha = 0.06) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           fill = PAL[5], alpha = 0.06) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = Inf,
           fill = PAL[3], alpha = 0.06) +
  geom_hline(yintercept = c(2, 10), linetype = c("dashed", "dashed"),
             colour = c(PAL[5], PAL[3]), linewidth = 0.7) +
  geom_point(aes(fill = preserved), shape = 21, size = 6, colour = "black", stroke = 0.6) +
  geom_text(aes(label = lab), fontface = "bold", size = 3.2) +
  annotate("text", x = max(plt1$size), y = 2,  label = "Z = 2",  hjust = 1.1, vjust = -0.5, size = 3, colour = PAL[5]) +
  annotate("text", x = max(plt1$size), y = 10, label = "Z = 10", hjust = 1.1, vjust = -0.5, size = 3, colour = PAL[3]) +
  scale_fill_manual(values = band_cols, name = "Preservation") +
  labs(title = "Module preservation (NetRep permutation Zsummary)",
       subtitle = "Cross-cohort: discovery -> test; permutation null = honest baseline",
       x = "Module size (genes present)", y = "Zsummary (mean Z over 7 statistics)") +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p1, file.path(ASSETS, "01_Zsummary_vs_size"), 7, 5.2)

# ---- 图2. 7 个保守性统计的 lollipop(分面 = 模块;长度=Z,颜色=p 显著)-------
zlong <- data.frame(
  module = rep(mods, times = length(stat_names)),
  stat   = rep(stat_names, each = length(mods)),
  Z      = as.vector(zmat),
  p      = as.vector(pres$p.values))
zlong$module <- factor(paste0("M", zlong$module), levels = paste0("M", mods))
zlong$stat <- factor(zlong$stat, levels = rev(stat_names))
zlong$sig <- ifelse(zlong$p < 0.05, "p < 0.05", "n.s.")
p2 <- ggplot(zlong, aes(Z, stat)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 2, linetype = "dashed", colour = PAL[5], linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = Z, y = stat, yend = stat, colour = sig), linewidth = 0.9) +
  geom_point(aes(colour = sig), size = 3) +
  facet_wrap(~ module, nrow = 1) +
  scale_colour_manual(values = c("p < 0.05" = PAL[1], "n.s." = PAL[8]), name = NULL) +
  labs(title = "Per-statistic preservation Z-scores",
       subtitle = "7 NetRep statistics; dashed = Z = 2; colour = permutation p",
       x = "Permutation Z-score", y = NULL) +
  theme_pub(base_size = 11, border = TRUE) +
  theme(legend.position = "top")
save_fig(p2, file.path(ASSETS, "02_perstat_lollipop"), 9.5, 4.6)

# ---- 图3. 置换 null 密度 + 观测值(诚实基线可视化)-------------------------
# 取代表统计 avg.cor(模块平均相关保守性);每模块画其 null 分布 + 观测竖线。
rep_stat <- "avg.cor"
si <- match(rep_stat, stat_names)
nd <- do.call(rbind, lapply(seq_along(mods), function(m) {
  data.frame(module = paste0("M", mods[m]), null = nu[m, si, ])
}))
nd$module <- factor(nd$module, levels = paste0("M", mods))
obs_df <- data.frame(module = factor(paste0("M", mods), levels = paste0("M", mods)),
                     obs = pres$observed[, si],
                     pval = pres$p.values[, si])
obs_df$lab <- sprintf("obs=%.2f\np=%.3g", obs_df$obs, obs_df$pval)
# 标签锚在各分面 null 分布的低值端(左上),避免被右侧观测线挤出画框
null_min <- tapply(nd$null, nd$module, min)
obs_df$lab_x <- null_min[as.character(obs_df$module)]
p3 <- ggplot(nd, aes(null)) +
  geom_density(aes(fill = module), alpha = 0.5, colour = NA) +
  geom_vline(data = obs_df, aes(xintercept = obs), colour = PAL[1], linewidth = 0.9) +
  geom_text(data = obs_df, aes(x = lab_x, y = Inf, label = lab),
            hjust = 0, vjust = 1.2, size = 2.7, colour = PAL[1], lineheight = 0.9) +
  facet_wrap(~ module, scales = "free", nrow = 1) +
  scale_fill_manual(values = pal_pub(length(mods), "npg"), guide = "none") +
  labs(title = sprintf("Permutation null vs observed: '%s' statistic", rep_stat),
       subtitle = "Red line = observed; null density = honest baseline. Preserved modules fall in the null tail.",
       x = paste0(rep_stat, " (preservation statistic)"), y = "Null density") +
  theme_pub(base_size = 11, border = TRUE)
save_fig(p3, file.path(ASSETS, "03_permutation_null_density"), 10, 4.2)

cat("完成。结果表见", normalizePath(args$outdir), "; 展示图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
