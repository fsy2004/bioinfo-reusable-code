# =============================================================================
# 编号       : 035
# 脚本名     : 多方法组合交集选择 (turnkey + 顶刊图)
# 分类       : 04_机器学习筛选特征基因
# 用途       : 在多种 ML 方法的特征列表中,枚举所有 k 方法组合的交集大小,排序选优,
#              输出最佳组合的交集基因 + UpSet 图 + 组合排行榜。
# 方法/包    : 集合运算 + 组合枚举;UpSetR;绘图 theme_pub.R
# 结果图     : Combo_ranking(top 组合交集大小);Combo_UpSet(选定组合的特征交集)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 035_method_combo_intersection.R
# 运行(自己): Rscript 035_method_combo_intersection.R --input data/method_sets --pick 5
# 可选参数 : --pick 5(组合的方法数) --methods "RF,Lasso,SVM,..."(指定组合,覆盖枚举)
# 输入规格 : --input 目录,内含多份方法特征列表(csv,列名 variable 或首列=基因)。
# 整理日期 : 2026-06-23(turnkey 重构;ggvenn→venn_pub/UpSet)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages(library(ggplot2)))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "method_sets"),
                      outdir = file.path(SCRIPT_DIR, "results"), pick = "5", methods = ""))
PICK <- as.integer(args$pick)
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- 读方法特征列表 ----
cat("Step 1/3: 读取各方法特征...\n")
files <- list.files(args$input, pattern = "\\.csv$", full.names = TRUE)
if (length(files) < PICK) stop(sprintf("方法数(%d)少于 --pick(%d)。", length(files), PICK))
sets <- lapply(files, function(f) { d <- read_table_smart(f); v <- if ("variable" %in% names(d)) d$variable else d[[1]]; unique(trimws(as.character(v))) })
names(sets) <- sub("^importanceGene\\.", "", tools::file_path_sans_ext(basename(files)))
cat("  方法数:", length(sets), "\n")

# ---- 枚举 k 组合交集 ----
cat("Step 2/3: 枚举", PICK, "方法组合交集...\n")
combos <- combn(names(sets), PICK, simplify = FALSE)
tab <- do.call(rbind, lapply(combos, function(cb) {
  data.frame(combination = paste(cb, collapse = " + "), n_intersect = length(Reduce(intersect, sets[cb])))
}))
tab <- tab[order(-tab$n_intersect), ]
write.csv(tab, file.path(args$outdir, "all_combinations.csv"), row.names = FALSE)
sel <- if (nzchar(args$methods)) trimws(strsplit(args$methods, ",")[[1]]) else combos[[which.max(vapply(combos, function(cb) length(Reduce(intersect, sets[cb])), 0))]]
sel <- intersect(sel, names(sets))
inter <- Reduce(intersect, sets[sel])
writeLines(inter, file.path(args$outdir, "selected_combo_intersection.txt"))
cat("  选定组合:", paste(sel, collapse = "+"), "→ 交集", length(inter), "基因\n")

# ---- 图 1: top 组合排行榜 ----
cat("Step 3/3: 绘图...\n")
topn <- head(tab, min(12, nrow(tab))); topn$combination <- factor(topn$combination, levels = rev(topn$combination))
p_rank <- ggplot(topn, aes(n_intersect, combination, fill = n_intersect)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_intersect), hjust = -0.3, size = 3.2, fontface = "bold") +
  scale_fill_viridis_c(option = "D", guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = paste0("Top ", PICK, "-method combinations"), x = "Intersection size", y = NULL) +
  theme_pub(base_size = 11, border = TRUE) + theme(axis.text.y = element_text(size = 8))
save_fig(p_rank, file.path(ASSETS, "Combo_ranking"), 7.5, 5); save_fig(p_rank, file.path(args$outdir, "Combo_ranking"), 7.5, 5)

# ---- 图 2: 选定组合 UpSet ----
if (requireNamespace("UpSetR", quietly = TRUE)) {
  allg <- unique(unlist(sets[sel])); mm <- as.data.frame(sapply(sets[sel], function(g) as.integer(allg %in% g))); rownames(mm) <- allg
  for (dest in c(file.path(ASSETS, "Combo_UpSet"), file.path(args$outdir, "Combo_UpSet"))) {
    grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 9, height = 5.5)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", point.size = 2.8, line.size = 1,
                        main.bar.color = "#3C5488", sets.bar.color = "#E64B35",
                        mainbar.y.label = "Intersection size", sets.x.label = "Features per method")); dev.off()
    grDevices::png(paste0(dest, ".png"), width = 9, height = 5.5, units = "in", res = 300)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", point.size = 2.8, line.size = 1,
                        main.bar.color = "#3C5488", sets.bar.color = "#E64B35",
                        mainbar.y.label = "Intersection size", sets.x.label = "Features per method")); dev.off()
  }
}
cat("完成。组合表/交集/图见", normalizePath(args$outdir), "\n")
