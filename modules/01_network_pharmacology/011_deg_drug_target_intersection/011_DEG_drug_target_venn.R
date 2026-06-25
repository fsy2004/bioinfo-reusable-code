# =============================================================================
# 编号       : 003 (通用靶点交集脚本,同用于 005/006/011)
# 脚本名     : 靶点列表交集/并集 Venn (turnkey + 顶刊图)
# 分类       : 01_network_pharmacology
# 用途       : 读取一个目录下的多份靶点/基因列表,求交集/并集,绘制 Venn(≤3集)、
#              集合大小柱状图、UpSet(≥3集),输出交集/并集表。
# 方法/包    : 集合运算 + venn_pub(零依赖,theme_pub.R)+ UpSetR
# 结果图     : Target_Venn;Set_size_bar;Target_UpSet
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 003_target_intersection_venn.R
# 运行(自己): Rscript 003_target_intersection_venn.R --input data/target_lists --outdir results/run1
# 输入规格 : --input 目录,内含多份靶点列表(csv/txt);自动识别 Gene / Gene Symbol /
#            首列为基因名。文件名 = 集合名(如 CTD / SwissTarget)。
# 整理日期 : 2026-06-23(turnkey 重构;ggvenn→venn_pub 零依赖)
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
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data"), outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# 自动取基因列
gene_col <- function(df) {
  for (c in c("Gene", "Gene.Symbol", "Gene Symbol", "gene", "Symbol")) if (c %in% names(df)) return(as.character(df[[c]]))
  as.character(df[[1]])
}
cat("Step 1/3: 读取靶点列表...\n")
files <- list.files(args$input, pattern = "\\.(csv|txt)$", full.names = TRUE, ignore.case = TRUE)
if (length(files) < 2) stop("--input 目录内需 ≥2 份靶点列表。")
sets <- lapply(files, function(f) {
  if (grepl("\\.csv$", f, ignore.case = TRUE)) v <- gene_col(read_table_smart(f))
  else { v <- trimws(readLines(f, warn = FALSE)); if (tolower(v[1]) %in% c("gene", "symbol")) v <- v[-1] }
  unique(trimws(v[!is.na(v) & v != ""]))
})
names(sets) <- tools::file_path_sans_ext(basename(files))
for (nm in names(sets)) cat("  ", nm, ":", length(sets[[nm]]), "靶点\n")

cat("Step 2/3: 交集 / 并集...\n")
uni <- Reduce(union, sets); inter <- Reduce(intersect, sets)
write.csv(data.frame(Gene = uni), file.path(args$outdir, "union_targets.csv"), row.names = FALSE)
write.csv(data.frame(Gene = inter), file.path(args$outdir, "intersection_targets.csv"), row.names = FALSE)
cat("  并集", length(uni), "· 交集", length(inter), "\n")

cat("Step 3/3: 绘图...\n")
if (length(sets) <= 3) {
  pv <- venn_pub(sets, title = "Target intersection")
  save_fig(pv, file.path(ASSETS, "Target_Venn"), 6, 6); save_fig(pv, file.path(args$outdir, "Target_Venn"), 6, 6)
}
ss <- data.frame(Set = names(sets), n = sapply(sets, length)); ss$Set <- factor(ss$Set, levels = ss$Set[order(ss$n)])
p_bar <- ggplot(ss, aes(n, Set)) +                            # lollipop(顶刊优于条形)
  geom_segment(aes(x = 0, xend = n, yend = Set, colour = Set), linewidth = 1.1) +
  geom_point(aes(colour = Set), size = 4.5) +
  geom_text(aes(label = n), hjust = -0.5, size = 3.4, fontface = "bold") +
  scale_colour_manual(values = pal_pub(nrow(ss), "npg"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Targets per set", x = "Target count", y = NULL) + theme_pub(base_size = 12, border = TRUE)
save_fig(p_bar, file.path(ASSETS, "Set_size_bar"), 6, 4); save_fig(p_bar, file.path(args$outdir, "Set_size_bar"), 6, 4)
if (length(sets) >= 3 && requireNamespace("UpSetR", quietly = TRUE)) {
  allg <- unique(unlist(sets)); mm <- as.data.frame(sapply(sets, function(g) as.integer(allg %in% g))); rownames(mm) <- allg
  for (dest in c(file.path(ASSETS, "Target_UpSet"), file.path(args$outdir, "Target_UpSet"))) {
    grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 8, height = 5)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", main.bar.color = "#3C5488", sets.bar.color = "#E64B35")); dev.off()
    grDevices::png(paste0(dest, ".png"), width = 8, height = 5, units = "in", res = 300)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", main.bar.color = "#3C5488", sets.bar.color = "#E64B35")); dev.off()
  }
}
cat("完成。交集/并集表 + 图见", normalizePath(args$outdir), "\n")
