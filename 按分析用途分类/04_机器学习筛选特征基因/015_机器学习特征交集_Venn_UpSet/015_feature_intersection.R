# =============================================================================
# 编号       : 015
# 脚本名     : 多方法特征交集 Venn / UpSet (turnkey + 顶刊图)
# 分类       : 04_机器学习筛选特征基因
# 用途       : 读取一个目录下的多份基因列表(各 ML 方法的特征),求交集并绘制
#              Venn(≤3 集)与 UpSet 图,输出全局交集与两两交集表。
# 方法/包    : 集合运算 + venn_pub(零依赖,theme_pub.R)+ UpSetR
# 结果图     : Feature_Venn(≤3集);Feature_UpSet(任意集合数)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 015_feature_intersection.R
# 运行(自己): Rscript 015_feature_intersection.R --input data/gene_sets_dir --outdir results/run1
# 输入规格 : --input 指向一个【目录】,内含多份基因列表(csv/txt,首列=基因名;
#            csv 首列名建议 gene)。每份文件名 = 集合名。
# 整理日期 : 2026-06-23(turnkey 重构;ggvenn→venn_pub 零依赖实现)
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
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "gene_sets"),
                      outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- 读取目录下所有基因列表 ----
cat("Step 1/3: 读取基因列表...\n")
files <- list.files(args$input, pattern = "\\.(csv|txt)$", full.names = TRUE, ignore.case = TRUE)
if (length(files) < 2) stop("--input 目录内需 ≥2 份基因列表。")
sets <- lapply(files, function(f) {
  if (grepl("\\.csv$", f, ignore.case = TRUE)) v <- as.character(read_table_smart(f)[[1]])
  else { v <- trimws(readLines(f, warn = FALSE)); if (tolower(v[1]) %in% c("gene", "variable")) v <- v[-1] }
  unique(trimws(v[v != ""]))
})
names(sets) <- tools::file_path_sans_ext(basename(files))
for (nm in names(sets)) cat("  ", nm, ":", length(sets[[nm]]), "基因\n")

# ---- 交集表 ----
cat("Step 2/3: 计算交集...\n")
inter <- Reduce(intersect, sets)
write.csv(data.frame(Gene = inter), file.path(args$outdir, "global_intersection.csv"), row.names = FALSE)
cat("  全局交集:", length(inter), "基因\n")
pw <- do.call(rbind, combn(names(sets), 2, simplify = FALSE) |> lapply(function(p) {
  g <- intersect(sets[[p[1]]], sets[[p[2]]])
  data.frame(SetA = p[1], SetB = p[2], n = length(g), genes = paste(g, collapse = ";"))
}))
write.csv(pw, file.path(args$outdir, "pairwise_intersection.csv"), row.names = FALSE)

# ---- 图 ----
cat("Step 3/3: 绘图...\n")
if (length(sets) <= 3) {
  pv <- venn_pub(sets, title = "Feature intersection (Venn)")
  save_fig(pv, file.path(ASSETS, "Feature_Venn"), 6, 6); save_fig(pv, file.path(args$outdir, "Feature_Venn"), 6, 6)
}
if (requireNamespace("UpSetR", quietly = TRUE) && length(sets) >= 2) {
  allg <- unique(unlist(sets)); mm <- as.data.frame(sapply(sets, function(g) as.integer(allg %in% g))); rownames(mm) <- allg
  for (dest in c(file.path(ASSETS, "Feature_UpSet"), file.path(args$outdir, "Feature_UpSet"))) {
    grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 9, height = 5.5)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", point.size = 2.8, line.size = 1,
                        main.bar.color = "#3C5488", sets.bar.color = "#E64B35",
                        mainbar.y.label = "Intersection size", sets.x.label = "Genes per set")); dev.off()
    grDevices::png(paste0(dest, ".png"), width = 9, height = 5.5, units = "in", res = 300)
    print(UpSetR::upset(mm, nsets = ncol(mm), order.by = "freq", point.size = 2.8, line.size = 1,
                        main.bar.color = "#3C5488", sets.bar.color = "#E64B35",
                        mainbar.y.label = "Intersection size", sets.x.label = "Genes per set")); dev.off()
  }
}
cat("完成。交集表/图见", normalizePath(args$outdir), "\n")
