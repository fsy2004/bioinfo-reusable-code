# =============================================================================
# 编号       : 008
# 脚本名     : GEO 表达矩阵整理(探针→基因)(turnkey)
# 分类       : 03_transcriptomics_deg
# 用途       : 读取 GEO 序列矩阵(GSE*.txt)与平台注释(GPL*.txt),把探针 ID 映射为
#              基因 Symbol,同一基因多探针按均值合并,输出基因级表达矩阵 geneMatrix.csv。
# 方法/包    : base R(无需第三方包);自动定位 ID_REF 起始行。
# 结果图     : 无(数据前处理模块,产出供 010/007 等下游使用)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 008_GEO_expr_matrix_tidy.R
# 运行(自己): Rscript 008_GEO_expr_matrix_tidy.R --gse GSExxx_series_matrix.txt --gpl GPLxxx.txt --symcol 11
#   或       : Rscript 008_GEO_expr_matrix_tidy.R --input data/geo_dir   (自动识别目录内 GSE*/GPL* )
# 输入规格 : GSE 序列矩阵 txt(含 "ID_REF" 表头行);GPL 平台 txt(制表符分隔,某列为基因 Symbol)。
#            --symcol = GPL 中基因 Symbol 所在列号(1 起始;不同平台不同,示例为 2)。
# 整理日期 : 2026-06-23(turnkey 重构;映射/聚合逻辑保持原状)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data"),
                      gse = "", gpl = "", symcol = "2",
                      outdir = file.path(SCRIPT_DIR, "results")))
SYMCOL <- as.integer(args$symcol)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# 自动识别 GSE / GPL 文件
if (args$gse == "" || args$gpl == "") {
  txts <- list.files(args$input, pattern = "\\.txt$", full.names = TRUE, ignore.case = TRUE)
  if (args$gse == "") args$gse <- grep("GSE", basename(txts), ignore.case = TRUE, value = FALSE) |>
    (\(i) txts[grep("GSE", basename(txts), ignore.case = TRUE)][1])()
  if (args$gpl == "") args$gpl <- txts[grep("GPL", basename(txts), ignore.case = TRUE)][1]
}
if (is.na(args$gse) || is.na(args$gpl) || args$gse == "" || args$gpl == "")
  stop("未找到 GSE/GPL 文件,请用 --gse/--gpl 指定或把文件放入 --input 目录。")
cat("Step 1/4: GSE =", basename(args$gse), " | GPL =", basename(args$gpl), "\n")

# ---- 读 GSE 序列矩阵(从 ID_REF 行起)----
lines <- readLines(args$gse)
idr <- which(grepl("ID_REF", lines))[1]
if (is.na(idr)) stop("GSE 文件中未找到 ID_REF 行。")
ed <- read.delim(args$gse, header = FALSE, sep = "\t", quote = "\"", skip = idr - 1, comment.char = "")
colnames(ed) <- as.character(unlist(ed[1, ])); ed <- ed[-1, , drop = FALSE]; colnames(ed)[1] <- "ProbeID"
samples <- setdiff(colnames(ed), "ProbeID")
for (c in samples) ed[[c]] <- as.numeric(as.character(ed[[c]]))
cat("Step 2/4: 表达矩阵", nrow(ed), "探针 x", length(samples), "样本\n")

# ---- 读 GPL 平台,建 探针→Symbol 映射 ----
pl <- read.delim(args$gpl, header = FALSE, sep = "\t", quote = "\"", comment.char = "#", stringsAsFactors = FALSE)
keep <- !grepl("^(ID|!)", pl[[1]]) & pl[[1]] != "" & ncol(pl) >= SYMCOL
sym <- ifelse(keep, sub("(.+?)///.*", "\\1", gsub('"', '', pl[[SYMCOL]])), "")
sym[grepl("\\s", sym)] <- ""          # 含空格的非单词 Symbol 丢弃(与原逻辑一致)
map <- data.frame(ProbeID = pl[[1]], geneSymbol = sym, stringsAsFactors = FALSE)
map <- map[map$geneSymbol != "", ]
cat("Step 3/4: 建立探针→基因映射", nrow(map), "条\n")

# ---- 合并 + 按基因取均值 ----
merged <- merge(ed, map, by = "ProbeID")
agg <- aggregate(merged[, samples, drop = FALSE], by = list(geneSymbol = merged$geneSymbol),
                 FUN = function(x) mean(x, na.rm = TRUE))
agg <- agg[order(agg$geneSymbol), ]
out <- file.path(args$outdir, "geneMatrix.csv")
write.csv(agg, out, row.names = FALSE)
cat("Step 4/4: 输出基因级矩阵", nrow(agg), "基因 →", normalizePath(out), "\n")
