# =============================================================================
# 编号       : 001 (通用靶点提取脚本,同用于 002/004)
# 脚本名     : 数据库靶点提取 (turnkey)
# 分类       : 01_network_pharmacology
# 用途       : 从数据库导出文件(CTD / SwissTargetPrediction / GeneCards 等)提取
#              靶点基因列表,可按评分列过滤,输出去重靶点表供下游 Venn/富集使用。
# 方法/包    : base R(自动识别基因列与评分列)
# 结果图     : 无(数据前处理;产出靶点列表)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 001_extract_targets.R
# 运行(自己): Rscript 001_extract_targets.R --input data/CTD_export.csv --score-min 0
# 可选参数 : --score-min 0(评分过滤阈值;SwissTarget 常用 0.1,GeneCards 常用 1)
# 输入规格 : 数据库导出 CSV;自动识别基因列(Gene / Gene Symbol)与评分列
#            (Probability / Relevance score / Score)。
# 整理日期 : 2026-06-23(turnkey 重构)
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
ex <- list.files(file.path(SCRIPT_DIR, "example_data"), pattern = "\\.csv$", full.names = TRUE)[1]
args <- bio_args(list(input = ex, outdir = file.path(SCRIPT_DIR, "results"), `score-min` = "0"))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
SMIN <- as.numeric(args$`score-min`)

cat("Step 1/2: 读取数据库导出:", basename(args$input), "\n")
df <- read_table_smart(args$input)
gcol <- intersect(c("Gene", "Gene.Symbol", "Gene Symbol", "Symbol"), names(df))[1]
if (is.na(gcol)) gcol <- names(df)[grep("gene|symbol", names(df), ignore.case = TRUE)][1]
if (is.na(gcol)) gcol <- names(df)[1]
scol <- intersect(c("Probability", "Relevance score", "Relevance.score", "Score", "Reference.Count"), names(df))[1]
genes <- as.character(df[[gcol]])
if (!is.na(scol) && SMIN > 0) { keep <- as.numeric(df[[scol]]) >= SMIN; genes <- genes[keep & !is.na(keep)]
  cat("  按", scol, ">=", SMIN, "过滤\n") }
genes <- unique(trimws(genes[!is.na(genes) & genes != ""]))

cat("Step 2/2: 输出靶点列表...\n")
write.csv(data.frame(Gene = genes), file.path(args$outdir, "targets.csv"), row.names = FALSE)
cat("  基因列=", gcol, if (!is.na(scol)) paste0(" · 评分列=", scol) else "", " → 去重靶点", length(genes), "个 →",
    normalizePath(file.path(args$outdir, "targets.csv")), "\n")
