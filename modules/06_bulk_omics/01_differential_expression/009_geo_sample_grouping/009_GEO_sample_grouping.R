# =============================================================================
# 编号       : 009
# 脚本名     : GEO 样本分组整理 + 归一化 (turnkey)
# 分类       : 03_transcriptomics_deg
# 用途       : 读取基因级表达矩阵(008 产出)与样本分组表,做重复值平均、自动 log2、
#              数组间归一化,并把分组类型作为后缀写入样本名,输出可直接进 010 的矩阵。
# 方法/包    : limma(avereps / normalizeBetweenArrays);分组表支持 csv 或 xlsx。
# 结果图     : 无(数据前处理;产出 "Sample_Type_Matrix.csv" 供 010 差异分析)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 009_GEO_sample_grouping.R
# 运行(自己): Rscript 009_GEO_sample_grouping.R --expr geneMatrix.csv --group sample_group.csv
# 输入规格 : --expr  基因级表达矩阵 CSV(首列基因名,其余样本列);
#            --group 分组表(csv/xlsx):第 1 列样本名(须与 expr 列名对应),第 2 列类型(如 con/tre)。
# 整理日期 : 2026-06-23(turnkey 重构;归一化逻辑保持原状)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages(library(limma)))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  expr  = file.path(SCRIPT_DIR, "example_data", "geneMatrix.csv"),
  group = file.path(SCRIPT_DIR, "example_data", "sample_group.csv"),
  outdir = file.path(SCRIPT_DIR, "results")))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 读表达矩阵 ----
cat("Step 1/5: 读取基因表达矩阵...\n")
em <- read_table_smart(args$expr, row_names = TRUE)
mat <- as.matrix(em); storage.mode(mat) <- "double"

# ---- Step 2. 重复值平均 + 自动 log2 + 归一化 ----
cat("Step 2/5: avereps + 自动 log2 + 数组间归一化...\n")
mat <- avereps(mat)
qx <- as.numeric(quantile(mat, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
if ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)) {
  mat[mat < 0] <- 0; mat <- log2(mat + 1); cat("  已自动 log2 转换\n")
}
mat <- normalizeBetweenArrays(mat); mat[is.na(mat)] <- 0

# ---- Step 3. 读分组表(csv / xlsx)----
cat("Step 3/5: 读取分组表...\n")
if (grepl("\\.xlsx?$", args$group, ignore.case = TRUE)) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("读取 xlsx 需要 readxl 包,或改用 csv 分组表。")
  gi <- as.data.frame(readxl::read_excel(args$group, col_names = TRUE))
} else gi <- read_table_smart(args$group)
if (ncol(gi) < 2) stop("分组表至少需两列:样本名、类型。")
snames <- as.character(gi[[1]]); stypes <- as.character(gi[[2]])

# ---- Step 4. 提取/排序样本 + 加类型后缀 ----
cat("Step 4/5: 对齐样本并加分组后缀...\n")
miss <- snames[!snames %in% colnames(mat)]
if (length(miss) > 0) { warning("表达矩阵缺失样本: ", paste(miss, collapse = ", "))
  keep <- snames %in% colnames(mat); snames <- snames[keep]; stypes <- stypes[keep] }
sub <- mat[, snames, drop = FALSE]; colnames(sub) <- paste0(snames, "_", stypes)
cat("  样本分组:", paste(names(table(stypes)), table(stypes), sep = "=", collapse = "  "), "\n")

# ---- Step 5. 输出 ----
cat("Step 5/5: 写出 Sample_Type_Matrix.csv...\n")
final <- cbind(GeneName = rownames(sub), as.data.frame(sub, check.names = FALSE))
write.csv(final, file.path(args$outdir, "Sample_Type_Matrix.csv"), row.names = FALSE)
writeLines(paste0("Number of ", names(table(stypes)), " samples: ", table(stypes)),
           file.path(args$outdir, "Sample_Summary.txt"))
cat("完成 → ", normalizePath(file.path(args$outdir, "Sample_Type_Matrix.csv")),
    "(可直接作为 010 的 --input)\n")
