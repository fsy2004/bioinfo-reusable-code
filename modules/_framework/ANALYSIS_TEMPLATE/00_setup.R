# =============================================================================
# 00_setup.R · 项目初始化:种子 / 路径 / 框架 / 断点续跑 / 可复现日志
# -----------------------------------------------------------------------------
# 用法:每个分析脚本顶部  source("00_setup.R")  即可。
# 提供:PROJ_ROOT、统一种子、目录、theme_pub、cache_step()、log_stat()、save_session()
# =============================================================================

suppressWarnings(suppressMessages({

## 1. 定位项目根 (铁律5:相对路径,禁 setwd 绝对路径) --------------------------
.this <- tryCatch({
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]])))
  else if (!is.null(sys.frame(1)$ofile)) dirname(normalizePath(sys.frame(1)$ofile))
  else getwd()
}, error = function(e) getwd())
PROJ_ROOT <<- .this

## 2. 读配置 ------------------------------------------------------------------
source(file.path(PROJ_ROOT, "config.R"))

## 3. 统一随机种子 (铁律1) ----------------------------------------------------
set.seed(SEED)
# 并行/抽样型工具(Seurat/uwot 等)在调用处也要再传 seed.use = SEED

## 4. 建标准目录 --------------------------------------------------------------
for (d in c(DIR_DATA, DIR_RESULTS, DIR_FIGURES, DIR_LOGS))
  dir.create(file.path(PROJ_ROOT, d), recursive = TRUE, showWarnings = FALSE)

## 5. 载入顶刊绘图框架 (铁律4/5:复用 theme_pub,不另写主题) -------------------
.fw <- FRAMEWORK_DIR
if (!file.exists(file.path(.fw, "theme_pub.R"))) {           # 向上自动搜 _framework
  p <- PROJ_ROOT; for (i in 1:6) { c <- file.path(p, "_framework", "theme_pub.R")
    if (file.exists(c)) { .fw <- dirname(c); break }; p <- dirname(p) }
}
if (file.exists(file.path(.fw, "theme_pub.R"))) {
  source(file.path(.fw, "theme_pub.R")); FRAMEWORK_DIR <<- .fw
  cat(sprintf("[setup] 已载入框架: %s\n", .fw))
} else cat("[setup][警告] 未找到 _framework/theme_pub.R,出图请改用矢量导出并设期刊配色\n")

}))  # end suppress

## 6. 断点续跑 (铁律5:耗时步骤幂等,产物在则跳过) -----------------------------
#' cache_step("qc", { ...重计算... })  产物存为 results/qc.rds;再次运行直接读取
cache_step <- function(name, expr, force = FALSE) {
  f <- file.path(PROJ_ROOT, DIR_RESULTS, paste0(name, ".rds"))
  if (!force && file.exists(f)) { cat(sprintf("[cache] 跳过 %s (已存在)\n", name)); return(readRDS(f)) }
  cat(sprintf("[run ] %s ...\n", name)); res <- eval.parent(substitute(expr))
  saveRDS(res, f); cat(sprintf("[done] %s -> %s\n", name, f)); res
}

## 7. 关键统计值落盘 (铁律6:数字由代码生成,不手填进文稿) ---------------------
log_stat <- function(key, value) {
  ln <- sprintf("%s\t%s", key, paste(format(value), collapse = ","))
  cat(ln, "\n", file = file.path(PROJ_ROOT, DIR_LOGS, "key_stats.tsv"), append = TRUE, sep = "")
  cat(sprintf("[stat] %s = %s\n", key, paste(format(value), collapse = ",")))
}

## 8. 可复现快照 (铁律6:记录环境) --------------------------------------------
save_session <- function() {
  writeLines(capture.output(sessionInfo()),
             file.path(PROJ_ROOT, DIR_LOGS, "sessionInfo.txt"))
  cat("[setup] sessionInfo 已写入 logs/sessionInfo.txt\n")
}

cat(sprintf("[setup] PROJ_ROOT=%s | SEED=%d | 目录就绪\n", PROJ_ROOT, SEED))
