# ==========================================================================
# 脚本名     : scFOCAL官方GUI输入准备.R
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 从压缩包 489.scFOCAL分析.rar 整理
# 原始文件   : 489.scFOCAL分析\scFOCAL分析.R
# 用途       : 读取 Seurat 对象与细胞ID-细胞类型映射表，生成 scFOCAL GUI 可上传的 RDS，并启动 scFOCAL Shiny交互界面。
# 结果图     : scFOCAL官方交互页面输出；本脚本主要准备RDS输入并启动GUI
# 非肿瘤消化适配: 适合。可作为非肿瘤消化系统单细胞功能定位/细胞类型解释的交互分析入口。
# 主要 R 包  : Seurat; scFOCAL; shiny
# 整理日期   : 2026-05-13
# 备注       : 保留bioinfo-reusable-code逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# -*- coding: utf-8 -*-
################################################################################
# scFOCAL 官方交互页面启动脚本
#
# 输入文件：
#   01_SingleCell_Analysis_Complete.RData
#   03_Cell_ID_to_CellType_Mapping.csv
#
# 功能：
#   1. 从 RData 中读取 Seurat 对象。
#   2. 将细胞注释表合并到 Seurat metadata。
#   3. 保存一个 scFOCAL 官方 GUI 可上传的 RDS 文件。
#   4. 启动 scFOCAL::runscFOCAL()，弹出官方 Shiny 交互页面。
################################################################################

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")

open_in_r_viewer <- function(url) {
  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()
  ) {
    rstudioapi::viewer(url)
  } else {
    utils::browseURL(url)
  }
}

workDir <- "H:/常用分析生信/489.scFOCAL分析"
setwd(workDir)

input_rdata <- file.path(workDir, "01_SingleCell_Analysis_Complete.RData")
cell_mapping_file <- file.path(workDir, "03_Cell_ID_to_CellType_Mapping.csv")
prepared_rds <- file.path(workDir, "scFOCAL_input_seurat_with_celltype.rds")

required_packages <- c("Seurat", "scFOCAL", "shiny")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop(
    "缺少以下 R 包，请先安装后重新运行：",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(scFOCAL)
  library(shiny)
})

message("读取 Seurat RData：", input_rdata)
env <- new.env(parent = emptyenv())
loaded_names <- load(input_rdata, envir = env)
seurat_names <- loaded_names[
  vapply(loaded_names, function(nm) inherits(env[[nm]], "Seurat"), logical(1))
]

if (length(seurat_names) != 1) {
  stop("RData 中必须且只能包含 1 个 Seurat 对象。当前检测到：", length(seurat_names), call. = FALSE)
}

seurat_obj <- env[[seurat_names]]
message("检测到 Seurat 对象：", seurat_names)

message("读取细胞注释：", cell_mapping_file)
cell_mapping <- read.csv(cell_mapping_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
required_cols <- c("Cell_ID", "Sample", "Cluster_ID", "Cell_Type", "Group")
missing_cols <- setdiff(required_cols, colnames(cell_mapping))
if (length(missing_cols) > 0) {
  stop("细胞注释文件缺少列：", paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(cell_mapping$Cell_ID) > 0) {
  stop("细胞注释文件 Cell_ID 存在重复值。", call. = FALSE)
}

cell_mapping <- cell_mapping[match(colnames(seurat_obj), cell_mapping$Cell_ID), , drop = FALSE]
if (anyNA(cell_mapping$Cell_ID)) {
  stop("Seurat 细胞 ID 与注释表 Cell_ID 不能完全匹配。", call. = FALSE)
}

rownames(cell_mapping) <- cell_mapping$Cell_ID

# 保留官方 GUI 最容易识别和选择的 metadata 列名。
metadata_to_add <- data.frame(
  Cell_ID = cell_mapping$Cell_ID,
  Sample = cell_mapping$Sample,
  Cluster_ID = cell_mapping$Cluster_ID,
  Cell_Type = cell_mapping$Cell_Type,
  Group = cell_mapping$Group,
  scFOCAL_Group_CellType = paste(cell_mapping$Group, cell_mapping$Cell_Type, sep = "__"),
  row.names = cell_mapping$Cell_ID,
  check.names = FALSE
)

seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_to_add)

message("保存 scFOCAL GUI 输入 RDS：", prepared_rds)
saveRDS(seurat_obj, prepared_rds)

options(shiny.launch.browser = open_in_r_viewer)

message("启动 scFOCAL 官方交互页面...")
scFOCAL::runscFOCAL()
