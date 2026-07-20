# ==========================================================================
# 脚本名     : 单细胞Scissor疾病相关细胞筛选.R
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 从压缩包 464.Scissor算法.rar 整理
# 原始文件   : 464.Scissor算法\Scissor算法.R
# 用途       : 将 bulk 表型信息与已注释 Seurat 对象关联，用 Scissor 识别疾病相关和对照相关细胞，并汇总细胞类型分布。
# 结果图     : Scissor+/- UMAP图；Scissor细胞类型分布柱状图；Scissor分组汇总表；细胞映射表
# 非肿瘤消化适配: 很适合。可把非肿瘤消化系统bulk队列与单细胞图谱连接起来，是有新意的桥接模块。
# 主要 R 包  : Seurat; ggplot2; Matrix; Scissor
# 整理日期   : 2026-05-13
# 备注       : 保留bioinfo-reusable-code逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
################################################################################
#                      Annotated scRNA + Scissor Analysis                      #
#                    已注释单细胞对象的 Scissor 分析流程                       #
################################################################################

# 设置工作目录
# 如需修改路径，只改下面这一行即可
workDir <- "H:\\常用分析生信\\464.Scissor算法"
workDir <- normalizePath(workDir, winslash = "\\", mustWork = TRUE)

# 加载所需R包
library(Seurat)
library(ggplot2)
library(Matrix)
library(Scissor)

infer_group_text <- function(values) {
  values <- as.character(values)
  values <- trimws(values)
  values[values == ""] <- NA_character_

  is_control <- grepl(
    "(^|[._-])(con|control|ctrl|normal|healthy)([._-]|$)",
    values,
    ignore.case = TRUE,
    perl = TRUE
  )
  is_disease <- grepl(
    "(^|[._-])(tre|treat|treated|disease|case|patient|tumou?r)([._-]|$)",
    values,
    ignore.case = TRUE,
    perl = TRUE
  )

  ifelse(is_control, "control", ifelse(is_disease, "Disease", NA_character_))
}

normalize_group_labels <- function(primary_values, fallback_values = NULL, context_label = "group") {
  normalized <- infer_group_text(primary_values)

  if (!is.null(fallback_values)) {
    missing_idx <- is.na(normalized)
    if (any(missing_idx)) {
      normalized[missing_idx] <- infer_group_text(fallback_values[missing_idx])
    }
  }

  if (any(is.na(normalized))) {
    bad_values <- unique(as.character(primary_values[is.na(normalized)]))
    stop(sprintf(
      "Failed to infer %s labels. Example values: %s / 无法识别%s分组，示例值：%s",
      context_label,
      paste(head(bad_values, 10), collapse = ", "),
      context_label,
      paste(head(bad_values, 10), collapse = ", ")
    ))
  }

  factor(normalized, levels = c("control", "Disease"))
}

collapse_duplicated_bulk_genes <- function(expr_matrix) {
  if (is.null(rownames(expr_matrix))) {
    stop("Bulk expression matrix must contain gene names in rownames / bulk表达矩阵必须在行名中包含基因名")
  }

  gene_names <- rownames(expr_matrix)
  duplicated_n <- sum(duplicated(gene_names))
  if (duplicated_n == 0) {
    return(expr_matrix)
  }

  message(sprintf("Detected %d duplicated genes in bulk matrix; collapsing duplicates by row mean.", duplicated_n))
  message(sprintf("bulk矩阵检测到 %d 个重复基因，按行均值进行合并。", duplicated_n))

  gene_counts <- table(gene_names)
  collapsed_matrix <- rowsum(expr_matrix, group = gene_names, reorder = FALSE)
  collapsed_matrix <- collapsed_matrix / as.numeric(gene_counts[rownames(collapsed_matrix)])
  collapsed_matrix
}

infer_bulk_phenotype <- function(sample_names) {
  phenotype_group <- normalize_group_labels(
    primary_values = sample_names,
    context_label = "bulk sample"
  )
  phenotype <- ifelse(phenotype_group == "control", 0L, 1L)
  names(phenotype) <- sample_names
  phenotype
}

normalize_scissor_selection <- function(selected_cells, all_cell_ids) {
  if (is.null(selected_cells) || length(selected_cells) == 0) {
    return(character(0))
  }

  if (is.numeric(selected_cells)) {
    selected_cells <- all_cell_ids[as.integer(selected_cells)]
  }

  selected_cells <- as.character(selected_cells)
  intersect(selected_cells, all_cell_ids)
}

extract_meta_column <- function(seurat_obj, column_name) {
  if (column_name %in% colnames(seurat_obj@meta.data)) {
    return(seurat_obj@meta.data[[column_name]])
  }
  rep(NA, ncol(seurat_obj))
}

make_scissor_compatible <- function(seurat_obj, assay_name = "RNA") {
  if (!assay_name %in% names(seurat_obj@assays)) {
    stop(sprintf(
      "Assay '%s' not found in Seurat object / Seurat对象中未找到assay：%s",
      assay_name,
      assay_name
    ))
  }

  current_assay <- seurat_obj[[assay_name]]
  if (!inherits(current_assay, "Assay5")) {
    return(seurat_obj)
  }

  message(sprintf("Converting %s from Assay5 to legacy Assay for Scissor compatibility...", assay_name))
  message(sprintf("将 %s 从 Assay5 转为旧版 Assay，以兼容 Scissor...", assay_name))

  counts_mat <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  data_mat <- GetAssayData(seurat_obj, assay = assay_name, layer = "data")

  legacy_assay <- SeuratObject::CreateAssayObject(counts = counts_mat)
  legacy_assay@data <- data_mat
  legacy_assay@key <- current_assay@key

  variable_features <- tryCatch(
    VariableFeatures(seurat_obj[[assay_name]]),
    error = function(e) character(0)
  )
  if (length(variable_features) > 0) {
    VariableFeatures(legacy_assay) <- variable_features
  }

  seurat_obj[[assay_name]] <- legacy_assay
  seurat_obj
}

summarize_selected_cells <- function(data_frame, group_col, output_col_name) {
  if (nrow(data_frame) == 0) {
    out <- data.frame(
      Scissor = character(0),
      Placeholder = character(0),
      Cell_Count = integer(0),
      Proportion = numeric(0),
      stringsAsFactors = FALSE
    )
    colnames(out)[2] <- output_col_name
    return(out)
  }

  tab <- as.data.frame(
    table(
      Scissor = data_frame$Scissor,
      Category = data_frame[[group_col]]
    ),
    stringsAsFactors = FALSE
  )
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  colnames(tab) <- c("Scissor", output_col_name, "Cell_Count")
  tab$Proportion <- ave(
    tab$Cell_Count,
    tab$Scissor,
    FUN = function(x) round(x / sum(x), 4)
  )
  tab
}

annotated_rdata_file <- file.path(workDir, "01_SingleCell_Analysis_Complete.RData")
mapping_file <- file.path(workDir, "03_Cell_ID_to_CellType_Mapping.csv")
bulk_input_file <- file.path(workDir, "Sample Type Matrix.csv")

output_dirs <- list(
  scissor_analysis = "08.Scissor_Analysis",
  final_results = "16.Final_Results"
)

message("Required packages loaded successfully / 必需R包加载完成")

if (!dir.exists(workDir)) {
  stop("Working directory does not exist / 工作目录不存在：", workDir)
}

setwd(workDir)
message("Working directory set to / 工作目录设置为：", getwd())

for (dir_name in output_dirs) {
  full_path <- file.path(workDir, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
  }
}

if (!file.exists(annotated_rdata_file)) {
  stop("Annotated RData file not found / 未找到已注释RData文件：", annotated_rdata_file)
}
if (!file.exists(mapping_file)) {
  stop("Cell annotation mapping file not found / 未找到细胞注释映射文件：", mapping_file)
}
if (!file.exists(bulk_input_file)) {
  stop("Bulk expression matrix not found / 未找到bulk表达矩阵文件：", bulk_input_file)
}

message("Loading annotated single-cell object...")
message("读取已注释单细胞对象...")

loaded_objects <- load(annotated_rdata_file)
seurat_object_names <- loaded_objects[vapply(
  loaded_objects,
  function(current_name) inherits(get(current_name), "Seurat"),
  logical(1)
)]

if (length(seurat_object_names) == 0) {
  stop("No Seurat object found in the RData file / RData中未找到Seurat对象")
}

seurat_object_name <- if ("Peripheral_Blood_Mononuclear_Cells" %in% seurat_object_names) {
  "Peripheral_Blood_Mononuclear_Cells"
} else {
  seurat_object_names[1]
}

Peripheral_Blood_Mononuclear_Cells <- get(seurat_object_name)

message(sprintf("Using Seurat object: %s", seurat_object_name))
message(sprintf("使用Seurat对象：%s", seurat_object_name))
message(sprintf("Seurat object dimensions: %d genes x %d cells", nrow(Peripheral_Blood_Mononuclear_Cells), ncol(Peripheral_Blood_Mononuclear_Cells)))
message(sprintf("Seurat对象维度：%d 基因 x %d 细胞", nrow(Peripheral_Blood_Mononuclear_Cells), ncol(Peripheral_Blood_Mononuclear_Cells)))

if (is.null(rownames(Peripheral_Blood_Mononuclear_Cells)) || nrow(Peripheral_Blood_Mononuclear_Cells) == 0) {
  stop("The Seurat object does not contain valid gene names / Seurat对象不包含有效基因名")
}

if (length(Peripheral_Blood_Mononuclear_Cells@graphs) == 0) {
  stop(
    "The Seurat object does not contain a cell-cell network required by Scissor / ",
    "Seurat对象中缺少Scissor所需的细胞网络信息"
  )
}

mapping_df <- read.csv(
  file = mapping_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_mapping_cols <- c("Cell_ID", "Cell_Type", "Group")
missing_mapping_cols <- setdiff(required_mapping_cols, colnames(mapping_df))
if (length(missing_mapping_cols) > 0) {
  stop(sprintf(
    "Mapping file is missing required columns: %s / 映射文件缺少必要列：%s",
    paste(missing_mapping_cols, collapse = ", "),
    paste(missing_mapping_cols, collapse = ", ")
  ))
}

if (anyDuplicated(mapping_df$Cell_ID) > 0) {
  stop("Duplicated Cell_ID values found in mapping file / 映射文件中存在重复Cell_ID")
}

obj_cells <- colnames(Peripheral_Blood_Mononuclear_Cells)
if (!setequal(obj_cells, mapping_df$Cell_ID)) {
  obj_only <- setdiff(obj_cells, mapping_df$Cell_ID)
  map_only <- setdiff(mapping_df$Cell_ID, obj_cells)
  stop(sprintf(
    paste0(
      "Cell IDs differ between Seurat object and mapping file. ",
      "Only in object: %s ; only in mapping: %s / ",
      "Seurat对象与映射文件的细胞ID不一致。对象独有：%s；映射独有：%s"
    ),
    paste(head(obj_only, 10), collapse = ", "),
    paste(head(map_only, 10), collapse = ", "),
    paste(head(obj_only, 10), collapse = ", "),
    paste(head(map_only, 10), collapse = ", ")
  ))
}

mapping_df <- mapping_df[match(obj_cells, mapping_df$Cell_ID), , drop = FALSE]

Peripheral_Blood_Mononuclear_Cells$cell_type <- as.character(mapping_df$Cell_Type)
Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- as.character(mapping_df$Cell_Type)
Peripheral_Blood_Mononuclear_Cells$cellType <- as.character(mapping_df$Cell_Type)
Peripheral_Blood_Mononuclear_Cells$Type <- normalize_group_labels(
  primary_values = mapping_df$Group,
  fallback_values = mapping_df$Cell_ID,
  context_label = "single-cell"
)
Peripheral_Blood_Mononuclear_Cells$group <- as.character(Peripheral_Blood_Mononuclear_Cells$Type)
Peripheral_Blood_Mononuclear_Cells$sample_id <- sub("_[^_]+$", "", obj_cells)
Peripheral_Blood_Mononuclear_Cells$Sample <- Peripheral_Blood_Mononuclear_Cells$sample_id

if ("Cluster_ID" %in% colnames(mapping_df)) {
  Peripheral_Blood_Mononuclear_Cells$seurat_clusters <- factor(as.character(mapping_df$Cluster_ID))
}

if (any(is.na(Peripheral_Blood_Mononuclear_Cells$cell_type)) ||
    any(trimws(Peripheral_Blood_Mononuclear_Cells$cell_type) == "")) {
  stop("Mapping file contains missing cell type annotations / 映射文件中存在缺失细胞类型注释")
}

Peripheral_Blood_Mononuclear_Cells <- make_scissor_compatible(
  seurat_obj = Peripheral_Blood_Mononuclear_Cells,
  assay_name = "RNA"
)
DefaultAssay(Peripheral_Blood_Mononuclear_Cells) <- "RNA"

Idents(Peripheral_Blood_Mononuclear_Cells) <- Peripheral_Blood_Mononuclear_Cells$cell_type

write.csv(
  mapping_df,
  file = file.path(output_dirs$scissor_analysis, "00_Input_Cell_Annotation_Mapping.csv"),
  row.names = FALSE
)

bulk_input <- read.csv(
  file = bulk_input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

if (ncol(bulk_input) < 2) {
  stop("Bulk matrix must contain one gene column and sample columns / bulk矩阵至少需要1列基因名和若干样本列")
}

bulk_gene_names <- trimws(as.character(bulk_input[[1]]))
bulk_gene_names <- sub("^\ufeff", "", bulk_gene_names)
valid_bulk_rows <- !is.na(bulk_gene_names) & bulk_gene_names != ""
bulk_input <- bulk_input[valid_bulk_rows, , drop = FALSE]
bulk_gene_names <- bulk_gene_names[valid_bulk_rows]

bulk_dataset <- as.matrix(bulk_input[, -1, drop = FALSE])
storage.mode(bulk_dataset) <- "numeric"
rownames(bulk_dataset) <- bulk_gene_names

if (anyNA(bulk_dataset)) {
  stop("Bulk matrix contains non-numeric values or missing values / bulk矩阵中存在非数值或缺失值")
}

bulk_dataset <- collapse_duplicated_bulk_genes(bulk_dataset)
phenotype <- infer_bulk_phenotype(colnames(bulk_dataset))
tag <- c("control", "disease")

bulk_group_info <- data.frame(
  Sample_ID = colnames(bulk_dataset),
  Group = ifelse(phenotype == 0L, "control", "disease"),
  Phenotype = phenotype,
  stringsAsFactors = FALSE
)

write.csv(
  bulk_group_info,
  file = file.path(output_dirs$scissor_analysis, "01_Bulk_Sample_Group_Info.csv"),
  row.names = FALSE
)

common_genes <- intersect(rownames(bulk_dataset), rownames(Peripheral_Blood_Mononuclear_Cells))
if (length(common_genes) < 1000) {
  stop(sprintf(
    "Too few common genes between bulk and single-cell datasets: %d / bulk与单细胞共同基因数过少：%d",
    length(common_genes),
    length(common_genes)
  ))
}

message(sprintf("Bulk matrix dimensions: %d genes x %d samples", nrow(bulk_dataset), ncol(bulk_dataset)))
message(sprintf("bulk矩阵维度：%d 基因 x %d 样本", nrow(bulk_dataset), ncol(bulk_dataset)))
message(sprintf("Common genes: %d", length(common_genes)))
message(sprintf("共同基因数：%d", length(common_genes)))

scissor_alpha <- 0.2
scissor_input_save <- file.path(output_dirs$scissor_analysis, "02_Scissor_Model_Input.RData")

# phenotype 编码:
# 0 = control
# 1 = disease
# 因此 Scissor+ 更偏向疾病组，Scissor- 更偏向对照组
message("Running Scissor...")
message("开始运行Scissor...")

scissor_result <- Scissor::Scissor(
  bulk_dataset,
  Peripheral_Blood_Mononuclear_Cells,
  phenotype,
  tag = tag,
  alpha = scissor_alpha,
  family = "binomial",
  Save_file = scissor_input_save
)

saveRDS(
  scissor_result,
  file = file.path(output_dirs$scissor_analysis, "03_Scissor_Result.rds")
)

all_cell_ids <- colnames(Peripheral_Blood_Mononuclear_Cells)
scissor_pos_cells <- normalize_scissor_selection(scissor_result$Scissor_pos, all_cell_ids)
scissor_neg_cells <- normalize_scissor_selection(scissor_result$Scissor_neg, all_cell_ids)

scissor_label <- setNames(rep("Background", length(all_cell_ids)), all_cell_ids)
scissor_label[scissor_pos_cells] <- "Scissor+"
scissor_label[scissor_neg_cells] <- "Scissor-"
scissor_label <- factor(scissor_label, levels = c("Background", "Scissor+", "Scissor-"))

scissor_association <- ifelse(
  scissor_label == "Scissor+",
  "Disease_Associated",
  ifelse(scissor_label == "Scissor-", "Control_Associated", "Background")
)

Peripheral_Blood_Mononuclear_Cells$scissor <- scissor_label
Peripheral_Blood_Mononuclear_Cells$scissor_association <- scissor_association

selected_cells_n <- length(scissor_pos_cells) + length(scissor_neg_cells)
selected_ratio <- round(selected_cells_n / ncol(Peripheral_Blood_Mononuclear_Cells) * 100, 3)

scissor_summary <- data.frame(
  Metric = c(
    "Bulk_Samples_Total",
    "Bulk_Control_Samples",
    "Bulk_Disease_Samples",
    "SingleCell_Cells_Total",
    "SingleCell_Control_Cells",
    "SingleCell_Disease_Cells",
    "Common_Genes",
    "Scissor_Pos_Cells",
    "Scissor_Neg_Cells",
    "Selected_Cells_Total",
    "Selected_Cell_Percentage",
    "Alpha",
    "Model_Family"
  ),
  Value = c(
    ncol(bulk_dataset),
    sum(phenotype == 0L),
    sum(phenotype == 1L),
    ncol(Peripheral_Blood_Mononuclear_Cells),
    sum(Peripheral_Blood_Mononuclear_Cells$Type == "control"),
    sum(Peripheral_Blood_Mononuclear_Cells$Type == "Disease"),
    length(common_genes),
    length(scissor_pos_cells),
    length(scissor_neg_cells),
    selected_cells_n,
    paste0(selected_ratio, "%"),
    scissor_alpha,
    "binomial"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scissor_summary,
  file = file.path(output_dirs$scissor_analysis, "04_Scissor_Summary.csv"),
  row.names = FALSE
)

cell_mapping <- data.frame(
  Cell_ID = all_cell_ids,
  Sample_ID = Peripheral_Blood_Mononuclear_Cells$sample_id,
  Cluster_ID = as.character(extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "seurat_clusters")),
  Cell_Type = as.character(Peripheral_Blood_Mononuclear_Cells$cell_type),
  Group = as.character(Peripheral_Blood_Mononuclear_Cells$Type),
  Scissor = as.character(Peripheral_Blood_Mononuclear_Cells$scissor),
  Scissor_Association = as.character(Peripheral_Blood_Mononuclear_Cells$scissor_association),
  nFeature_RNA = extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "nFeature_RNA"),
  nCount_RNA = extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "nCount_RNA"),
  percent_mt = extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "percent.mt"),
  decontX_contamination = extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "decontX_contamination"),
  decontX_threshold = extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "decontX_threshold"),
  Doublet_Classification = extract_meta_column(Peripheral_Blood_Mononuclear_Cells, "Doublet_Classification"),
  stringsAsFactors = FALSE
)

write.csv(
  cell_mapping,
  file = file.path(output_dirs$scissor_analysis, "05_Scissor_Cell_Results.csv"),
  row.names = FALSE
)

selected_cell_mapping <- cell_mapping[cell_mapping$Scissor != "Background", , drop = FALSE]
scissor_celltype_summary <- summarize_selected_cells(selected_cell_mapping, "Cell_Type", "Cell_Type")
scissor_group_summary <- summarize_selected_cells(selected_cell_mapping, "Group", "Group")

write.csv(
  scissor_celltype_summary,
  file = file.path(output_dirs$scissor_analysis, "06_Scissor_CellType_Summary.csv"),
  row.names = FALSE
)

write.csv(
  scissor_group_summary,
  file = file.path(output_dirs$scissor_analysis, "07_Scissor_Group_Summary.csv"),
  row.names = FALSE
)

if ("umap" %in% names(Peripheral_Blood_Mononuclear_Cells@reductions)) {
  scissor_colors <- c(
    "Background" = "grey85",
    "Scissor+" = "#E64B35",
    "Scissor-" = "#4DBBD5"
  )

  umap_embeddings <- Embeddings(Peripheral_Blood_Mononuclear_Cells, reduction = "umap")
  umap_plot_df <- data.frame(
    UMAP_1 = umap_embeddings[, 1],
    UMAP_2 = umap_embeddings[, 2],
    Scissor = as.character(Peripheral_Blood_Mononuclear_Cells$scissor),
    Type = as.character(Peripheral_Blood_Mononuclear_Cells$Type),
    stringsAsFactors = FALSE
  )
  umap_plot_df$Scissor <- factor(
    umap_plot_df$Scissor,
    levels = c("Background", "Scissor+", "Scissor-")
  )
  draw_order <- c("Background", "Scissor-", "Scissor+")
  umap_plot_df <- umap_plot_df[order(match(umap_plot_df$Scissor, draw_order)), , drop = FALSE]

  pdf(file = file.path(output_dirs$scissor_analysis, "08_Scissor_UMAP.pdf"),
      width = 10, height = 8)
  print(
    ggplot(umap_plot_df, aes(x = UMAP_1, y = UMAP_2, color = Scissor)) +
    geom_point(size = 1, alpha = 0.9) +
    scale_color_manual(
      values = scissor_colors,
      breaks = c("Background", "Scissor+", "Scissor-")
    ) +
    ggtitle("Scissor Selection on UMAP") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    )
  )
  dev.off()
} else {
  message("UMAP reduction not found, skipping UMAP plots / 未检测到UMAP结果，跳过UMAP绘图")
}

if (nrow(scissor_celltype_summary) > 0) {
  pdf(file = file.path(output_dirs$scissor_analysis, "10_Scissor_CellType_Barplot.pdf"),
      width = 10, height = 7)
  print(
    ggplot(scissor_celltype_summary, aes(x = Cell_Type, y = Cell_Count, color = Scissor)) +
      geom_linerange(aes(ymin = 0, ymax = Cell_Count),                       # lollipop(顶刊优于条形)
                     position = position_dodge(width = 0.6), linewidth = 0.9) +
      geom_point(position = position_dodge(width = 0.6), size = 2.6) +
      coord_flip() +
      scale_color_manual(values = c("Scissor+" = "#E64B35", "Scissor-" = "#4DBBD5")) +
      labs(
        title = "Cell Type Distribution of Scissor-Selected Cells",
        x = "Cell Type",
        y = "Cell Count"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
  dev.off()
}

save(
  Peripheral_Blood_Mononuclear_Cells,
  scissor_result,
  phenotype,
  bulk_group_info,
  cell_mapping,
  file = file.path(output_dirs$final_results, "01_AnnotatedCell_Scissor_Analysis_Complete.RData")
)

write.csv(
  cell_mapping,
  file = file.path(output_dirs$final_results, "02_Cell_Annotation_Scissor_Mapping.csv"),
  row.names = FALSE
)

write.csv(
  scissor_summary,
  file = file.path(output_dirs$final_results, "03_Scissor_Analysis_Summary.csv"),
  row.names = FALSE
)

message("\n================================================================================")
message("                  Annotated Single-Cell Scissor Analysis Complete!              ")
message("                      已注释单细胞对象的Scissor分析完成！                        ")
message("================================================================================")
message(sprintf("Scissor+ cells (disease-associated): %d", length(scissor_pos_cells)))
message(sprintf("Scissor+ 细胞（疾病相关）: %d", length(scissor_pos_cells)))
message(sprintf("Scissor- cells (control-associated): %d", length(scissor_neg_cells)))
message(sprintf("Scissor- 细胞（对照相关）: %d", length(scissor_neg_cells)))
message(sprintf("Selected cell percentage: %.3f%%", selected_ratio))
message(sprintf("筛选细胞占比: %.3f%%", selected_ratio))
message(sprintf("Scissor results saved to: %s", file.path(workDir, output_dirs$scissor_analysis)))
message(sprintf("Scissor结果保存至: %s", file.path(workDir, output_dirs$scissor_analysis)))
