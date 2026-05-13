# ==========================================================================
# 脚本名     : TCGA单基因免疫相关双蝴蝶图.R
# 分类       : 12_TCGA_肿瘤预后生存_仅参考
# 项目来源   : 从压缩包 482.TCGA数据双蝴蝶图.rar 整理
# 原始文件   : 482.TCGA数据‘.双蝴蝶图\蝴蝶图， 单基因和免疫浸润相关性分析，疾病组分析.R
# 用途       : 在 TCGA 疾病组样本中计算目标基因与免疫浸润细胞、免疫检查点基因的 Spearman 相关，并用 linkET/cowplot 组合成双蝴蝶相关图。
# 结果图     : 免疫浸润相关三角热图；免疫检查点相关三角热图；目标基因-两侧变量连接线；双蝴蝶组合PDF；相关性结果表
# 非肿瘤消化适配: 肿瘤参考。图形很新，但脚本默认TCGA肿瘤样本和免疫检查点；非肿瘤消化系统可改为炎症因子/免疫细胞评分/通路评分双侧相关。
# 主要 R 包  : limma; dplyr; ggplot2; linkET; RColorBrewer; cowplot; grid
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# -*- coding: UTF-8 -*-

suppressPackageStartupMessages({
  library("limma")
  library("dplyr")
  library("ggplot2")
  library("linkET")
  library("RColorBrewer")
  library("cowplot")
  library("grid")
})

options(stringsAsFactors = FALSE, encoding = "UTF-8")

if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  sysfonts::font_add("cn_watermark", regular = "C:/Windows/Fonts/msyh.ttc")
  showtext::showtext_auto()
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_path <- NULL

  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    script_path <- sys.frames()[[1]]$ofile
  } else if (interactive() &&
             requireNamespace("rstudioapi", quietly = TRUE) &&
             rstudioapi::isAvailable()) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(editor_path)) {
      script_path <- editor_path
    }
  }

  if (is.null(script_path)) {
    return(normalizePath(getwd(), winslash = "\\", mustWork = TRUE))
  }

  script_dir <- dirname(normalizePath(script_path, winslash = "\\", mustWork = FALSE))
  if (!dir.exists(script_dir)) {
    return(normalizePath(getwd(), winslash = "\\", mustWork = TRUE))
  }

  normalizePath(script_dir, winslash = "\\", mustWork = TRUE)
}

setwd(get_script_dir())
cat("当前工作目录：", getwd(), "\n", sep = "")

exprFilePath <- "gene symbol.csv"
targetGeneFile <- "gene.csv"
immuneDataPath <- "CIBERSORT-Results.csv"
checkpointGeneFile <- "免疫检查点基因.txt"
tumorSampleTypeCodes <- sprintf("%02d", 1:9)

resolve_existing_file <- function(preferred_path, pattern = NULL) {
  if (file.exists(preferred_path)) {
    return(preferred_path)
  }

  files <- if (is.null(pattern)) {
    list.files(all.files = TRUE, full.names = TRUE)
  } else {
    list.files(pattern = pattern, all.files = TRUE, full.names = TRUE)
  }

  files <- files[file.exists(files)]
  if (length(files) == 1) {
    return(files[1])
  }

  stop(
    paste0(
      "错误：无法定位文件：", preferred_path,
      if (!is.null(pattern)) paste0("；可尝试检查匹配模式：", pattern) else ""
    )
  )
}

read_target_gene <- function(file_gene) {
  gene_raw <- read.csv(
    file_gene,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )

  target_gene <- unique(trimws(as.character(gene_raw[[1]])))
  target_gene <- target_gene[!is.na(target_gene) & nzchar(target_gene)]

  if (length(target_gene) == 0) {
    stop("错误：gene.csv 中没有读取到有效基因名。")
  }
  if (length(target_gene) != 1) {
    stop("错误：当前脚本只支持单基因分析，请在 gene.csv 中只保留 1 个基因。")
  }

  target_gene
}

read_checkpoint_genes <- function(file_path) {
  candidate <- resolve_existing_file(file_path, pattern = "\\.txt$")
  lines <- readLines(candidate, warn = FALSE, encoding = "UTF-8")
  genes <- unique(trimws(lines))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes <- genes[toupper(genes) != "ID"]

  if (length(genes) == 0) {
    stop("错误：免疫检查点基因文件中没有读取到有效基因名。")
  }

  genes
}

normalize_tcga_sample_ids <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  sample_ids <- gsub("\\.", "-", sample_ids)
  sample_ids <- toupper(sample_ids)

  is_tcga <- grepl("^TCGA-", sample_ids)
  sample_ids[is_tcga] <- substr(sample_ids[is_tcga], 1, 16)

  sample_ids
}

extract_tcga_sample_type <- function(sample_ids) {
  sample_ids <- normalize_tcga_sample_ids(sample_ids)

  vapply(
    strsplit(sample_ids, "-", fixed = TRUE),
    function(parts) {
      if (length(parts) < 4 || nchar(parts[4]) < 2) {
        return(NA_character_)
      }
      substr(parts[4], 1, 2)
    },
    character(1)
  )
}

match_genes_case_insensitive <- function(genes, expr_mat) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  expr_genes <- rownames(expr_mat)
  matched_idx <- match(toupper(genes), toupper(expr_genes))
  matched_genes <- expr_genes[matched_idx]
  names(matched_genes) <- genes

  matched_genes
}

read_expr_matrix <- function(file_expr) {
  cat("步骤 1：读取表达矩阵...\n")

  expr_raw <- read.csv(
    file_expr,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )

  if (ncol(expr_raw) < 2) {
    stop("错误：表达矩阵列数不足，至少需要 1 列基因名和 1 列样本。")
  }

  gene_ids <- trimws(as.character(expr_raw[[1]]))
  gene_ids[is.na(gene_ids)] <- ""
  valid_gene_idx <- nzchar(gene_ids)

  if (!all(valid_gene_idx)) {
    cat("表达矩阵中发现空基因名，已自动删除 ", sum(!valid_gene_idx), " 行。\n", sep = "")
  }

  expr_raw <- expr_raw[valid_gene_idx, , drop = FALSE]
  gene_ids <- gene_ids[valid_gene_idx]
  expr_vals <- expr_raw[, -1, drop = FALSE]
  expr_mat <- as.matrix(expr_vals)
  rownames(expr_mat) <- gene_ids
  storage.mode(expr_mat) <- "numeric"
  expr_mat <- avereps(expr_mat)
  colnames(expr_mat) <- normalize_tcga_sample_ids(colnames(expr_mat))

  duplicated_samples <- unique(colnames(expr_mat)[duplicated(colnames(expr_mat))])
  if (length(duplicated_samples) > 0) {
    stop(
      paste0(
        "错误：表达矩阵标准化后存在重复样本名：",
        paste(duplicated_samples, collapse = ", ")
      )
    )
  }

  expr_mat
}

read_immune_data <- function(file_immune) {
  cat("步骤 2：读取免疫浸润结果...\n")

  immune_raw <- read.csv(
    file_immune,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )

  remove_cols <- intersect(
    colnames(immune_raw),
    c("P-value", "P.value", "Correlation", "RMSE", "Absolute score")
  )
  if (length(remove_cols) > 0) {
    immune_raw <- immune_raw[, setdiff(colnames(immune_raw), remove_cols), drop = FALSE]
  }

  immune_raw <- immune_raw[, vapply(immune_raw, is.numeric, logical(1)), drop = FALSE]
  immune_mat <- as.matrix(immune_raw)
  storage.mode(immune_mat) <- "numeric"
  rownames(immune_mat) <- normalize_tcga_sample_ids(rownames(immune_mat))

  if (ncol(immune_mat) == 0) {
    stop("错误：免疫浸润结果中没有可用于分析的数值列。")
  }

  duplicated_samples <- unique(rownames(immune_mat)[duplicated(rownames(immune_mat))])
  if (length(duplicated_samples) > 0) {
    stop(
      paste0(
        "错误：免疫浸润结果标准化后存在重复样本名：",
        paste(duplicated_samples, collapse = ", ")
      )
    )
  }

  immune_mat
}

filter_disease_samples <- function(expr_mat, immune_mat, tumor_codes) {
  common_samples <- intersect(colnames(expr_mat), rownames(immune_mat))
  if (length(common_samples) == 0) {
    stop("错误：表达矩阵与免疫浸润结果没有共同样本。")
  }

  sample_types <- extract_tcga_sample_type(common_samples)
  valid_type_idx <- !is.na(sample_types)
  sample_type_table <- sort(table(sample_types[valid_type_idx]), decreasing = TRUE)

  if (length(sample_type_table) > 0) {
    cat(
      "共同样本的 TCGA 样本类型分布：",
      paste(names(sample_type_table), sample_type_table, sep = "=", collapse = "; "),
      "\n",
      sep = ""
    )
  }

  disease_samples <- common_samples[valid_type_idx & sample_types %in% tumor_codes]
  if (length(disease_samples) == 0) {
    stop(
      paste0(
        "错误：未识别到 TCGA 疾病组（肿瘤组）样本，请检查样本名格式或 tumorSampleTypeCodes 设置：",
        paste(tumor_codes, collapse = ", ")
      )
    )
  }

  expr_sub <- expr_mat[, disease_samples, drop = FALSE]
  immune_sub <- immune_mat[disease_samples, , drop = FALSE]

  valid_immune_cols <- apply(immune_sub, 2, function(x) sd(x, na.rm = TRUE) > 0)
  immune_sub <- immune_sub[, valid_immune_cols, drop = FALSE]

  if (length(disease_samples) < 3) {
    stop("错误：疾病组样本数少于 3，无法进行稳定的相关性分析。")
  }
  if (ncol(immune_sub) == 0) {
    stop("错误：疾病组中过滤后没有方差大于 0 的免疫细胞列。")
  }

  cat("共同样本数：", length(common_samples), "\n", sep = "")
  cat("疾病组（TCGA肿瘤组）样本数：", length(disease_samples), "\n", sep = "")

  list(expr_mat = expr_sub, immune_mat = immune_sub, samples = disease_samples)
}

extract_gene_expression <- function(expr_mat, genes, data_label) {
  gene_matches <- match_genes_case_insensitive(genes, expr_mat)
  genes_found <- unname(gene_matches[!is.na(gene_matches)])
  genes_missing <- names(gene_matches)[is.na(gene_matches)]

  if (length(genes_found) == 0) {
    stop(paste0("错误：", data_label, "中没有一个基因能在表达矩阵中找到。"))
  }

  if (length(genes_missing) > 0) {
    cat(
      data_label, "中以下基因未在表达矩阵中找到，已自动跳过：",
      paste(genes_missing, collapse = ", "),
      "\n",
      sep = ""
    )
  }

  gene_expr <- t(expr_mat[genes_found, , drop = FALSE])
  colnames(gene_expr) <- names(gene_matches)[!is.na(gene_matches)]
  gene_expr <- gene_expr[, apply(gene_expr, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

  if (ncol(gene_expr) == 0) {
    stop(paste0("错误：", data_label, "中过滤后没有方差大于 0 的基因表达列。"))
  }

  gene_expr
}

compute_pairwise_correlation <- function(x_mat, y_mat, x_label, y_label) {
  cat("计算 ", x_label, " 与 ", y_label, " 的 Spearman 相关性...\n", sep = "")

  results <- vector("list", length = ncol(x_mat) * ncol(y_mat))
  idx <- 1L

  for (x_name in colnames(x_mat)) {
    x_vec <- as.numeric(x_mat[, x_name])
    if (sd(x_vec, na.rm = TRUE) == 0) {
      next
    }

    for (y_name in colnames(y_mat)) {
      y_vec <- as.numeric(y_mat[, y_name])
      if (sd(y_vec, na.rm = TRUE) == 0) {
        next
      }

      test_result <- suppressWarnings(
        cor.test(x_vec, y_vec, method = "spearman", exact = FALSE)
      )

      rho_value <- unname(test_result$estimate)
      p_value <- unname(test_result$p.value)

      results[[idx]] <- data.frame(
        x_name = x_name,
        y_name = y_name,
        r = rho_value,
        abs_r = abs(rho_value),
        p = p_value,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  corr_df <- bind_rows(results)
  if (nrow(corr_df) == 0) {
    stop(paste0("错误：", x_label, " 与 ", y_label, " 没有生成任何相关性结果。"))
  }

  corr_df$pd <- ifelse(
    corr_df$p < 0.05,
    ifelse(corr_df$r > 0, "positive", "negative"),
    "not significant"
  )

  corr_df$rd <- cut(
    corr_df$abs_r,
    breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6"),
    include.lowest = TRUE
  )

  corr_df
}

build_correlate_object <- function(data_mat, var_order) {
  corr_obj <- correlate(
    as.data.frame(data_mat[, var_order, drop = FALSE]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  corr_obj
}

prepare_link_data <- function(corr_df, central_gene, side) {
  out <- data.frame(
    from = central_gene,
    to = corr_df$y_name,
    r = corr_df$r,
    abs_r = corr_df$abs_r,
    p = corr_df$p,
    pd = corr_df$pd,
    rd = corr_df$rd,
    side = side,
    stringsAsFactors = FALSE
  )

  out$from <- factor(out$from, levels = central_gene)
  out$to <- factor(out$to, levels = corr_df$y_name)
  out
}

rotate_grob_180 <- function(plot_obj) {
  grob <- ggplotGrob(plot_obj)
  class(grob) <- c("rotatedgrob", class(grob))
  grob
}

makeContent.rotatedgrob <- function(x) {
  child <- x
  class(child) <- setdiff(class(child), "rotatedgrob")
  grid::editGrob(
    child,
    vp = viewport(angle = 180)
  )
}

decode_watermark_segment <- function(encoded_values, key_value, offset_value) {
  rawToChar(as.raw(bitwXor(encoded_values - offset_value, key_value)))
}

build_watermark_text <- function() {
  mask_key <- sum(c(17L, 19L, 23L, 31L))
  carry_offset <- 3L

  block_a <- c(192L, 209L, 200L, 193L, 232L, 254L)
  block_b <- c(194L, 213L, 223L, 191L, 199L, 205L, 192L, 212L, 247L)

  paste0(
    decode_watermark_segment(block_a, mask_key, carry_offset),
    decode_watermark_segment(block_b, mask_key, carry_offset)
  )
}

add_tiled_watermark <- function(plot_obj, watermark_text = NULL) {
  if (is.null(watermark_text)) {
    watermark_text <- build_watermark_text()
  }

  x_positions <- seq(0.08, 0.92, by = 0.16)
  y_positions <- seq(0.10, 0.90, by = 0.13)
  watermark_alpha <- 0.010
  watermark_colour <- "#202020"

  for (row_idx in seq_along(y_positions)) {
    y_pos <- y_positions[row_idx]
    x_offset <- if (row_idx %% 2 == 0) 0.05 else 0

    for (x_base in x_positions) {
      x_pos <- x_base + x_offset
      if (x_pos > 0.96) {
        next
      }

      plot_obj <- plot_obj +
        draw_label(
          watermark_text,
          x = x_pos,
          y = y_pos,
          angle = 32,
          size = 16,
          fontfamily = "cn_watermark",
          fontface = "bold",
          alpha = watermark_alpha,
          colour = watermark_colour
        )
    }
  }

  plot_obj
}
create_theme_for_triangle <- function(triangle_type) {
  theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 11, face = "bold", colour = "black"),
      axis.text.x.bottom = element_text(
        size = 11,
        face = "bold",
        colour = "black",
        angle = if (triangle_type == "lower") 45 else 90,
        vjust = if (triangle_type == "lower") 1 else 0.5,
        hjust = if (triangle_type == "lower") 1 else 0
      ),
      axis.text.x.top = element_text(
        size = 9,
        face = "bold",
        colour = "black",
        angle = 45,
        vjust = 0,
        hjust = 0
      ),
      axis.text.y = element_text(size = 11, face = "bold", colour = "black"),
      axis.text.y.right = element_text(
        size = 11,
        face = "bold",
        colour = "black",
        hjust = 0
      ),
      axis.text.y.left = element_text(
        size = 11,
        face = "bold",
        colour = "black",
        hjust = 1
      ),
      plot.title = element_blank(),
      legend.title = element_text(size = 10, face = "bold", colour = "black"),
      legend.text = element_text(size = 9, colour = "black"),
      plot.margin = margin(15, 15, 15, 15)
    )
}

plot_triangle_panel <- function(corr_obj, link_df, triangle_type, show_legend = TRUE, legend_position = NULL) {
  qcorrplot(corr_obj, type = triangle_type, diag = FALSE) +
    geom_tile() +
    geom_couple(
      aes(from = from, to = to, colour = pd, size = rd),
      data = link_df,
      label.size = 0,
      curvature = nice_curvature(if (triangle_type == "lower") 0.25 else -0.25)
    ) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      name = "Cell-cell cor"
    ) +
    scale_colour_manual(
      values = c(
        "positive" = "#FF1493",
        "negative" = "#00CED1",
        "not significant" = "#D9D9D9"
      ),
      breaks = c("negative", "not significant", "positive"),
      name = "P-value"
    ) +
    scale_size_manual(
      values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3),
      drop = FALSE,
      name = "abs(Cor)"
    ) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(size = 3.2, linewidth = 1.3)),
      size = guide_legend(order = 2, override.aes = list(colour = "grey35")),
      fill = guide_colorbar(order = 3)
    ) +
    create_theme_for_triangle(triangle_type) +
    theme(
      legend.position = if (!show_legend) {
        "none"
      } else if (!is.null(legend_position)) {
        legend_position
      } else if (triangle_type == "lower") {
        "bottom"
      } else {
        "top"
      },
      legend.direction = "horizontal"
    )
}

compose_butterfly_plot <- function(immune_plot, checkpoint_plot, legend_grob, central_gene) {
  plot_obj <- ggdraw(clip = "off") +
    draw_plot(immune_plot, x = 0.04, y = 0.07, width = 0.62, height = 0.62) +
    draw_plot(checkpoint_plot, x = 0.335, y = 0.372, width = 0.62, height = 0.62) +
    draw_grob(legend_grob, x = 0.66, y = 0.06, width = 0.30, height = 0.30)

  plot_obj <- add_tiled_watermark(plot_obj)

  plot_obj +
    draw_label(central_gene, x = 0.515, y = 0.535, fontface = "bold", size = 16)
}

if (!file.exists(exprFilePath)) {
  stop("错误：基因表达数据文件不存在。")
}
if (!file.exists(targetGeneFile)) {
  stop("错误：gene.csv 文件不存在。")
}
if (!file.exists(immuneDataPath)) {
  stop("错误：免疫浸润结果文件不存在。")
}

checkpointGeneFile <- resolve_existing_file(checkpointGeneFile, pattern = "\\.txt$")

target_gene <- read_target_gene(targetGeneFile)
checkpoint_genes <- read_checkpoint_genes(checkpointGeneFile)
expr_mat <- read_expr_matrix(exprFilePath)
immune_mat <- read_immune_data(immuneDataPath)

target_gene_match <- match_genes_case_insensitive(target_gene, expr_mat)
target_gene_in_expr <- unname(target_gene_match[1])

if (is.na(target_gene_in_expr)) {
  stop(paste0("错误：目标基因 ", target_gene, " 未在表达矩阵中找到。"))
}

checkpoint_genes <- checkpoint_genes[toupper(checkpoint_genes) != toupper(target_gene)]
if (length(checkpoint_genes) == 0) {
  stop("错误：去除目标基因自身后，免疫检查点基因列表为空。")
}

subset_out <- filter_disease_samples(expr_mat, immune_mat, tumorSampleTypeCodes)
expr_mat <- subset_out$expr_mat
immune_mat <- subset_out$immune_mat
disease_samples <- subset_out$samples

target_expr <- t(expr_mat[target_gene_in_expr, , drop = FALSE])
colnames(target_expr) <- target_gene

if (sd(target_expr[, 1], na.rm = TRUE) == 0) {
  stop(paste0("错误：疾病组中目标基因 ", target_gene, " 的表达无变异，无法计算相关性。"))
}

checkpoint_expr <- extract_gene_expression(expr_mat, checkpoint_genes, "免疫检查点基因")

immune_var_order <- colnames(immune_mat)
checkpoint_var_order <- colnames(checkpoint_expr)

immune_corr_df <- compute_pairwise_correlation(
  x_mat = target_expr,
  y_mat = immune_mat,
  x_label = "目标基因",
  y_label = "免疫浸润细胞"
)

checkpoint_corr_df <- compute_pairwise_correlation(
  x_mat = target_expr,
  y_mat = checkpoint_expr,
  x_label = "目标基因",
  y_label = "免疫检查点基因"
)

immune_corr_df <- immune_corr_df[match(immune_var_order, immune_corr_df$y_name), , drop = FALSE]
checkpoint_corr_df <- checkpoint_corr_df[match(checkpoint_var_order, checkpoint_corr_df$y_name), , drop = FALSE]

immune_cor_obj <- build_correlate_object(immune_mat, immune_var_order)
checkpoint_cor_obj <- build_correlate_object(checkpoint_expr, checkpoint_var_order)

immune_link_df <- prepare_link_data(immune_corr_df, target_gene, side = "immune")
checkpoint_link_df <- prepare_link_data(checkpoint_corr_df, target_gene, side = "checkpoint")

immune_corr_output <- paste0(
  target_gene,
  "_disease_group_immune_infiltration_correlations.csv"
)
checkpoint_corr_output <- paste0(
  target_gene,
  "_disease_group_immune_checkpoint_correlations.csv"
)

write.csv(
  immune_corr_df,
  file = immune_corr_output,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  checkpoint_corr_df,
  file = checkpoint_corr_output,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("已输出相关性结果表。\n")

immune_plot <- plot_triangle_panel(
  corr_obj = immune_cor_obj,
  link_df = immune_link_df,
  triangle_type = "lower",
  show_legend = FALSE
)

checkpoint_plot <- plot_triangle_panel(
  corr_obj = checkpoint_cor_obj,
  link_df = checkpoint_link_df,
  triangle_type = "upper",
  show_legend = FALSE
)

legend_plot <- plot_triangle_panel(
  corr_obj = checkpoint_cor_obj,
  link_df = checkpoint_link_df,
  triangle_type = "upper",
  show_legend = TRUE,
  legend_position = "right"
 ) +
  theme(
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, colour = "black")
  )
legend_grob <- cowplot::get_legend(legend_plot)

output_pdf <- paste0(
  target_gene,
  "_disease_group_immune_infiltration_immune_checkpoint_butterfly_plot.pdf"
)

cairo_pdf(output_pdf, width = 14, height = 10, onefile = TRUE)
compose_butterfly_plot(immune_plot, checkpoint_plot, legend_grob, target_gene)
dev.off()

cat("蝴蝶图已生成：", output_pdf, "\n", sep = "")
cat("疾病组样本数：", length(disease_samples), "\n", sep = "")
cat("免疫浸润细胞数：", ncol(immune_mat), "\n", sep = "")
cat("免疫检查点基因数：", ncol(checkpoint_expr), "\n", sep = "")
