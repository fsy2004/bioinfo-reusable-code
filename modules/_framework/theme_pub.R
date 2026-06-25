# =============================================================================
# theme_pub.R  ·  顶刊级绘图共享主题库 (Top-journal figure toolkit, R)
# -----------------------------------------------------------------------------
# 用途   : 全库统一的 ggplot2 主题、期刊配色、矢量导出与多panel合成工具。
#          每个模块 source() 本文件即可获得一致的 Nature/Cell 级图风格。
# 提供   : theme_pub()  统一主题(Helvetica/Arial、去网格、黑轴线、期刊字号)
#          pal_pub()    期刊离散配色 (NPG/AAAS/Lancet/NEJM/JAMA/Set 等)
#          scale_*_pub  离散/连续(viridis)填充与着色快捷封装
#          save_fig()   一次性导出矢量 PDF + 300dpi PNG(README 预览图)
#          compose_panels()  patchwork 多panel合成并加 A/B/C 角标
#          read_table_smart()  稳健读表 (自动分隔符/编码/行名)
# 依赖   : ggplot2 (必需); ggsci scales patchwork systemfonts ggrepel viridisLite
#          (可选,缺失时自动降级,不报错)
# 约定   : 图中文字一律英文(投稿规范),代码注释中文。
# =============================================================================

suppressWarnings(suppressMessages({
  library(ggplot2)
  .has <- function(p) requireNamespace(p, quietly = TRUE)
  for (.p in c("ggsci", "scales", "patchwork", "systemfonts", "ggrepel",
               "viridisLite", "grid")) if (.has(.p)) suppressMessages(library(.p, character.only = TRUE))
}))

# ---- 字体:优先 Arial/Helvetica(期刊标准无衬线),不可用则降级 sans ----
.pick_font <- function() {
  pref <- c("Arial", "Helvetica", "Helvetica Neue", "Liberation Sans", "DejaVu Sans")
  if (.has("systemfonts")) {
    fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
    hit <- pref[pref %in% fams]
    if (length(hit)) return(hit[1])
  }
  ""  # "" = 设备默认 sans,任何平台都安全
}
PUB_FONT <- .pick_font()

# ---- 期刊离散配色板 ----------------------------------------------------------
# 自带一套不依赖 ggsci 的高级配色,保证缺包也能用;ggsci 在则提供更多选择。
PUB_PALETTES <- list(
  npg    = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4",
             "#91D1C2","#DC0000","#7E6148","#B09C85"),
  aaas   = c("#3B4992","#EE0000","#008B45","#631879","#008280","#BB0021",
             "#5F559B","#A20056","#808180","#1B1919"),
  lancet = c("#00468B","#ED0000","#42B540","#0099B4","#925E9F","#FDAF91",
             "#AD002A","#ADB6B6","#1B1919"),
  nejm   = c("#BC3C29","#0072B5","#E18727","#20854E","#7876B1","#6F99AD",
             "#FFDC91","#EE4C97"),
  jama   = c("#374E55","#DF8F44","#00A1D5","#B24745","#79AF97","#6A6599","#80796B"),
  vivid  = c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
             "#E377C2","#7F7F7F","#BCBD22","#17BECF")
)

#' 取 n 个期刊配色
#' @param n 需要的颜色数 (超出板长自动插值扩展)
#' @param name 配色板名: npg/aaas/lancet/nejm/jama/vivid
pal_pub <- function(n = NULL, name = "npg") {
  base <- PUB_PALETTES[[name]]
  if (is.null(base)) base <- PUB_PALETTES$npg
  if (is.null(n)) return(base)
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)  # 类别过多时平滑扩展
}

# ---- 统一主题 ----------------------------------------------------------------
#' Nature/Cell 风格 ggplot 主题
#' @param base_size 基础字号 (论文单栏建议 8-11)
#' @param grid 是否保留极淡网格线 (默认 FALSE,更干净)
#' @param border 是否加面板外框 (多panel合成时常用)
theme_pub <- function(base_size = 11, base_family = PUB_FONT,
                      grid = FALSE, border = FALSE, legend = "right") {
  half <- base_size / 2
  th <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(size = base_size + 1, face = "bold", hjust = 0,
                                     margin = ggplot2::margin(b = half)),
      plot.subtitle   = element_text(size = base_size - 1, colour = "grey30",
                                     margin = ggplot2::margin(b = half)),
      plot.caption    = element_text(size = base_size - 3, colour = "grey45"),
      axis.title      = element_text(size = base_size, colour = "black"),
      axis.text       = element_text(size = base_size - 1, colour = "black"),
      axis.ticks      = element_line(colour = "black", linewidth = 0.4),
      axis.line       = element_line(colour = "black", linewidth = 0.5),
      panel.border    = if (border) element_rect(colour = "black", fill = NA, linewidth = 0.6)
                        else element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.title    = element_text(size = base_size - 1, face = "bold"),
      legend.text     = element_text(size = base_size - 2),
      legend.key      = element_blank(),
      legend.background = element_blank(),
      legend.position = legend,
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text      = element_text(size = base_size - 1, face = "bold", colour = "black",
                                     margin = ggplot2::margin(half/2, half/2, half/2, half/2)),
      plot.margin     = ggplot2::margin(half, half, half, half)
    )
  if (grid) th <- th + theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                             panel.grid.minor = element_blank())
  else      th <- th + theme(panel.grid = element_blank())
  th
}

# ---- 配色 scale 快捷封装 -----------------------------------------------------
scale_fill_pub  <- function(name = "npg", ...) ggplot2::scale_fill_manual(values = pal_pub(name = name), ...)
scale_color_pub <- function(name = "npg", ...) ggplot2::scale_colour_manual(values = pal_pub(name = name), ...)
scale_colour_pub <- scale_color_pub
# 连续量统一 viridis(色盲友好、印刷稳健)
scale_fill_cont  <- function(option = "D", ...) ggplot2::scale_fill_viridis_c(option = option, ...)
scale_color_cont <- function(option = "D", ...) ggplot2::scale_colour_viridis_c(option = option, ...)

# ---- 一次导出 矢量PDF + 300dpi PNG ------------------------------------------
#' @param plot ggplot/grob 对象
#' @param file 输出路径前缀(不含扩展名)或含扩展名;两种格式都会生成
#' @param width,height 英寸
#' @param dpi PNG 分辨率
save_fig <- function(plot, file, width = 7, height = 6, dpi = 300) {
  stem <- sub("\\.(pdf|png)$", "", file)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  # 矢量 PDF(Cairo 优先,支持系统字体)
  ok_pdf <- tryCatch({
    ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
                    device = grDevices::cairo_pdf); TRUE
  }, error = function(e) tryCatch({
    ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height); TRUE
  }, error = function(e2) FALSE))
  # 300dpi PNG(README 预览;ragg/默认设备即高质量抗锯齿)
  ok_png <- tryCatch({
    ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi); TRUE
  }, error = function(e) FALSE)
  invisible(c(pdf = ok_pdf, png = ok_png))
}

# ---- 多panel合成(patchwork)+ A/B/C 角标 -----------------------------------
#' @param plots ggplot 列表
#' @param ncol,nrow 布局
#' @param tag "A" 大写字母 / "1" 数字 / NULL 不加角标
compose_panels <- function(plots, ncol = NULL, nrow = NULL, tag = "A",
                           widths = NULL, heights = NULL, guides = "keep") {
  if (!.has("patchwork")) { warning("patchwork 缺失,返回首个panel"); return(plots[[1]]) }
  p <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow,
                             widths = widths, heights = heights, guides = guides)
  if (!is.null(tag)) p <- p + patchwork::plot_annotation(
    tag_levels = tag,
    theme = theme(plot.tag = element_text(size = 16, face = "bold", family = PUB_FONT)))
  p
}

# ---- 稳健读表(自动分隔符 / 编码 / 可选首列为行名)---------------------------
read_table_smart <- function(path, row_names = FALSE) {
  if (!file.exists(path)) stop(sprintf("输入文件不存在: %s", path))
  sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
  df <- tryCatch(
    utils::read.csv(path, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, fileEncoding = "UTF-8"),
    error = function(e) utils::read.csv(path, sep = sep, header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE))
  if (isTRUE(row_names) && ncol(df) >= 2) { rn <- df[[1]]; df <- df[, -1, drop = FALSE]; rownames(df) <- make.unique(as.character(rn)) }
  df
}

# ---- 极简 CLI 参数解析: --key value -----------------------------------------
#' 返回 named list;支持 --input x --outdir y --flag(无值=TRUE)
bio_args <- function(defaults = list()) {
  a <- commandArgs(trailingOnly = TRUE); out <- defaults; i <- 1
  while (i <= length(a)) {
    if (grepl("^--", a[i])) {
      k <- sub("^--", "", a[i])
      if (i + 1 <= length(a) && !grepl("^--", a[i + 1])) { out[[k]] <- a[i + 1]; i <- i + 2 }
      else { out[[k]] <- TRUE; i <- i + 1 }
    } else i <- i + 1
  }
  out
}

# ---- 定位当前脚本所在目录(Rscript 下可靠)-----------------------------------
bio_script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", a[m[1]]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}

# ---- 零依赖 Venn 图(2-3 集合;更多请用 UpSet)-------------------------------
#' @param sets 命名 list,每元素为一组元素向量(如基因名)
#' @param fill 各集合填充色(默认期刊配色)
venn_pub <- function(sets, fill = NULL, title = NULL, base_size = 13) {
  n <- length(sets); nm <- names(sets); if (is.null(nm)) nm <- paste0("Set", seq_len(n))
  if (!n %in% 2:3) stop("venn_pub 仅支持 2-3 集合;更多请用 UpSetR。")
  if (is.null(fill)) fill <- pal_pub(n, "npg")
  circle <- function(x, y, r, k) { t <- seq(0, 2 * pi, length.out = 200); data.frame(x = x + r * cos(t), y = y + r * sin(t), grp = k) }
  cnt <- function(inc, exc) length(setdiff(Reduce(intersect, sets[inc]), if (length(exc)) Reduce(union, sets[exc]) else character(0)))
  if (n == 2) {
    cx <- c(-0.6, 0.6); cy <- c(0, 0); r <- 1.15
    polys <- rbind(circle(cx[1], cy[1], r, 1), circle(cx[2], cy[2], r, 2))
    lab <- data.frame(x = c(-1.1, 1.1, 0), y = c(0, 0, 0),
                      l = c(cnt(1, 2), cnt(2, 1), cnt(c(1, 2), integer(0))))
    setlab <- data.frame(x = cx, y = c(1.5, 1.5), l = nm)
  } else {
    cx <- c(-0.6, 0.6, 0); cy <- c(0.5, 0.5, -0.6); r <- 1.25
    polys <- rbind(circle(cx[1], cy[1], r, 1), circle(cx[2], cy[2], r, 2), circle(cx[3], cy[3], r, 3))
    lab <- data.frame(
      x = c(-1.2, 1.2, 0, 0, -0.85, 0.85, 0),
      y = c(1.1, 1.1, -1.2, 0.85, -0.35, -0.35, 0.1),
      l = c(cnt(1, c(2, 3)), cnt(2, c(1, 3)), cnt(3, c(1, 2)),
            cnt(c(1, 2), 3), cnt(c(1, 3), 2), cnt(c(2, 3), 1), cnt(c(1, 2, 3), integer(0))))
    setlab <- data.frame(x = c(-1.3, 1.3, 0), y = c(1.7, 1.7, -1.8), l = nm)
  }
  ggplot() +
    geom_polygon(data = polys, aes(x, y, group = grp, fill = factor(grp)), alpha = 0.45, colour = "grey30") +
    geom_text(data = lab, aes(x, y, label = l), fontface = "bold", size = base_size / 3) +
    geom_text(data = setlab, aes(x, y, label = l), fontface = "bold", size = base_size / 2.6, colour = fill) +
    scale_fill_manual(values = fill, guide = "none") +
    coord_equal() + labs(title = title) +
    theme_void(base_family = PUB_FONT) +
    theme(plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5))
}

invisible(message("[theme_pub] 顶刊主题已加载 · 字体=", ifelse(PUB_FONT == "", "sans(default)", PUB_FONT)))
