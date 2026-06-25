# =============================================================================
# 编号       : 007
# 脚本名     : GO / KEGG 富集分析 (turnkey + 顶刊图)
# 分类       : 02_enrichment
# 用途       : 对一组候选基因做 GO(BP/CC/MF)与 KEGG 通路富集,输出通路表与
#              顶刊级合成图(GO 分面点图 + KEGG 棒棒糖图 + 多panel合成图)。
# 方法/包    : clusterProfiler(enrichGO/enrichKEGG) + org.Hs.eg.db + DOSE;
#              绘图 ggplot2 + 共享主题 theme_pub.R(viridis 期刊配色)。
# 结果图     : GO_enrichment_dotplot;KEGG_enrichment_lollipop;Fig_enrichment_composite
# -----------------------------------------------------------------------------
# 运行(零改动跑示例): Rscript 007_GO_KEGG_enrichment.R
# 运行(自己的数据)  : Rscript 007_GO_KEGG_enrichment.R --input data/genes.csv --outdir results/run1
# 可选参数: --top 10(每类展示数) --pvalue 0.05 --padjust 0.05 --organism hsa --keytype SYMBOL
# 输入规格: 单列 CSV,列名 Gene,每行一个人类基因 Symbol(或 ENTREZ,见 --keytype)。
# 整理日期 : 2026-06-23(turnkey 重构;分析逻辑保持 clusterProfiler 标准流程)
# =============================================================================

# ---- turnkey preamble: 定位共享框架并加载顶刊主题 ----------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R,请确认模块位于库目录结构内。")
}
source(.find_fw())

suppressWarnings(suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(dplyr); library(ggplot2)
}))

# ---- 参数(默认指向 example_data,可被 --key value 覆盖)----------------------
SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(
  input    = file.path(SCRIPT_DIR, "example_data", "gene_list.csv"),
  outdir   = file.path(SCRIPT_DIR, "results"),
  top      = "8",      # 每个本体/通路类别展示条目数
  pvalue   = "0.05",
  padjust  = "0.05",
  organism = "hsa",    # KEGG organism
  keytype  = "SYMBOL"  # 输入基因类型: SYMBOL 或 ENTREZID
))
TOP <- as.integer(args$top); PCUT <- as.numeric(args$pvalue); QCUT <- as.numeric(args$padjust)
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1. 读基因并转 ENTREZ ----------------------------------------------
cat("Step 1/5: 读取基因列表...\n")
genes_raw <- read_table_smart(args$input)
gcol <- if ("Gene" %in% names(genes_raw)) "Gene" else names(genes_raw)[1]
gene_symbols <- unique(na.omit(as.character(genes_raw[[gcol]])))
if (length(gene_symbols) == 0) stop("未读到任何基因,请检查输入文件(需列名 Gene)。")
cat("  读入基因数:", length(gene_symbols), "\n")

if (toupper(args$keytype) == "SYMBOL") {
  conv <- suppressWarnings(bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  entrez <- unique(conv$ENTREZID)
  cat("  成功转换 ENTREZ:", length(entrez), "/", length(gene_symbols), "\n")
} else { entrez <- unique(gene_symbols) }
if (length(entrez) < 3) stop("有效 ENTREZ 基因过少,无法富集。")

# ---- Step 2. GO 富集(离线,org.Hs.eg.db)-----------------------------------
cat("Step 2/5: GO 富集 (BP/CC/MF)...\n")
ego <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "ALL",
                pvalueCutoff = PCUT, qvalueCutoff = QCUT, readable = TRUE)
go_df <- if (!is.null(ego)) as.data.frame(ego) else data.frame()
if (nrow(go_df)) write.csv(go_df, file.path(args$outdir, "GO_results.csv"), row.names = FALSE)
cat("  显著 GO 条目:", nrow(go_df), "\n")

# ---- Step 3. KEGG 富集(需联网;失败则跳过,不中断)--------------------------
cat("Step 3/5: KEGG 富集...\n")
kegg_df <- data.frame()
ekegg <- tryCatch(
  enrichKEGG(gene = entrez, organism = args$organism, pvalueCutoff = PCUT, qvalueCutoff = QCUT),
  error = function(e) { cat("  KEGG 在线查询失败(可能无网络),跳过:", conditionMessage(e), "\n"); NULL })
if (!is.null(ekegg)) {
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(ekegg)
  if (nrow(kegg_df)) write.csv(kegg_df, file.path(args$outdir, "KEGG_results.csv"), row.names = FALSE)
}
cat("  显著 KEGG 通路:", nrow(kegg_df), "\n")

# 通用: GeneRatio "a/b" → 数值
ratio_num <- function(x) vapply(strsplit(as.character(x), "/"), function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))

# 基因–通路概念网络(cnetplot 风,ggraph 手绘可主题化)
build_concept_network <- function(df, top_n = 5, title = "Gene–pathway network") {
  if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("tidygraph", quietly = TRUE) ||
      is.null(df) || nrow(df) == 0 || !"geneID" %in% names(df)) return(NULL)
  top <- df %>% slice_min(p.adjust, n = top_n, with_ties = FALSE)
  edges <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    gs <- strsplit(top$geneID[i], "/")[[1]]
    data.frame(from = top$Description[i], to = gs, stringsAsFactors = FALSE)
  }))
  if (is.null(edges) || nrow(edges) == 0) return(NULL)
  path_nodes <- data.frame(name = top$Description, type = "Pathway",
                           size = top$Count, col = pal_pub(nrow(top), "npg"), stringsAsFactors = FALSE)
  gene_nodes <- data.frame(name = unique(edges$to), type = "Gene", size = 1, col = "grey60", stringsAsFactors = FALSE)
  nodes <- rbind(path_nodes, gene_nodes); nodes <- nodes[!duplicated(nodes$name), ]
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  set.seed(42)
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(colour = "grey80", alpha = 0.5, linewidth = 0.35) +
    ggraph::geom_node_point(aes(size = size, colour = col, filter = type == "Pathway")) +
    ggraph::geom_node_point(aes(filter = type == "Gene"), size = 1.8, colour = "grey60") +
    ggraph::geom_node_text(aes(label = ifelse(type == "Gene", name, "")),
                           repel = TRUE, size = 2.5, colour = "grey25", max.overlaps = 80,
                           segment.colour = "grey85", segment.size = 0.2) +
    ggraph::geom_node_label(aes(label = ifelse(type == "Pathway", name, ""), colour = col),
                            repel = TRUE, size = 3.2, fontface = "bold", label.size = 0.4,
                            fill = "white", alpha = 0.9, max.overlaps = 50, show.legend = FALSE) +
    scale_colour_identity() +
    scale_size_continuous(range = c(6, 14), name = "Gene count") +
    labs(title = title) +
    ggraph::theme_graph(base_family = PUB_FONT) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0),
          legend.position = "right")
}

# ---- Step 4. 顶刊级图 --------------------------------------------------------
cat("Step 4/5: 绘制顶刊级图...\n")

## (A) GO 分面点图:BP/CC/MF 各取 top,点大小=Count,颜色=-log10(p.adjust),viridis
p_go <- NULL
if (nrow(go_df)) {
  go_top <- go_df %>% group_by(ONTOLOGY) %>% slice_min(p.adjust, n = TOP, with_ties = FALSE) %>% ungroup() %>%
    mutate(GeneRatio = ratio_num(GeneRatio),
           Description = factor(Description, levels = rev(unique(Description[order(ONTOLOGY, p.adjust)]))),
           neglog10P = -log10(p.adjust))
  p_go <- ggplot(go_top, aes(GeneRatio, Description)) +
    geom_segment(aes(x = 0, xend = GeneRatio, yend = Description), colour = "grey80", linewidth = 0.4) +
    geom_point(aes(size = Count, colour = neglog10P)) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_colour_viridis_c(option = "D", name = expression(-log[10]~italic(P)[adj])) +
    scale_size_continuous(range = c(2.5, 7), name = "Gene count") +
    labs(x = "Gene ratio", y = NULL, title = "GO enrichment (BP / CC / MF)") +
    theme_pub(base_size = 11, border = TRUE) +
    theme(axis.text.y = element_text(size = 9))
  save_fig(p_go, file.path(ASSETS, "GO_enrichment_dotplot"), width = 8.2, height = 8)
  save_fig(p_go, file.path(args$outdir, "GO_enrichment_dotplot"), width = 8.2, height = 8)
}

## (B) KEGG 棒棒糖图:top 通路,颜色梯度=-log10(p.adjust)
p_kegg <- NULL
if (nrow(kegg_df)) {
  kg_top <- kegg_df %>% slice_min(p.adjust, n = max(TOP, 12), with_ties = FALSE) %>%
    mutate(Description = factor(Description, levels = rev(Description)), neglog10P = -log10(p.adjust))
  p_kegg <- ggplot(kg_top, aes(Count, Description)) +
    geom_segment(aes(x = 0, xend = Count, yend = Description), colour = "grey80", linewidth = 0.5) +
    geom_point(aes(size = Count, colour = neglog10P)) +
    scale_colour_viridis_c(option = "C", name = expression(-log[10]~italic(P)[adj])) +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    labs(x = "Gene count", y = NULL, title = "KEGG pathway enrichment") +
    theme_pub(base_size = 11, border = TRUE) +
    theme(axis.text.y = element_text(size = 9))
  save_fig(p_kegg, file.path(ASSETS, "KEGG_enrichment_lollipop"), width = 7.5, height = 6.5)
  save_fig(p_kegg, file.path(args$outdir, "KEGG_enrichment_lollipop"), width = 7.5, height = 6.5)
}

## (C) 基因–通路概念网络(优先 KEGG,无则用 GO-BP)
net_src <- if (nrow(kegg_df)) kegg_df else go_df[go_df$ONTOLOGY == "BP", , drop = FALSE]
net_title <- if (nrow(kegg_df)) "Gene–KEGG pathway network (top 5)" else "Gene–GO:BP network (top 5)"
p_net <- tryCatch(build_concept_network(net_src, top_n = 5, title = net_title),
                  error = function(e) { cat("  网络图绘制跳过:", conditionMessage(e), "\n"); NULL })
if (!is.null(p_net)) {
  save_fig(p_net, file.path(ASSETS, "GenePathway_network"), width = 9, height = 7.5)
  save_fig(p_net, file.path(args$outdir, "GenePathway_network"), width = 9, height = 7.5)
}

# 说明:按规范只输出独立单图(GO 点图 / KEGG 棒棒糖 / 概念网络),投稿拼版自理;
#       如需自动合成,可调用 compose_panels()(theme_pub.R 提供),此处默认不生成。

# ---- Step 5. 完成 ------------------------------------------------------------
cat("Step 5/5: 完成。结果表与图见:", normalizePath(args$outdir), "\n")
cat("  展示图(README 用)见:", normalizePath(ASSETS), "\n")
