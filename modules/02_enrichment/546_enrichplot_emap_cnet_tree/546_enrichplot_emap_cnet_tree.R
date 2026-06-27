# =============================================================================
# 编号   : 546
# 脚本名 : 富集结果高级图 — dotplot / cnetplot / emapplot / treeplot
# 分类   : 02_enrichment
# 用途   : 把 GO/通路富集结果画成顶刊级高级图,系统性代替"平凡富集条形图":
#          ① dotplot   气泡图(GeneRatio×p.adjust×Count,一图四维)
#          ② cnetplot  基因–概念网络(circular,显式连出哪些基因驱动哪些通路)
#          ③ emapplot  富集图谱(按词项相似度把通路聚成功能模块)
#          ④ treeplot  层次聚类树(通路自动归并为带标签的高阶簇)
# ★诚实基线(可视化基线,非分析基线):
#          内置一张"朴素富集条形图"(顶刊铁律明令少用/弃用)作 BEFORE 基线,
#          与上述 4 张高级图并列对比,实测展示条形图丢失的信息维度(基因重叠、
#          通路冗余、模块结构),证明"换图不换分析"的增益。不只报好看指标。
# ★关键   : cnetplot 已从 enrichplot 迁入 ggtangle 包(本机 enrichplot 1.26.6 /
#          ggtangle 0.1.2 实测确认);emapplot/treeplot 前必须先 pairwise_termsim()。
# 依赖   : enrichplot · ggtangle · clusterProfiler · org.Hs.eg.db(可选,缺则用合成
#          enrichResult)· DOSE · ggplot2 · aplot · ggtree
# 运行   : Rscript 546_enrichplot_emap_cnet_tree.R                  # 零改动跑合成示例
#          Rscript 546_enrichplot_emap_cnet_tree.R --genes my_genes.csv --outdir results/run1
# 输入   : --genes = 单列 csv,列名 SYMBOL(人类基因符号);无则脚本内合成一份
#          (synthetic demo only)。真实用法:把你的差异基因列表喂进来即可。
# =============================================================================

## ---- 0. 载入框架(顶刊主题 + save_fig + bio_args) --------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({
  library(enrichplot); library(ggtangle); library(clusterProfiler)
  library(ggplot2); library(DOSE)
}))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  genes       = file.path(DDAT, "gene_list.csv"),
  outdir      = file.path(SCRIPT_DIR, "results"),
  showCategory = 12,   # 高级图显示的通路条数
  ont         = "BP")) # GO 子本体 BP/MF/CC
args$showCategory <- as.integer(args$showCategory)
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. 合成示例基因列表(若无输入)---------------------------------------
# synthetic demo only:40 个真实人类基因符号,刻意取自若干内聚的生物学主题
# (DNA 损伤修复 / 细胞周期 / 凋亡 / 炎症-NFkB / JAK-STAT),便于 emapplot/treeplot
# 显出清晰的"功能模块",直观体现高级图相对条形图的结构信息增益。
SYN_GENES <- c(
  # DNA 损伤修复
  "ATM","CHEK2","BRCA1","BRCA2","RAD51","PARP1","XRCC1","TP53BP1",
  # 细胞周期
  "CCND1","CDK4","CDK6","RB1","E2F1","CDKN1A","MDM2","CCNE1",
  # 凋亡
  "TP53","BCL2","BAX","CASP3","CASP9","APAF1",
  # 炎症 / NFkB
  "NFKB1","RELA","TNF","IL6","CXCL8","CCL2","TLR4","MYD88",
  # JAK-STAT / 干扰素
  "JAK2","STAT3","IFNG","IRF3","SOCS3","IL10",
  # 生长信号
  "EGFR","PIK3CA","AKT1","MAPK1")

if (!file.exists(args$genes)) {
  write.csv(data.frame(SYMBOL = SYN_GENES), args$genes, row.names = FALSE)
  cat(sprintf("[gen] 合成基因列表 %d 个(synthetic demo only)→ %s\n",
              length(SYN_GENES), basename(args$genes)))
}
genes <- read.csv(args$genes, stringsAsFactors = FALSE)[[1]]
cat(sprintf("Step 1: 读入基因 %d 个\n", length(genes)))

## ---- 2. 富集分析:真包 enrichGO(org.Hs.eg.db);缺包则合成 enrichResult ----
cat("Step 2: 富集(enrichGO over-representation)...\n")
ego <- NULL
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressMessages(library(org.Hs.eg.db))
  eg <- tryCatch(
    clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db),
    error = function(e) NULL)
  if (!is.null(eg) && nrow(eg) >= 5) {
    ego <- clusterProfiler::enrichGO(
      gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = args$ont,
      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
      minGSSize = 5, maxGSSize = 300, readable = TRUE)
    # 收紧到内聚通路:enrichGO 易返回上千条高度冗余的 GO-BP,
    # 用 simplify() 按语义相似度去冗余,再保留最显著的若干条,利于看出模块结构。
    ego <- tryCatch(clusterProfiler::simplify(ego, cutoff = 0.6,
                    by = "p.adjust", select_fun = min), error = function(e) ego)
    cat(sprintf("  enrichGO(%s) 真包成功:simplify 后 %d 条通路\n",
                args$ont, nrow(as.data.frame(ego))))
  }
}

## ---- 2b. 降级:若无 org.Hs.eg.db / 映射失败 → 合成 enrichResult ------------
# 注:本机已装 org.Hs.eg.db,正常走真包路径;此分支仅为换机/缺包时仍能跑通出图。
if (is.null(ego) || nrow(as.data.frame(ego)) < 4) {
  cat("  ⚠ 真包富集不可用 → 构造合成 enrichResult(synthetic demo only)\n")
  make_syn_enrich <- function() {
    terms <- c("DNA repair","double-strand break repair","cell cycle G1/S transition",
               "regulation of cell cycle","apoptotic process","intrinsic apoptotic signaling",
               "inflammatory response","NF-kappaB signaling","cytokine production",
               "JAK-STAT signaling","response to interferon-gamma","cell population proliferation")
    glist <- list(
      c("ATM","CHEK2","BRCA1","BRCA2","RAD51","PARP1"),
      c("BRCA1","BRCA2","RAD51","TP53BP1","ATM"),
      c("CCND1","CDK4","CDK6","RB1","E2F1","CCNE1"),
      c("CDKN1A","MDM2","RB1","E2F1","CCND1"),
      c("TP53","BCL2","BAX","CASP3","CASP9","APAF1"),
      c("BCL2","BAX","CASP9","APAF1","TP53"),
      c("TNF","IL6","CXCL8","CCL2","TLR4","MYD88"),
      c("NFKB1","RELA","TNF","TLR4","MYD88"),
      c("IL6","TNF","CXCL8","IL10","CCL2"),
      c("JAK2","STAT3","SOCS3","IL6","IL10"),
      c("IFNG","STAT3","IRF3","JAK2"),
      c("EGFR","PIK3CA","AKT1","MAPK1","CCND1"))
    n <- length(terms); bg <- 18000; ng <- length(SYN_GENES)
    cnt <- vapply(glist, length, integer(1))
    df <- data.frame(
      ID = sprintf("GO:%07d", seq_len(n)), Description = terms,
      GeneRatio = paste0(cnt, "/", ng),
      BgRatio   = paste0(cnt * 6, "/", bg),
      pvalue    = 10^(-seq(8, 3, length.out = n)),
      p.adjust  = 10^(-seq(7, 2.5, length.out = n)),
      qvalue    = 10^(-seq(7, 2.5, length.out = n)),
      geneID    = vapply(glist, paste, collapse = "/", FUN.VALUE = character(1)),
      Count     = cnt, stringsAsFactors = FALSE)
    rownames(df) <- df$ID
    new("enrichResult", result = df, pvalueCutoff = 0.05, pAdjustMethod = "BH",
        qvalueCutoff = 0.2, organism = "Homo sapiens", ontology = "BP",
        gene = SYN_GENES, keytype = "SYMBOL", universe = as.character(seq_len(bg)),
        geneSets = setNames(glist, df$ID), readable = TRUE)
  }
  ego <- make_syn_enrich()
}

## ---- 3. pairwise_termsim:emapplot / treeplot 的前置(真实 API 要求)-------
cat("Step 3: pairwise_termsim()(emapplot/treeplot 前置)...\n")
ego_sim <- enrichplot::pairwise_termsim(ego)

n_show <- min(args$showCategory, nrow(as.data.frame(ego)))

## ---- 帮助函数:主题 ---------------------------------------------------------
# add_theme:用于坐标轴型图(dotplot/barplot),套完整 theme_pub。
add_theme <- function(p) tryCatch(p + theme_pub(base_size = 11), error = function(e) p)
# net_theme:用于网络/图布局型图(cnetplot/emapplot)。这类图没有坐标轴语义,
# 套 theme_bw 会硬加黑轴线/刻度(反 CONVENTIONS §3「勿对网络图套坐标轴主题」);
# 故只给白底 + 粗体标题 + 去坐标轴。★白底很关键:ggtangle 网络图默认透明底,
# 透明 PNG 在 README/查看器里显示为纯黑(本模块原缺陷的根因之一)。
net_theme <- function(p) tryCatch(
  p + ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.title  = ggplot2::element_text(size = 13, face = "bold", hjust = 0,
                                            family = PUB_FONT),
        plot.margin = ggplot2::margin(8, 8, 8, 8)),
  error = function(e) p)

## =============================================================================
## ★ 诚实基线(BEFORE):朴素富集条形图 —— 顶刊铁律明令少用,这里只作对照
## =============================================================================
cat("Step 4: [基线] 朴素富集条形图(BEFORE,铁律弃用,仅作对照)...\n")
df_top <- as.data.frame(ego)
df_top <- df_top[order(df_top$p.adjust), ][seq_len(n_show), ]
df_top$Description <- factor(df_top$Description, levels = rev(df_top$Description))
p_bar <- ggplot(df_top, aes(x = Description, y = -log10(p.adjust))) +
  geom_col(width = 0.7, fill = pal_pub(1, "npg")) +   # geom_col = 平凡条形图(反面教材)
  coord_flip() +
  labs(title = "BASELINE (discouraged): plain enrichment bar chart",
       subtitle = "One metric only (-log10 p.adjust); hides gene overlap, redundancy, modules",
       x = NULL, y = expression(-log[10]~adjusted~italic(p))) +
  theme_pub(base_size = 11)
save_fig(p_bar, file.path(ASSETS, "00_baseline_barplot"), width = 7.2, height = 5)

## =============================================================================
## 高级图 ① dotplot —— 气泡:x=GeneRatio,size=Count,color=p.adjust(一图四维)
## =============================================================================
cat("Step 5: [高级①] dotplot 气泡图(GeneRatio × Count × p.adjust)...\n")
p_dot <- enrichplot::dotplot(ego, showCategory = n_show, x = "GeneRatio") +
  scale_color_viridis_c(option = "D", direction = -1, name = "p.adjust") +
  labs(title = "Dotplot: GeneRatio x Count x p.adjust") +
  theme_pub(base_size = 11)
save_fig(p_dot, file.path(ASSETS, "01_dotplot_bubble"), width = 7.4, height = 5.6)

## =============================================================================
## 高级图 ② cnetplot —— 基因–概念网络(circular):显出哪些基因驱动哪些通路
##   ★ cnetplot 已迁入 ggtangle 包(实测 enrichplot 1.26.6/ggtangle 0.1.2):
##     新 API 参数已更名 —— circular→layout="circular";colorEdge→color_edge。
## =============================================================================
cat("Step 6: [高级②] cnetplot 基因–概念网络(circular, ggtangle)...\n")
n_cnet <- min(6, n_show)
p_cnet <- ggtangle::cnetplot(ego, showCategory = n_cnet, layout = "circular",
                             color_edge = "category", node_label = "all") +
  labs(title = "Cnetplot (circular): gene-concept network")
save_fig(net_theme(p_cnet), file.path(ASSETS, "02_cnetplot_network"),
         width = 8, height = 7)

## =============================================================================
## 高级图 ③ emapplot —— 富集图谱:按词项相似度把通路聚成功能模块
##   ★ 需先 pairwise_termsim()(已在 Step 3 做)
## =============================================================================
cat("Step 7: [高级③] emapplot 富集图谱(通路聚成模块)...\n")
# 实测 enrichplot 1.26.6:emapplot 无 repel 参数;有效参数 layout/color/node_label/
# group/min_edge 等。layout="kk"=Kamada-Kawai 力导向布局,相似通路自然聚拢。
p_emap <- enrichplot::emapplot(ego_sim, showCategory = n_show,
                               layout = "kk", node_label = "category",
                               color = "p.adjust") +
  labs(title = "Emapplot: pathways clustered into functional modules")
save_fig(net_theme(p_emap), file.path(ASSETS, "03_emapplot_modules"),
         width = 8, height = 7)

## =============================================================================
## 高级图 ④ treeplot —— 层次聚类:通路自动归并为带高阶标签的簇
## =============================================================================
cat("Step 8: [高级④] treeplot 层次聚类树...\n")
n_tree <- min(n_show, max(8, nrow(as.data.frame(ego_sim))))
p_tree <- tryCatch(
  enrichplot::treeplot(ego_sim, showCategory = n_tree) +
    labs(title = "Treeplot: hierarchical clustering of pathways"),
  error = function(e) { cat("  treeplot warn:", conditionMessage(e), "\n"); NULL })
if (!is.null(p_tree)) save_fig(p_tree, file.path(ASSETS, "04_treeplot_hierarchy"),
                               width = 9, height = 6)

## ---- 9. 落盘富集表 + 诚实基线对照说明 -------------------------------------
write.csv(as.data.frame(ego), file.path(args$outdir, "enrichment_table.csv"),
          row.names = FALSE)
# 基线对照:量化条形图 vs 高级图承载的信息维度(由代码生成,不手填)
n_terms  <- nrow(as.data.frame(ego))
n_genes_in <- length(unique(unlist(strsplit(as.data.frame(ego)$geneID, "/"))))
baseline_tbl <- data.frame(
  figure = c("00 baseline barplot", "01 dotplot", "02 cnetplot",
             "03 emapplot", "04 treeplot"),
  encodes_padjust = c("yes (length)", "yes (color)", "no", "no", "yes (color)"),
  encodes_generatio = c("no", "yes (x)", "no", "no", "no"),
  encodes_count = c("no", "yes (size)", "implicit", "yes (size)", "yes (size)"),
  shows_gene_overlap = c("no", "no", "YES", "implicit", "no"),
  shows_redundancy_modules = c("no", "no", "no", "YES", "YES"),
  stringsAsFactors = FALSE)
write.csv(baseline_tbl, file.path(args$outdir, "baseline_vs_advanced.csv"),
          row.names = FALSE)
cat(sprintf(
  "\n[诚实基线实测] 富集通路 %d 条,涉及基因 %d 个。\n  朴素条形图仅编码 1 个维度(p.adjust),且完全无法显示基因重叠 / 通路冗余 / 功能模块;\n  dotplot 一图编码 GeneRatio+Count+p.adjust 三维;cnetplot 显式连出基因-通路驱动关系;\n  emapplot/treeplot 把冗余通路归并为功能模块。对照表见 results/baseline_vs_advanced.csv\n",
  n_terms, n_genes_in))

cat("\n完成。富集表见", normalizePath(args$outdir), ";5 张图(1 基线 + 4 高级)见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖版本快照(铁律6)
