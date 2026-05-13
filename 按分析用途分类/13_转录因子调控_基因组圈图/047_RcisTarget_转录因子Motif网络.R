# ==========================================================================
# 脚本名     : 转录因子基序调控网络.R
# 分类       : 13_转录因子调控_基因组圈图
# 项目来源   : 从压缩包 107.转录因子（基序转录调控）网络.rar 整理
# 原始文件   : 107.转录因子（基序转录调控）网络\TF.R
# 用途       : 基于 RcisTarget 对候选基因集进行 motif/转录因子富集分析，并构建 motif-TF-gene 调控网络。
# 结果图     : AUC直方图；显著基因富集曲线；incidence matrix曲线；igraph多布局网络图；visNetwork交互网络HTML；Sankey网络HTML
# 非肿瘤消化适配: 适合。可用于非肿瘤消化系统差异基因、WGCNA模块基因、单细胞marker的转录调控解释。
# 主要 R 包  : RcisTarget; igraph; visNetwork; DT; reshape2; htmlwidgets; data.table; networkD3
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# 设置工作目录
setwd("H:\\常用分析生信\\107.转录因子（基序转录调控）网络")

# 加载必要的库
library(RcisTarget)      # 用于 Motif 富集分析
library(igraph)          # 用于绘制网络图
library(visNetwork)      # 用于创建交互式网络可视化
library(DT)              # 用于创建动态数据表
library(reshape2)        # 用于数据重塑
library(htmlwidgets)     # 用于保存 HTML 小部件
library(data.table)      # 提供高效的数据处理功能

# 设置文件路径
featherFile <- "hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"  # Motif 排名数据库文件
geneFile <- "gene.txt"                                                                    # 基因列表文件
motifAnnotationFile <- "motifAnnotations_hgnc.RData"                                      # Motif 注释文件

# 导入 Motif 排名数据库
motifRankings <- importRankings(featherFile)  # 导入排名数据库

# 读取基因列表
geneList1 <- read.table(geneFile, stringsAsFactors = FALSE)[, 1]  # 从文件中读取基因列表的第一列
geneLists <- list(geneSet = geneList1)                            # 将基因列表存储在列表中，名称为 geneSet

# 加载 Motif 注释文件
load(motifAnnotationFile)  # 加载 Motif 注释对象

# 验证 Motif 注释对象已加载
if (!exists("motifAnnotations")) {
  stop("motifAnnotations 对象未找到，请检查加载的 RData 文件。")  # 如果对象不存在，停止运行并提示错误
}

# 进行 Motif 富集分析
motifEnrichmentTable_wGenes <- cisTarget(
  geneLists,                  # 基因列表
  motifRankings,              # Motif 排名数据库
  motifAnnot = motifAnnotations  # Motif 注释对象
)

# 计算 AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)  # 计算 Motif 的 AUC 值
auc <- getAUC(motifs_AUC)[1, ]                   # 获取第一个基因集的 AUC 值

# 绘制并保存 AUC 直方图
pdf("AUC_Histogram.pdf")  # 打开 PDF 设备，开始保存图形
hist(auc, main = "AUC Histogram", xlab = "AUC", breaks = 100,    # 绘制 AUC 直方图
     col = "#ff000050", border = "darkred")
nesThreshold <- mean(auc) + 3 * sd(auc)     # 计算 NES 阈值
abline(v = nesThreshold, col = "blue")      # 添加垂直线表示阈值
dev.off()                                   # 关闭 PDF 设备

# 添加 Motif 注释
motifEnrichmentTable <- addMotifAnnotation(
  motifs_AUC,             # AUC 结果
  nesThreshold = 3.5,     # NES 阈值
  motifAnnot = motifAnnotations,    # Motif 注释对象
  highlightTFs = list(gene = "ACBD3")  # 高亮显示的转录因子
)

# 查看结果
print(class(motifEnrichmentTable))                   # 打印结果表的类
print(dim(motifEnrichmentTable))                     # 打印结果表的维度
print(head(motifEnrichmentTable[, -"TF_lowConf", with = FALSE]))  # 查看结果表的前几行

# 添加显著基因信息
motifEnrichmentTable_wGenes <- addSignificantGenes(
  motifEnrichmentTable,    # 富集结果表
  rankings = motifRankings,  # Motif 排名数据库
  geneSets = geneLists       # 基因列表
)
print(dim(motifEnrichmentTable_wGenes))  # 打印包含显著基因的结果表维度

# 设置基因集名称
geneSetName <- "geneSet"  # 指定要使用的基因集名称

# 绘制显著基因的曲线并保存为 PDF
selectedMotifs <- c("cisbp__M6275", sample(motifEnrichmentTable$motif, 2))  # 选择要分析的 Motif
pdf("SignificantGenes_Curves.pdf")   # 打开 PDF 设备
par(mfrow = c(2, 2))                 # 设置绘图区域为 2 行 2 列
getSignificantGenes(
  geneLists[[geneSetName]],   # 基因集
  motifRankings,              # Motif 排名数据库
  signifRankingNames = selectedMotifs,  # 选择的 Motif 名称
  plotCurve = TRUE,           # 绘制累积分布曲线
  maxRank = 5000,             # 最大排名
  genesFormat = "none",       # 不返回基因列表
  method = "aprox"            # 使用近似方法
)
dev.off()  # 关闭 PDF 设备

# 添加 Motif 标识的 Logo（需要相关数据库）
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable)  # 为结果添加 Motif Logo

# 提取前 12 个结果用于展示
resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:12, ]  # 获取前 12 条结果

# 创建并保存包含 Logo 的数据表
datatable_file <- "MotifEnrichmentTable.html"  # 指定保存文件的名称
datatable(
  resultsSubset[, -c("enrichedGenes", "TF_lowConf"), with = FALSE],  # 删除不需要的列
  escape = FALSE,    # 允许 HTML 标签（用于显示 Logo）
  filter = "top",    # 在表格顶部添加过滤器
  options = list(pageLength = 5)  # 设置每页显示 5 行
) %>%
  saveWidget(file = datatable_file)  # 保存为 HTML 文件

# 提取注释的转录因子
annotatedTfs <- lapply(
  split(motifEnrichmentTable_wGenes$TF_highConf, motifEnrichmentTable$geneSet),  # 按基因集分组
  function(x) {
    genes <- gsub(" \\(.*\\). ", "; ", x, fixed = FALSE)  # 清理转录因子名称字符串
    genesSplit <- unique(unlist(strsplit(genes, "; ")))   # 分割并去除重复的转录因子
    return(genesSplit)  # 返回转录因子列表
  }
)
print(annotatedTfs$geneSet)  # 打印注释的转录因子列表

# 获取显著的 Motif 名称
signifMotifNames <- motifEnrichmentTable$motif[1:3]  # 选择前 3 个显著的 Motif

# 获取关联矩阵并绘制曲线，保存为 PDF
pdf("IncidenceMatrix_Curves.pdf")  # 打开 PDF 设备
incidenceMatrixResult <- getSignificantGenes(
  geneLists[[geneSetName]],       # 基因集
  motifRankings,                  # Motif 排名数据库
  signifRankingNames = signifMotifNames,  # 显著的 Motif 名称
  plotCurve = TRUE,               # 绘制累积分布曲线
  maxRank = 5000,                 # 最大排名
  genesFormat = "incidMatrix",    # 返回关联矩阵
  method = "aprox"                # 使用近似方法
)
dev.off()  # 关闭 PDF 设备
incidenceMatrix <- incidenceMatrixResult$incidMatrix  # 提取关联矩阵

# 创建边的数据框
edges <- melt(incidenceMatrix)  # 将矩阵转换为长格式数据
edges <- edges[which(edges[, 3] == 1), 1:2]  # 过滤出关联值为 1 的边
colnames(edges) <- c("from", "to")  # 重命名列名为 'from' 和 'to'

# 创建 igraph 对象
g <- graph_from_incidence_matrix(incidenceMatrix, directed = FALSE)  # 从关联矩阵创建无向图

# 设置节点颜色
V(g)$color <- ifelse(V(g)$type, "red", "green")  # Motif 节点为红色，基因节点为绿色

# 根据 Motif 名称设置边的颜色
edge_colors <- sapply(E(g), function(e) {
  motif_name <- V(g)$name[get.edges(g, e)[1]]  # 获取边关联的 Motif 名称
  if (grepl("cisbp__M6275", motif_name)) {
    return("red")     # 如果是指定的 Motif，设置相应的颜色
  } else if (grepl("cisbp__M6279", motif_name)) {
    return("blue")
  } else if (grepl("cisbp__M0062", motif_name)) {
    return("green")
  } else if (grepl("cisbp__M4575", motif_name)) {
    return("yellow")
  } else if (grepl("cisbp__M4476", motif_name)) {
    return("pink")
  } else {
    return("grey")    # 其他 Motif 的边为灰色
  }
})
E(g)$color <- edge_colors  # 应用边颜色设置

# 绘制不同布局的网络图并保存为 PDF
layouts <- grep("^layout_", ls("package:igraph"), value = TRUE)[-1]  # 获取所有布局函数名称
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]  # 排除不适用的布局

pdf("Graph_Layouts.pdf")  # 打开 PDF 设备
par(mfrow = c(3, 5), mar = c(1, 1, 1, 1))  # 设置绘图区域和边距
for (layout in layouts) {
  l <- do.call(layout, list(g))  # 计算布局
  plot(g, edge.arrow.mode = 0, layout = l, main = layout, vertex.label.cex = 0.7)  # 绘制网络图
}
dev.off()  # 关闭 PDF 设备

# 准备 visNetwork 的节点和边数据
motifs <- unique(as.character(edges$from))  # 获取唯一的 Motif 节点
genes <- unique(as.character(edges$to))     # 获取唯一的基因节点
nodes <- data.frame(
  id = c(motifs, genes),                               # 节点 ID
  label = c(motifs, genes),                            # 节点标签
  title = c(motifs, genes),                            # 鼠标悬停时显示的标题
  shape = c(rep("diamond", length(motifs)), rep("ellipse", length(genes))),  # 节点形状
  color = c(rep("purple", length(motifs)), rep("skyblue", length(genes)))    # 节点颜色
)

# 创建 visNetwork 可视化并保存为 HTML
visNet <- visNetwork(nodes, edges) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)  # 设置交互选项
saveWidget(visNet, "Network.html")  # 保存网络图为 HTML 文件
library(networkD3)
library(htmlwidgets)
# 你的代码
links <- data.frame(source=edges$from, target=edges$to, value=1)
nodes <- data.frame(name=unique(c(links$source, links$target)))
links$source <- match(links$source, nodes$name)-1
links$target <- match(links$target, nodes$name)-1
sank <- sankeyNetwork(Links=links, Nodes=nodes, Source="source", Target="target", Value="value", NodeID="name")
saveWidget(sank, "sankey.html")
# 获取布局名
layouts <- grep("^layout_", ls("package:igraph"), value = TRUE)[-1] 
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)] 

# 对每个布局，单独画图保存PDF
for (layout in layouts) {
  l <- do.call(layout, list(g))
  pdf(paste0(layout, ".pdf"), width=8, height=6)  # 每个布局生成单独PDF
  plot(g, edge.arrow.mode = 0, layout = l, main = layout, vertex.label.cex = 0.7)
  dev.off()
}
#