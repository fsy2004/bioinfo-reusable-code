# ==========================================================================
# 脚本名     : 基因染色体位置圈图.R
# 分类       : 13_转录因子调控_基因组圈图
# 项目来源   : 从压缩包 330.基因在染色体位置分布圈图.rar 整理
# 原始文件   : 330.基因在染色体位置分布圈图\基因在染色体位置分布圈图.R
# 用途       : 读取目标基因及基因组坐标，使用 circlize 绘制基因在染色体上的圈图分布。
# 结果图     : 染色体ideogram圈图；基因标签圈图；基因位置信息表
# 非肿瘤消化适配: 适合但偏展示。非肿瘤消化系统候选基因集可用作补充可视化。
# 主要 R 包  : circlize; RColorBrewer
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================

# 设置工作目录
setwd("H:\\常用分析生信\\330.基因在染色体位置分布圈图")
#https://www.ncbi.nlm.nih.gov/gene/3
# 加载 circlize 包，用于绘制圆形可视化图表
library(circlize)

# 加载 RColorBrewer 包，用于获取专业的配色方案
library(RColorBrewer)

# 定义包含基因表达信息的输入文件路径
geneFile="gene.csv"

# 定义包含基因位置信息的输入文件路径
posFile="geneREF.csv"


# ============ 数据读取和预处理部分 ============

# 读取基因位置信息文件
# header=T: 第一行是列名
# sep=",": 使用逗号作为分隔符
# check.names=F: 保持原始的列名不进行修改
genepos=read.table(posFile, header=T, sep=",", check.names=F)

# 重命名列，使用更简洁的名称
# 原始列名: Gene, Chr, Start, End -> 新列名: genename, chr, start, end
colnames(genepos)=c('genename','chr','start','end')

# 重新排列列的顺序，将基因名称放在最后
genepos=genepos[,c('chr','start','end','genename')]

# 设置行名称为基因名称，便于后续按基因名称进行索引
row.names(genepos)=genepos[,'genename']

# 读取要绘制的基因列表文件
# header=F: 没有列名
# sep=",": 使用逗号作为分隔符
geneRT=read.table(geneFile, header=F, sep=",", check.names=F)

# 按照 geneFile 中列出的基因顺序，筛选 genepos 中对应的数据
# 这样可以只绘制指定的基因在圆形图上的位置
genepos=genepos[as.vector(geneRT[,1]),]

# 将处理后的基因位置数据保存为 bed0 格式，用于后续绘制
bed0=genepos


# ============ 颜色配置部分 ============

# 利用 colorRampPalette 从 RColorBrewer 的 Dark2 调色板中扩展出足够数量的颜色
# Dark2 提供深色且易区分的颜色，更适合生物学数据可视化，共8种基础颜色
# 这里扩展为24种颜色，用于显示24条染色体
nColors <- 24
myColors <- colorRampPalette(brewer.pal(8, "Dark2"))(nColors)


# ============ 圆形图绘制部分 ============

# 创建并保存为 PDF 文件，设置图表大小为 6x6 英寸
pdf(file="circlize.pdf", width=6, height=6)

# 清空之前的圆形图设置
circos.clear()

# 使用 hg38 人类基因组初始化圆形图
# plotType=NULL: 不绘制轨道，仅显示染色体框架
circos.initializeWithIdeogram(species="hg38", plotType=NULL)

# 绘制内层轨道：显示每条染色体的颜色和名称标签
# ylim = c(0, 1): Y轴范围
# panel.fun: 定义轨道内的绘制函数
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  # 获取当前扇形区域对应的染色体名称
  chr = CELL_META$sector.index

  # 获取当前扇形区域的X轴范围
  xlim = CELL_META$xlim

  # 获取当前扇形区域的Y轴范围
  ylim = CELL_META$ylim

  # 根据染色体编号获取对应的颜色（从自定义的调色板 myColors 中选择）
  # 处理标准染色体编号（chr1-chr22）以及性染色体（chrX->23, chrY->24）
  if (chr == "chrX") {
    color_idx <- 23
  } else if (chr == "chrY") {
    color_idx <- 24
  } else {
    # 提取数字部分进行转换
    chr_index <- as.numeric(gsub("chr", "", chr))
    if (is.na(chr_index)) {
      # 对于其他非标准染色体名称，使用循环方式分配颜色
      chr_index <- (match(chr, unique(CELL_META$all.sector.index)) - 1) %% nColors + 1
    }
    color_idx <- min(chr_index, nColors)
  }

  # 在轨道上绘制矩形，使用 Dark2 调色板中的颜色填充（深色系，显示效果更佳）
  circos.rect(xlim[1], 0, xlim[2], 1, col=myColors[color_idx])

  # 在矩形中心添加染色体名称标签（白色，居内侧）
  circos.text(mean(xlim), mean(ylim), chr, cex=0.6, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height=0.15, bg.border = NA)

# 绘制基因组图谱（中层轨道）：显示标准的人类基因组理想图
# species="hg38": 使用hg38人类基因组版本
# track.height=mm_h(6): 设置轨道高度为6毫米
circos.genomicIdeogram(species = "hg38", track.height=mm_h(6))

# 绘制基因标签（外层轨道）：将选定的基因名称显示在圆形图外侧
# bed0: 包含基因位置和名称的数据
# labels.column=4: 使用第4列（基因名称）作为标签
# side="inside": 标签显示在圆形内侧
# cex=0.8: 设置文字大小为正常大小的80%
circos.genomicLabels(bed0, labels.column=4, side = "inside", cex=0.8)

# 清空圆形图的所有设置和布局，准备下一个绘制任务
circos.clear()

# 保存PDF文件
dev.off()


# ============ 导出数据部分 ============

# 将基因位置信息重新整理为便于阅读的格式
# bed0 中的列顺序为: chr, start, end, genename
# 将其调整为: genename, chr, start, end
output_data <- bed0[, c('genename', 'chr', 'start', 'end')]

# 重置行名称，使其与基因名称对应
rownames(output_data) <- NULL

# 输出为 CSV 文件
# row.names=FALSE: 不输出行号
# quote=FALSE: 不添加引号
# 输出文件名: gene_position_info.csv
write.csv(output_data, file="gene_position_info.csv", row.names=FALSE, quote=FALSE)

# 打印确认信息
cat("基因位置信息已导出到: gene_position_info.csv\n")
cat("共包含", nrow(output_data), "个基因\n")

# 完成