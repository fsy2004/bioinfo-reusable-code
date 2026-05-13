# =============================================================================
# 编号       : R008
# 脚本名     : 01GEO数据处理.有gene symbol情况.R
# 分类       : 03_GEO转录组整理与差异分析
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 整理 GEO 转录组表达矩阵与基因注释，生成后续分析所需矩阵。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : 
# 整理时间   : 2026-05-10
# =============================================================================
# 1. 设置工作目录
setwd("C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/12.GEO数据整理/常见情况..转录组数据整理-")  # 设置工作目录

# 2. 设置参数：用户指定基因信息所在行数（例如第11列）并转换为0起始索引
inputRow <- 11                              # 用户指定使用第?行数据（1起始）
targetColIdx <- inputRow - 1                  # 转换为0起始索引，用于平台文件中目标列定位

# 3. 自动识别GSE和GPL开头的文件
cat("正在搜索GSE和GPL文件...\n")
allFiles <- list.files(pattern = "\\.(txt|TXT)$")  # 获取所有txt文件

# 查找GSE开头的文件
gseFiles <- allFiles[grepl("^GSE", allFiles, ignore.case = TRUE)]
if (length(gseFiles) == 0) {
  stop("错误：未找到GSE开头的文件！")
} else if (length(gseFiles) > 1) {
  cat("找到多个GSE文件：\n")
  print(gseFiles)
  exprFilePath <- gseFiles[1]  # 默认使用第一个
  cat("使用第一个文件：", exprFilePath, "\n")
} else {
  exprFilePath <- gseFiles[1]
  cat("找到GSE文件：", exprFilePath, "\n")
}

# 查找GPL开头的文件
gplFiles <- allFiles[grepl("^GPL", allFiles, ignore.case = TRUE)]
if (length(gplFiles) == 0) {
  stop("错误：未找到GPL开头的文件！")
} else if (length(gplFiles) > 1) {
  cat("找到多个GPL文件：\n")
  print(gplFiles)
  platformFilePath <- gplFiles[1]  # 默认使用第一个
  cat("使用第一个文件：", platformFilePath, "\n")
} else {
  platformFilePath <- gplFiles[1]
  cat("找到GPL文件：", platformFilePath, "\n")
}

# 定义输出文件的路径
outputFilePath <- "geneMatrix.csv"            # 最终输出文件路径（改为CSV格式）

# 4. 读取表达数据文件（自动从ID_REF行开始）
cat("正在加载表达数据文件：", exprFilePath, "...\n")  # 输出提示信息

# 先读取所有行，找到ID_REF所在的行号
allLines <- readLines(exprFilePath)
# 查找包含ID_REF的行（可能在第一列的任意位置，前面可能有引号或其他字符）
idRefLine <- which(grepl("ID_REF", allLines, ignore.case = FALSE))[1]  # 找到第一个包含ID_REF的行

if (is.na(idRefLine)) {
  stop("错误：未在文件中找到包含ID_REF的行！")
}

cat("找到ID_REF所在行：第", idRefLine, "行\n")
cat("将跳过前", idRefLine - 1, "行样本信息\n")

# 从ID_REF行开始读取数据
exprData <- read.delim(exprFilePath,          # 读取文件
                       header = FALSE,         # 暂不使用自动表头
                       sep = "\t",             # 使用制表符分隔
                       quote = "\"",           # 使用双引号包裹字符型数据
                       skip = idRefLine - 1,   # 跳过ID_REF之前的所有行
                       comment.char = "")      # 不忽略任何字符

# 将第一行（ID_REF行）设置为列名
colnames(exprData) <- exprData[1, ]
exprData <- exprData[-1, ]                    # 删除第一行（因为已经作为列名）
colnames(exprData)[1] <- "ProbeID"            # 将第一列重命名为 "ProbeID"

# 重置行名
rownames(exprData) <- NULL

cat("表达数据加载完成，共", ncol(exprData)-1, "个样本。\n")  # 输出样本数（减去探针ID列）

# 5. 读取平台文件（GPL.txt）
cat("正在加载平台文件：", platformFilePath, "...\n")  # 输出提示信息
platformData <- read.delim(platformFilePath,    # 读取平台文件
                           header = FALSE,        # 无表头
                           sep = "\t",            # 以制表符分隔
                           quote = "\"",          # 使用双引号包裹字符型数据
                           comment.char = "#",    # 忽略以"#"开头的注释行
                           stringsAsFactors = FALSE)  # 不将字符转换为因子

# 6. 构建探针与基因符号映射（逐行处理并显示进度）
cat("开始处理平台文件数据...\n")         # 输出提示信息
totalRows <- nrow(platformData)               # 获取平台文件总行数
geneMapping <- list()                         # 初始化空列表用于存储映射：探针ID -> 基因符号
pbPlat <- txtProgressBar(min = 0, max = totalRows, style = 3)  # 创建进度条，范围从0到总行数

for (row in 1:totalRows) {                     # 遍历平台文件的每一行
  # 若整行全为空，或第一列以 "ID" 或 "!" 开头，则跳过该行
  if (all(platformData[row, ] == "") || grepl("^(ID|\\!)", platformData[row, 1])) {
    setTxtProgressBar(pbPlat, row)             # 更新进度条
    next
  }
  
  # 检查当前行是否包含目标列（targetColIdx + 1）
  if (ncol(platformData) >= (targetColIdx + 1)) {
    rawGene <- platformData[row, targetColIdx + 1]  # 提取目标列的原始基因信息
    # 如果原始信息非空且不包含空格（只包含一个单词）
    if (rawGene != "" && !grepl(".+\\s+.+", rawGene)) {
      # 若包含"///"，只取左侧部分；同时去除引号
      cleanGene <- sub("(.+?)///(.+)", "\\1", rawGene)
      cleanGene <- gsub('"', '', cleanGene)
      # 建立映射：将当前行第一列（探针ID）映射到处理后的基因符号
      geneMapping[[ as.character(platformData[row, 1]) ]] <- cleanGene
    }
  }
  setTxtProgressBar(pbPlat, row)             # 更新进度条
}
close(pbPlat)                                 # 关闭进度条
cat("\n平台文件处理完成，成功建立映射数：", length(geneMapping), "\n")  # 输出映射建立完成信息

# 7. 合并表达数据与平台映射数据
cat("正在合并表达数据与平台映射数据...\n")  # 输出提示信息
# 从表达数据中提取探针ID列并保存
probeIDs <- exprData[["ProbeID"]]
# 取出除探针ID以外的表达数据部分
exprDataCore <- exprData[, -1, drop = FALSE]

# 将表达数据列转换为数值型（非常重要！）
cat("正在将表达数据转换为数值型...\n")
for (col in colnames(exprDataCore)) {
  exprDataCore[[col]] <- as.numeric(as.character(exprDataCore[[col]]))
}

# 重新组合ProbeID和数值型表达数据
exprDataWithID <- data.frame(ProbeID = probeIDs, exprDataCore, stringsAsFactors = FALSE)

# 将平台映射转换为数据框，包含探针ID和对应的基因符号
mappingDF <- data.frame(ProbeID = names(geneMapping),
                        geneSymbol = unlist(geneMapping),
                        stringsAsFactors = FALSE)
# 利用 merge 函数按"ProbeID"匹配，合并表达数据与基因映射信息
mergedData <- merge(exprDataWithID, mappingDF, by = "ProbeID")
cat("合并完成，共获得", nrow(mergedData), "行有效数据。\n")  # 输出合并后数据行数

# 8. 按基因分组，计算各样本的平均表达值（逐基因处理显示进度）
cat("正在按基因分组并计算各样本平均表达值...\n")  # 输出提示信息
# 获取所有样本列名称（排除ProbeID和新合并的geneSymbol列）
sampleCols <- setdiff(colnames(mergedData), c("ProbeID", "geneSymbol"))
# 根据基因符号对数据进行分组
geneGroups <- split(mergedData, mergedData$geneSymbol)
# 获取所有不同基因名称
geneNamesUnique <- names(geneGroups)
totalGenes <- length(geneNamesUnique)         # 总基因数
resultList <- vector("list", totalGenes)        # 初始化列表保存计算结果
pbAgg <- txtProgressBar(min = 0, max = totalGenes, style = 3)  # 创建聚合进度条

for (i in seq_along(geneNamesUnique)) {         # 遍历每个基因组
  geneName <- geneNamesUnique[i]                # 当前基因名称
  groupData <- geneGroups[[i]]                  # 取出该基因对应的所有探针数据
  # 对样本列（除ProbeID和geneSymbol）计算均值，忽略缺失值
  means <- colMeans(groupData[, sampleCols, drop = FALSE], na.rm = TRUE)
  # 将当前基因名称与计算的均值合并为一行
  resultList[[i]] <- c(geneSymbol = geneName, means)
  setTxtProgressBar(pbAgg, i)                   # 更新进度条
}
close(pbAgg)                                   # 关闭聚合进度条
# 将结果列表合并为数据框
aggData <- do.call(rbind, resultList)
# 转换为数据框，并保持字符型基因名称
aggData <- as.data.frame(aggData, stringsAsFactors = FALSE)
# 对除geneSymbol外的数值型列转换为数值型（因do.call返回字符型）
for (col in colnames(aggData)[-1]) {
  aggData[[col]] <- as.numeric(aggData[[col]])
}
# 对结果按照基因名称排序
aggData <- aggData[order(aggData$geneSymbol), ]
cat("数据聚合完成，共计算出", nrow(aggData), "个基因的平均表达值。\n")  # 输出聚合结果

# 9. 写入输出文件（geneMatrix.csv）
cat("正在写入输出文件：", outputFilePath, "...\n")  # 输出写入提示信息
# 利用 write.csv 写出数据，使用逗号分隔、不写行号
write.csv(aggData, file = outputFilePath, row.names = FALSE)
cat("输出文件生成成功：", outputFilePath, "\n")  # 输出完成提示
