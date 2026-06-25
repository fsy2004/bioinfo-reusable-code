# =============================================================================
# 编号       : R019
# 脚本名     : 来源.R
# 分类       : 06_immune_infiltration
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 免疫分析所需来源/辅助函数或数据说明脚本。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : 
# 内容重复   : 与 R017 字节级相同（任选其一即可）
# 整理时间   : 2026-05-10
# =============================================================================
##这段代码来自于 CIBERSORT R script v1.03，由 Aaron M. Newman（斯坦福大学）开发。详细的使用说明、依赖要求以及许可条款可以在 CIBERSORT 官方网站上找到，例如在 CIBERSORT 官网 和 许可说明 中。
#原始来源是 CIBERSORT R script v1.03，其详细信息如下：

#作者与机构：
#该脚本由 Aaron M. Newman 编写，隶属于斯坦福大学（联系方式：amnewman@stanford.edu）。

#功能简介：
#CIBERSORT 是一种用于从混合基因表达数据中估计各个细胞类型比例的计算方法。脚本主要利用支持向量回归（SVR）对输入数据进行分解。其核心算法会尝试多个 nu 参数（例如 0.25、0.5、0.75），以寻找最优模型，并进一步计算相关系数和均方根误差（RMSE），从而评估模型的表现。
#将该 R 脚本保存到本地目录中。

#perm 参数控制置换检验的次数（用于计算 p 值，建议设置为 100 或以上），

#QN 参数控制是否对混合数据进行量化归一化（默认 TRUE）。

#其他说明：

#签名矩阵构建功能当前未包含在该 R 脚本中，若需要完整功能（包括签名矩阵构建），建议使用 CIBERSORT 的 Java 版本。

#脚本会将分析结果输出文件，并附带一个矩阵对象，其中包含各混合样本的细胞比例估计、p 值、相关系数和 RMSE。

#许可证信息可参考 CIBERSORT 许可证。

#这些细节可以在 CIBERSORT 的官方网站（http://cibersort.stanford.edu）以及相关教程中找到更详细的描述。
#' CIBERSORT R script v1.03
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export