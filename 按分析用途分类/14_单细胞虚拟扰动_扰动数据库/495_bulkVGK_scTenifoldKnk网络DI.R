# =============================================================================
# 495_bulkVGK_scTenifoldKnk网络DI.R
# 用途    : 在 bulk 共表达网络上做"虚拟基因敲除"——无需 scRNA 也能算 in-silico KO 影响
# 来源    : scTenifoldKnk https://github.com/cailab-tamu/scTenifoldKnk + igraph 网络内聚
# 补库依据 : cat14 的虚拟扰动均面向 scRNA；论文2 (Biomedicines 2026, 糖尿病足 KIF13A)
#           把虚拟敲除套到 WGCNA/igraph bulk 网络上，用"差异影响 DI"在两组(愈合/不愈合)
#           各自敲除每个基因、比较网络内聚下降 → 找拓扑驱动基因。对共病 bulk 数据极有用。
# 依赖    : install.packages(c("scTenifoldKnk","igraph","WGCNA"))   # 先确认再装
# 输入    : expr_A, expr_B —— 两组(如 nonhealed / healed) 基因×样本表达矩阵
# 输出    : 每基因 Differential Impact (DI) 排名 + 驱动发现图(常规DE显著性 vs DI)
# =============================================================================
suppressMessages({library(igraph)})

cohesion <- function(adj){                         # 网络内聚 = 边密度 与 全局聚类系数 的均值
  g <- graph_from_adjacency_matrix(abs(adj) > 0.6, mode = "undirected", diag = FALSE)
  mean(c(edge_density(g), transitivity(g, type = "global")), na.rm = TRUE)
}
vgk_impact <- function(expr, top = 2000){          # 逐基因敲除 → 内聚下降量
  v   <- head(order(apply(expr, 1, var), decreasing = TRUE), top)
  m   <- expr[v, ]; base <- cohesion(cor(t(m)))
  sapply(rownames(m), function(g) base - cohesion(cor(t(m[setdiff(rownames(m), g), ]))))
}
# impA <- vgk_impact(expr_A); impB <- vgk_impact(expr_B)
# DI   <- impA[names(impB)] - impB                  # 差异影响：组特异拓扑重要性 (DI 高 = 驱动)
# 驱动发现图(论文 Fig1C): plot(-log10(deg_pval), DI)  左上角 = 高影响但常规DE漏掉的 hub
#
# 更严格的张量版(纯单组 KO): scTenifoldKnk::scTenifoldKnk(countMatrix, gKO = "KIF13A")
