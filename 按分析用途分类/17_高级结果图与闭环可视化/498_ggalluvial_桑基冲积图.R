# =============================================================================
# 498_ggalluvial_桑基冲积图.R
# 用途    : 多列流向图(药物→靶点→通道 / 配体→受体→细胞)——补 cat17 缺的桑基/冲积图
# 来源    : ggalluvial  https://github.com/corybrunson/ggalluvial   (CRAN)
# 补库依据 : advanced_figure_tools.csv 无桑基/alluvial；论文1 Fig4G(药物→hub→离子通道)、
#           论文4 Fig7A/B(信号 river) 都需要。
# 依赖    : install.packages("ggalluvial")
# 输入    : df —— 长表, 每行一条流, 列 = 各层(axis1/axis2/axis3) + freq
# 输出    : 三列冲积图(PDF)
# =============================================================================
library(ggplot2); library(ggalluvial)

# df <- data.frame(drug = ..., hub = ..., channel = ..., freq = 1)
# ggplot(df, aes(axis1 = drug, axis2 = hub, axis3 = channel, y = freq)) +
#   geom_alluvium(aes(fill = drug), width = 1/12, alpha = .8) +
#   geom_stratum(width = 1/12, fill = "grey90", color = "grey40") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
#   scale_x_discrete(limits = c("Drug","Hub","Ion channel"), expand = c(.05,.05)) +
#   theme_minimal() + theme(legend.position = "none")
#
# 替代方案:
#   networkD3::sankeyNetwork(...)         # 交互式 HTML 桑基图
#   CellChat::netAnalysis_river(object)   # 单细胞信号"河流图"(论文4 Fig7A/B 即此)
