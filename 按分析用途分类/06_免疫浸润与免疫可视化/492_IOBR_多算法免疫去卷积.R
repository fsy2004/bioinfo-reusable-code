# =============================================================================
# 492_IOBR_多算法免疫去卷积.R
# 用途    : 一站式多算法免疫细胞去卷积 + 跨方法交叉验证（巨噬 M1/M2 等一致性论证）
# 来源    : IOBR  https://github.com/IOBR/IOBR  (Zeng et al., Front Immunol 2021)
# 补库依据 : 覆盖矩阵 cat06 原仅 CIBERSORT；论文4 (IJMS 2026, Alzheimer M2 髓系) 用 IOBR
#           七法 (xCell/CIBERSORT/EPIC/MCPcounter/quanTIseq/TIMER/ESTIMATE) 交叉验证
#           M2 巨噬随 Braak 分期富集 —— 多法一致 = 可信度论证。一个包补齐 cat06 缺口。
# 依赖    : remotes::install_github("IOBR/IOBR")   # 首次安装；按惯例先确认再装
# 输入    : eset  —— 基因(行, HGNC symbol) × 样本(列) 表达矩阵, 建议 TPM 或 log2(TPM+1)
#           pdata —— data.frame(ID, group)   group = 分组/分期
# 输出    : 合并去卷积矩阵 + 堆叠柱状图 + 组间比较 + 多法一致性
# =============================================================================
suppressMessages({library(IOBR); library(tidyverse)})

# eset  <- readRDS("expr_tpm.rds")          # rownames = symbol
# pdata <- read.csv("pdata.csv")            # 列: ID, group

## 1) 多算法去卷积 ----------------------------------------------------------
methods <- c("cibersort","epic","mcpcounter","xcell","quantiseq","timer","estimate")
decon <- lapply(methods, function(m)
  tryCatch(deconvo_tme(eset = eset, method = m, arrays = FALSE, perm = 200),
           error = function(e){message(m, " 失败: ", conditionMessage(e)); NULL}))
names(decon) <- methods
tme <- purrr::reduce(Filter(Negate(is.null), decon),
                     ~ dplyr::inner_join(.x, .y, by = "ID"))   # 按样本(ID)合并

## 2) 堆叠柱（以 CIBERSORT LM22 为例）---------------------------------------
p_bar <- cell_bar_plot(input = decon$cibersort, title = "CIBERSORT (LM22)")

## 3) 组间比较（以巨噬 M2 为例，跨方法 facet）-------------------------------
m2_cols <- grep("M2", colnames(tme), value = TRUE, ignore.case = TRUE)
df <- tme %>% inner_join(pdata, by = "ID") %>%
  pivot_longer(any_of(m2_cols), names_to = "method_cell", values_to = "score")
p_box <- ggplot(df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = .4) +
  facet_wrap(~ method_cell, scales = "free_y") +
  ggpubr::stat_compare_means(label = "p.signif") + theme_bw()

## 4) 跨方法一致性（同一细胞类型不同算法的 Spearman 相关）------------------
# cors <- cor(tme[, m2_cols], method = "spearman", use = "pairwise")
# pheatmap::pheatmap(cors, display_numbers = TRUE, main = "M2 跨算法一致性")

# ggsave("IOBR_M2_box.pdf", p_box, width = 9, height = 6)
# 备注: ESTIMATE 额外给 immune/stromal/purity 评分; deconvo_tme 亦支持 method="ips"
