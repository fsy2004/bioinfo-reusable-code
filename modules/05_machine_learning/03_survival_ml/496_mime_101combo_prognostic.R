# =============================================================================
# 496_mime_101combo_prognostic.R
# 用途    : 10 种算法 × 101 组合 LOOCV 构建最优"生存"预后签名(按平均 C-index 选最优)
# 来源    : Mime  https://github.com/l-magnificence/Mime   ->  library(Mime1)
# 补库依据 : cat04 的"15种ML×175组合"是分类模型；论文3 (JCMM 2026, 肝癌 genistein) 用
#           Mime 式 101 组合做"生存"签名 (C-index 0.763)。补 Cox 生存集成这一缺口。
# 依赖    : devtools::install_github("l-magnificence/Mime")     # 先确认再装
# 输入    : list_train_vali_Data —— named list, 每队列 data.frame(ID, OS.time, OS, gene1, gene2 ...)
#                                   (如 TCGA 训练 + ICGC/GEO 验证；表达 log2)
#           genes —— 候选基因向量
# 输出    : 各算法组合 C-index 表 + 最优模型 + 风险评分 + 101 组合可视化(论文 Fig1k)
# =============================================================================
library(Mime1)
# res <- ML.Dev.Prog.Sig(
#   train_data           = list_train_vali_Data[[1]],
#   list_train_vali_Data = list_train_vali_Data,
#   candidate_genes      = genes,
#   mode                 = "all",                 # 10 算法的全部 101 组合
#   unicox.filter.for.candi = TRUE, unicox_p_cutoff = 0.05,
#   nodesize = 5, seed = 5201314)
#
# cindex_dis_all(res)            # 101 组合 C-index 棒棒糖/热图 = 论文 Fig1k 高级图
# cindex_dis_select(res, model = "最优组合名")
# survplot <- rs_sur(res, model_name = "最优组合名", dataset = names(list_train_vali_Data))
#
# 算法: RSF / Enet / Lasso / Ridge / StepCox / CoxBoost / plsRcox / SuperPC / GBM / survival-SVM
# 同框架可做: 免疫治疗应答(ML.Dev.Pred.Category) / 8 法特征选择取交集
