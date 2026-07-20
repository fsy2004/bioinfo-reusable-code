# =============================================================================
# 494_genki_vgae_ko.py
# 用途    : 基于图变分自编码器(VGAE)的单细胞虚拟基因敲除——给现有 scTenifoldKnk/
#           CellOracle 体系再加一种"图自编码器"范式 (用 KL 散度量化扰动强度)
# 来源    : GenKI  https://github.com/yjgeno/GenKI  (Yang et al., NAR 2023)
# 补库依据 : 覆盖矩阵 cat14 有 scTenifoldKnk/CellOracle/GEARS/Squidiff，无 GenKI；
#           论文1 (IJMS 2026, 离子通道CRC) 用 GenKI 做 hub→离子通道 in-silico KO。
# 依赖    : pip install GenKI    # 自带 CPU 版 torch/torch-geometric；GPU 需先装匹配 CUDA 版
#           Python>=3.10；输入需 .h5ad (AnnData)
# 输入    : adata.h5ad —— 单细胞表达(野生型/对照)
# 输出    : 每基因 KL 散度 + 百分位排名 + bootstrap 显著性
# =============================================================================
from GenKI import GenKI

# ---- 高阶一步式 ----
ranked = GenKI.from_h5ad(
    "adata.h5ad",
    target_gene=["KCNQ2"],            # TODO: 要敲除的基因
).run(epochs=100, seed=8096, n_permutations=100)
print(ranked.head(20))               # 受该基因敲除影响最大的基因(按 KL/距离排名)

# ---- 低阶可控流程(自定义 GRN/阈值时) ----
# from GenKI.dataLoader import DataLoader
# from GenKI.train import VGAE_trainer
# from GenKI import utils
# dl   = DataLoader(adata, target_gene=["KCNQ2"], rebuild_GRN=True, GRN_file_dir="GRN", n_cpus=4)
# d_wt = dl.load_data();        d_ko = dl.load_kodata()
# tr   = VGAE_trainer(d_wt, epochs=100, seed=8096); tr.train()
# z_wt = tr.get_latent_vars(d_wt);  z_ko = tr.get_latent_vars(d_ko)
# dis  = utils.get_distance(z_ko, z_wt, by="KL")     # 每基因 KL 散度
# res  = utils.get_generank(adata, dis, rank=True)   # 排名表
#
# 提示: KL 分布重尾——论文用"百分位排名 + 100 次置换 bagging(>=95%一致)"而非参数检验；
#       负对照取 50 个随机非 hub 基因过同一冻结模型，Mann-Whitney U 做组间比较。
