synthetic, for demo only
------------------------
本目录所有文件由 594_colocboost_colocalization.R 的 Step 0 用固定种子 (SEED=42) 合成,
不是真实人类遗传数据,仅用于冒烟测试与生成 README 展示图。
删除本目录后重跑脚本会按同一种子逐字节重建。

sumstat_<trait>.csv : variant,pos,maf,beta,se,z,n  (4 个性状)
region_ld.csv       : 200x200 Pearson 相关矩阵,首行/首列为 variant 名
true_causal.csv     : 合成时植入的真因果变异 (rs060, rs150) —— 真实分析中不存在此文件
