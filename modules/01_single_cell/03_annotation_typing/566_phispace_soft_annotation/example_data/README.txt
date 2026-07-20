synthetic, for demo only —— 本目录 4 个 CSV 全部由 set.seed(566) 的泊松模型生成,
不是真实测序数据,仅用于冒烟测试与生成 assets/ 展示图。

reference_counts.csv    200 genes x 240 cells,第一列 gene,其余列为参考集细胞计数
reference_metadata.csv  cell, cell_type (TypeA/TypeB/TypeC,每类 80 个纯细胞)
query_counts.csv        200 genes x 300 cells,第一列 gene
query_metadata.csv      cell, true_group, true_TypeA, true_TypeB, true_TypeC
                        真实组成比例已知(纯 A/B/C 各 60 个 + A-B 连续过渡态 120 个),
                        用于评估连续得分能否还原真实混合比例。真实数据无此列时脚本自动跳过评估。
