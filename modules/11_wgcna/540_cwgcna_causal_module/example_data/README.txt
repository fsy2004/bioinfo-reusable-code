synthetic, for demo only — 由主脚本 540_cwgcna_causal_module.R 自动生成(set.seed(42))。

expr.csv   : 150 样本 × 160 基因表达;4 个对照模块,基因前缀标注真身:
             CAU_* = 真上游因模块(driver → 模块 → trait)
             EFF_* = 真下游果模块(trait → 模块)
             CON_* = 仅混杂相关模块(confounder 同时影响模块与 trait)
             NUL_* = 阴性背景模块(与 trait 无关)
traits.csv : sample / trait(目标性状) / driver(外生工具,锚定因果方向) / confounder

真实数据请按 README ① 的列规格准备;删除这两个 CSV 后重跑脚本会重新生成。
合成数据仅用于冒烟测试与生成展示图,勿用于任何真实结论。
