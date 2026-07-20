synthetic, for demo only
========================
synthetic_counts.npz —— 由 588_sccausalvi_causal_perturbation.py 的 make_example() 生成
(seed=42),不是真实实验数据,仅用于冒烟测试与生成展示图。

数组:
  X              float32 (960, 300)  Poisson 原始计数
  gene_names     str     (300,)      G000..G299
  cell_type      str     (960,)      CT_A / CT_B / CT_C
  condition      str     (960,)      ctrl / stim
  batch          int8    (960,)      0 / 1(与 condition 正交,不构成混杂)
  true_responder bool    (960,)      真值:stim 且属于会响应的细胞类型(CT_A/CT_B)

设计意图:处理效应是**细胞类型特异**的 —— CT_A 强响应(基因 0-39)、CT_B 中等响应
(基因 25-59)、CT_C 完全不响应(阴性对照细胞类型)。这样"全局差异表达"必然失真,
正好用来检验因果解耦方法与线性基线。true_responder 只用于事后评估,不参与任何打分。
