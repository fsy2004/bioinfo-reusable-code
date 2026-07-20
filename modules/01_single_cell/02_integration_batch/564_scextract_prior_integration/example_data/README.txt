synthetic, for demo only — NOT real biological data.

synthetic_3batch.h5ad
  由 564_scextract_prior_integration.py 的 make_synthetic() 以固定种子 (seed=0) 生成;
  删除本文件后重跑主脚本会逐位重建同一份数据。

  1890 cells x 300 genes, 稠密 float32, 数值形如 log-normalized 表达。
  obs['batch']      : batch1 / batch2 / batch3
  obs['cell_type']  : Tcell / Bcell / Myeloid / Fibroblast (三批共享)
                      + RareStromal (仅存在于 batch3)

  设计意图(不是随机凑数):
   1. 批次效应(乘性 scale + 加性 shift + 基因级漂移)幅度 **大于** 细胞类型信号,
      模拟真实 scRNA 的处境;若信号压过噪声,所有方法都满分,基线就失去鉴别力。
   2. RareStromal 与 Fibroblast 共享 2/3 标志基因,且只出现在一个批次 ——
      过度校正会把它并进 Fibroblast。这是专门用来检验
      "aligning without flattening real differences" 的探针。
