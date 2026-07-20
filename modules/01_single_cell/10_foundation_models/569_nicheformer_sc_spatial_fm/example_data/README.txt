synthetic, for demo only —— 本目录 4 个 CSV 全部为脚本生成的合成数据,仅用于冒烟测试与展示图,
不来自任何真实组织/公开数据集,不得用于任何生物学结论。

spatial_query_expression.csv        700 细胞 × 60 基因,原始计数(行=cell_id,列=gene)
spatial_query_meta.csv              cell_id, x, y, niche, cell_type
dissociated_reference_expression.csv 500 细胞 × 60 基因,原始计数(含模拟平台偏移)
dissociated_reference_meta.csv      cell_id, cell_type

生成逻辑(见 README.md「输入数据」段):4 个细胞类型程序 + 3 个空间生态位;
生态位对表达的直接影响很弱(+0.35 log),但生态位之间的细胞类型组成不同 ——
因此 niche 标签主要要靠**邻域上下文**而非细胞自身表达才能恢复。
