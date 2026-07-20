# synthetic, for demo only

`transcripts.csv` 由 `573_proseg_cell_segmentation.py` 的 `simulate_transcripts()` 生成
(种子 573),**不是真实实验数据**,仅用于冒烟测试与生成 README 展示图。

模拟内容:120 个细胞的核质心(抖动网格,含彼此靠近的边界难例)、每细胞按细胞类型特异
的表达谱抽 60–160 条转录本、转录本半径按「核内均匀 + 胞质半正态长尾」两段分布,
另加 12% 均匀分布的环境/扩散噪声转录本(`true_cell_id = -1`)。

`cell_id` 只对落在核内的转录本赋值(其余为 -1),对应 proseg 官方 README 的输入要求:
转录本表须自带「初步的核/细胞分配」。`true_cell_id` / `true_cell_type` 是 ground truth,
真实数据没有这两列,仅本模块用于给基线打分。

文件删除后重跑脚本会自动重新生成。
