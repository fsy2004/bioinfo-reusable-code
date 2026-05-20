# 17_高级结果图与闭环可视化

本模块收集一区/顶刊风格结果图的可复用工具入口和模板思路，重点用于把多组学证据组织成闭环，而不是简单美化柱状图。

推荐闭环：

`phenotype -> cell state -> communication / spatial niche -> perturbable target -> validation / prediction`

## 优先图型

1. Patient-level risk landscape  
   将机器学习预测概率、真实标签、数据集来源、细胞状态比例和 SHAP 贡献放到同一张患者级热图里。

2. Response circuit circos  
   用 circular tracks 展示 cell-state shift、SHAP、Scissor enrichment，用 chord links 展示通讯重排。

3. Target evidence wheel  
   把 virtual perturbation target、cell-state specificity、drug support、external validation 串成治疗证据图。

4. Ligand-target-response circuit  
   把 sender cell、ligand、receiver cell、target gene program 和 response phenotype 连起来。

5. Spatial niche response map  
   将单细胞状态映射到空间切片，并叠加 niche、communication vector 或 pathway activity。

## 文件

- `advanced_figure_tools.csv`: 推荐工具、用途、论文和 GitHub 地址。
- `download_advanced_figure_tools.ps1`: 一键 clone/update 相关工具到本地 `external_tools`。
- `literature_download_links_for_fdm.txt`: FDM 批量下载链接。
- `templates/`: 可移植模板说明。

## 使用建议

- 不建议把第三方工具源码直接提交进本仓库；用 `download_advanced_figure_tools.ps1` 拉取即可。
- 新项目只复制需要的模板脚本，保留本模块作为工具索引。
- 高级图必须服务于证据闭环；如果只是把条形图换成环形图，通常不值得放主图。

