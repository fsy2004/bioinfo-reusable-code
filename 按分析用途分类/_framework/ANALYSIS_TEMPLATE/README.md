# ANALYSIS_TEMPLATE · 标准分析项目骨架

开一个**新的多步分析项目**（不是单模块出图）时，复制本文件夹做起点。它把 12 条质量铁律（见 `../QUALITY_CHECKLIST.md`）落到代码结构里：统一种子、相对路径、断点续跑、关键统计落盘、依赖快照，并复用 `_framework` 的顶刊绘图工具。

## 怎么用
```bash
# 1) 复制骨架到你的新项目目录(例:SSc 某个子分析)
cp -r _framework/ANALYSIS_TEMPLATE  ~/Desktop/我的新分析

# 2) 改 config(R 用 config.R / Python 用 config.py):写 SEED、路径、关键参数依据
#    阈值类(min_genes/mt% 等)先留空,在 run_pipeline 里按本数据分布定后写回

# 3) 把数据放进 data/,逐步实现 run_pipeline 的 Step,跑:
Rscript run_pipeline.R          # 或在 VSCode 工作区里逐块运行
#   Python 脚本里:  from setup_env import *

# 4) 收尾:对照 ../QUALITY_CHECKLIST.md 逐条自查 + 机检
python ../qc_lint.py .
```

## 文件
| 文件 | 作用 |
|---|---|
| `config.R` / `config.py` | 集中配置:随机种子、路径、关键参数(逐项写依据) |
| `00_setup.R` / `setup_env.py` | 初始化:解析项目根、设种子、建目录、载 `_framework`、提供 `cache_step()`/`log_stat()`/`save_session()` |
| `run_pipeline.R` | 主流程骨架:分步 + 断点续跑 + 每步内嵌【铁律自查点】+ 矢量出图 |

## 标准目录(运行后自动建)
```
data/      输入        results/   中间产物(.gitignore)
figures/   矢量图      logs/      sessionInfo / key_stats.tsv
```

## 关键工具(setup 提供)
- `cache_step("名", { 计算 })` — 产物存在则跳过，断点续跑（铁律5）
- `log_stat("键", 值)` — 关键统计值写 `logs/key_stats.tsv`，文稿数字从这里取（铁律6）
- `save_session()` — 写 `sessionInfo` / `pip freeze`，锁定依赖（铁律6）
- 绘图直接用 `_framework` 的 `theme_pub()`/`pal_pub()`/`save_fig()`（R）、`set_pub_style()`/`pal()`/`save_fig()`（Py）

> 框架定位：setup 优先用环境变量 `BIOFW_DIR`，否则按 config 默认路径，再否则从项目向上自动搜 `_framework`。换机器时设 `BIOFW_DIR` 即可。
