# 统一框架规范 (Conventions)

本库所有模块遵循同一套 **turnkey(开箱即跑)** 与 **顶刊图** 规范。本文件是改造与新增模块的唯一标准。

## 1. 模块目录结构

每个模块是一个自带文件夹的"即跑单元":

```
<类别>/<NNN_模块名>/
├─ <NNN_模块名>.R 或 .py     # 主脚本(turnkey)
├─ README.md                 # 5 段式文档(见模板)
├─ example_data/             # 小型合成示例输入(可直接跑通)
│  └─ <input>.csv
├─ assets/                   # 提交进库的顶刊级展示图 PNG(README 引用)
│  └─ <fig>.png
└─ results/                  # 运行时生成,.gitignore 忽略
```

`assets/` **必须提交**(README 配图来源);`results/` 永不提交。

## 2. Turnkey 运行约定

- **零改动即跑**:`Rscript <模块>.R` 或 `python <模块>.py` —— 默认读 `example_data/`,写 `results/`。
- **换数据即跑**:`Rscript <模块>.R --input data/你的.csv --outdir results/run1`。
- **禁止** `setwd("绝对路径")`;一律用脚本相对路径定位(R 见 `bio_script_dir()`,Py 见 `ROOT = Path(__file__).parent`)。
- 顶部保留 **参数区**(原有风格),但默认值必须指向 `example_data/`,且关键参数支持 `--key value` 覆盖。
- 保留原有 **中文头部元信息块** 与 **`cat("Step X")` 进度**;**分析逻辑一字不改**,只标准化 I/O + 升级出图。

## 3. 顶刊图风格规范

统一 `source(".../_framework/theme_pub.R")`(R)或 `from pubstyle import *`(Py)。

- **主题**:`theme_pub()` / `set_pub_style()` —— Arial/Helvetica、去网格、黑轴线、刊级字号。
- **配色**:离散用期刊板 `pal_pub("npg"/"lancet"/"aaas"/...)`;连续量用 **viridis**;发散用 RdBu。禁用 Set3/Pastel 等"软"配色。
- **导出**:`save_fig()` 一次出 **矢量 PDF + 300dpi PNG**;README 用 PNG。
- **单图为主**:每张图**独立成文件**(便于投稿时自行在 AI/Illustrator 拼版);`compose_panels()`(R)/`panel_labels()`(Py)仅作**可选**工具,默认**不自动输出合成图**。
- **图中文字英文**(投稿规范),代码注释中文。
- **不要简陋默认图**:包默认图(enrichplot barplot、TwoSampleMR mr_*_plot、base plot.roc、FeaturePlot 灰红等)优先重绘为精修版,底层分析结果不变。

## 4. 文档规范

每模块 `README.md` 必含 5 段(见 `TEMPLATE_module_README.md`):
①输入数据(规格表+样例) ②方法/原理 ③用途 ④特点/亮点 ⑤输出结果图(清单+内嵌预览图);
并附 **运行命令** 与 **依赖安装**。主页 README 汇总所有模块并提供"数据→图"索引。

## 5. 示例数据规范

`example_data/` 内为 **小型合成数据**,仅用于冒烟测试与生成展示图,文件头/README 注明 `synthetic, for demo only`。严格匹配 README 的输入规格(列名、类型、命名约定)。
