# =============================================================================
# config.R · 分析项目集中配置  (复用 _framework 规范)
# -----------------------------------------------------------------------------
# 把所有"会变的东西"集中在这里:随机种子、路径、关键参数与阈值。
# 铁律(见 _framework/QUALITY_CHECKLIST.md):
#   · 不在分析脚本里硬编码路径/阈值;
#   · 阈值按【本数据集分布】定,勿套通用值(如 min.cells=200/mt<20% 一把过);
#   · 每个关键参数注明依据(文献/诊断图)。
# =============================================================================

## ---- 随机种子 (铁律1:可复现) ------------------------------------------------
SEED <- 42L            # 所有随机过程统一用它;Python 端 config.py 用【同一值】

## ---- 路径 (铁律5:不硬编码;全部相对项目根,由 00_setup.R 解析 PROJ_ROOT) -----
DIR_DATA    <- "data"        # 原始/中间输入
DIR_RESULTS <- "results"     # 中间产物/表格 (.gitignore 忽略)
DIR_FIGURES <- "figures"     # 顶刊级图 (矢量 PDF + 300dpi PNG)
DIR_LOGS    <- "logs"        # sessionInfo / 运行日志 / 关键统计值

## ---- 复用代码库框架 (铁律5:复用不重写) -------------------------------------
# _framework 位置:优先环境变量 BIOFW_DIR,否则用下面默认;00_setup 还会向上自动搜。
FRAMEWORK_DIR <- Sys.getenv("BIOFW_DIR",
  unset = "C:/Users/fsy/Desktop/bioinfo-reusable-code/modules/_framework")

## ---- 关键参数 (铁律1:不盲用默认值,逐项写依据) ------------------------------
# 下为单细胞示例;换你的分析类型时改写,但"留空待按数据定"的纪律保持。
PARAMS <- list(
  ## QC:阈值先留 NA,在 run_pipeline 里看【本数据】小提琴/分位数分布后再定,勿套通用值
  min_genes    = NA,    # nFeature 下尾(按分布拐点),勿默认 200
  max_mito_pct = NA,    # mt% 上尾(按分布拐点),勿默认 20
  ## 降维聚类:默认值仅起点,须用诊断图核实后写回依据
  n_pcs        = 30L,   # 依据:ElbowPlot/JackStraw(运行后核实并改)
  cluster_res  = 0.5,   # 依据:clustree 聚类稳定性(勿默认 0.8 一把过)
  ## marker / 注释
  markers_top  = 5L,    # 每类展示的 top marker 数
  annot_method = "manual+marker"  # 铁律2:注释必须 marker 验证,不裸信自动注释
)

## ---- 依赖锁定提示 (铁律6) ---------------------------------------------------
# 00_setup.R 运行后会把 sessionInfo() 写入 logs/;严格复现建议 renv::init()+snapshot()。
