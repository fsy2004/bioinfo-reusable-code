# =============================================================================
# config.py · 分析项目集中配置 (Python; 与 config.R 对齐,SEED 同值)
# 铁律见 _framework/QUALITY_CHECKLIST.md:不硬编码路径/阈值;阈值按本数据分布定。
# =============================================================================
import os

# ---- 随机种子 (铁律1) ----  与 config.R 用【同一值】
SEED = 42

# ---- 路径 (铁律5:相对项目根,由 setup_env.py 解析 PROJ_ROOT) ----
DIR_DATA, DIR_RESULTS, DIR_FIGURES, DIR_LOGS = "data", "results", "figures", "logs"

# ---- 复用代码库框架 (铁律5) ----  优先环境变量 BIOFW_DIR
FRAMEWORK_DIR = os.environ.get(
    "BIOFW_DIR", r"C:/Users/fsy/Desktop/bioinfo-reusable-code/modules/_framework")

# ---- 关键参数 (铁律1:不盲用默认,逐项写依据) ----  单细胞示例
PARAMS = dict(
    min_genes=None,      # 按本数据 nFeature 分布定,勿默认 200
    max_mito_pct=None,   # 按本数据 mt% 分布定,勿默认 20
    n_pcs=30,            # 依据:方差累计/elbow(运行后核实改)
    cluster_res=0.5,     # 依据:聚类稳定性(勿默认一把过)
    markers_top=5,
)
