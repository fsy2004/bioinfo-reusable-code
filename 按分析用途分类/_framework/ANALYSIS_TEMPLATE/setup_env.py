# =============================================================================
# setup_env.py · 项目初始化 (Python; 对应 R 的 00_setup.R)
# 用法:分析脚本顶部  from setup_env import *
# 提供:PROJ_ROOT、统一种子、目录、pubstyle、cache_step()、log_stat()、save_session()
# =============================================================================
import os, sys, random, pickle, subprocess
from pathlib import Path

# 1. 项目根 (铁律5:相对路径,禁绝对硬编码)
PROJ_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJ_ROOT))
from config import SEED, DIR_DATA, DIR_RESULTS, DIR_FIGURES, DIR_LOGS, FRAMEWORK_DIR, PARAMS

# 2. 统一随机种子 (铁律1) —— 各随机函数调用处仍要显式传 random_state=SEED
random.seed(SEED)
try:
    import numpy as np; np.random.seed(SEED)
except Exception:
    pass

# 3. 建标准目录
for _d in (DIR_DATA, DIR_RESULTS, DIR_FIGURES, DIR_LOGS):
    (PROJ_ROOT / _d).mkdir(parents=True, exist_ok=True)

# 4. 载入顶刊绘图框架 (铁律4/5:复用 pubstyle,不另写) —— 向上自动搜 _framework
_fw = Path(FRAMEWORK_DIR)
if not (_fw / "pubstyle.py").exists():
    _p = PROJ_ROOT
    for _ in range(6):
        if (_p / "_framework" / "pubstyle.py").exists():
            _fw = _p / "_framework"; break
        _p = _p.parent
if (_fw / "pubstyle.py").exists():
    sys.path.insert(0, str(_fw))
    try:
        from pubstyle import set_pub_style, pal, save_fig, panel_labels  # noqa: F401
        set_pub_style()
        print(f"[setup] 已载入框架: {_fw}")
    except Exception as e:
        print(f"[setup][警告] pubstyle 载入失败: {e}")
else:
    print("[setup][警告] 未找到 _framework/pubstyle.py,出图请手动设矢量+期刊配色")

# 5. 断点续跑 (铁律5:幂等,产物在则跳过)
def cache_step(name, fn, force=False):
    f = PROJ_ROOT / DIR_RESULTS / f"{name}.pkl"
    if f.exists() and not force:
        print(f"[cache] 跳过 {name} (已存在)"); return pickle.loads(f.read_bytes())
    print(f"[run ] {name} ..."); res = fn()
    f.write_bytes(pickle.dumps(res)); print(f"[done] {name} -> {f}"); return res

# 6. 关键统计值落盘 (铁律6:数字由代码生成,不手填)
def log_stat(key, value):
    with open(PROJ_ROOT / DIR_LOGS / "key_stats.tsv", "a", encoding="utf-8") as h:
        h.write(f"{key}\t{value}\n")
    print(f"[stat] {key} = {value}")

# 7. 可复现快照 (铁律6:锁定依赖版本)
def save_session():
    out = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                         capture_output=True, text=True).stdout
    (PROJ_ROOT / DIR_LOGS / "requirements_freeze.txt").write_text(out, encoding="utf-8")
    print("[setup] pip freeze -> logs/requirements_freeze.txt")

print(f"[setup] PROJ_ROOT={PROJ_ROOT} | SEED={SEED} | 目录就绪")
