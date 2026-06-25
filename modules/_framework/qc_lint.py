#!/usr/bin/env python3
# =============================================================================
# qc_lint.py · 分析脚本质量机检 (12 铁律里可机检的子集)
# -----------------------------------------------------------------------------
# 用法:  python qc_lint.py <文件或目录>
#        扫 .R/.r/.py,报告:① 硬编码绝对路径/ setwd  ② 用了随机却没固定种子
#                          ③ 出图疑似非矢量(只 png 无 pdf/save_fig)  ④ 缺环境快照
# 退出码:发现 [高危] 返回 1,否则 0(可接入 git hook / CI)。
# 这是辅助,不替代人工对照 QUALITY_CHECKLIST.md。
# =============================================================================
import re, sys
from pathlib import Path

# 强制 UTF-8 输出,避免 Windows 控制台(GBK)下中文/中文路径乱码
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ABS_PATH = re.compile(r'["\'][A-Za-z]:[\\/]|["\'](?:/Users/|/home/|/mnt/|/data1?/)')
SETWD    = re.compile(r'\bsetwd\s*\(')
RAND_R   = re.compile(r'\b(sample|rnorm|runif|rbinom|kmeans|Rtsne|RunUMAP|RunTSNE|RunPCA|FindClusters|FindNeighbors|FindAllMarkers)\b')
SEED_R   = re.compile(r'set\.seed|seed\.use|random\.seed')
RAND_PY  = re.compile(r'\b(train_test_split|KMeans|leiden|louvain)\b|np\.random|umap\.|sc\.tl\.(umap|leiden|louvain)|random\.(sample|shuffle|choice|random)')
SEED_PY  = re.compile(r'random_state|np\.random\.seed|random\.seed|seed\.use|seed=')
GGSAVE   = re.compile(r'\bggsave\s*\(|\b(png|jpeg|tiff)\s*\(')
SAVEFIG  = re.compile(r'\bplt\.savefig\s*\(|\bfig\.savefig\s*\(')
SAVE_FIG_FW = re.compile(r'\bsave_fig\s*\(')
PDF      = re.compile(r'\.pdf|\.eps|\.svg|cairo_pdf')
ENVSNAP  = re.compile(r'sessionInfo|save_session|pip freeze|session_info')

def lint(p: Path):
    out = []
    txt = p.read_text(encoding="utf-8", errors="ignore")
    lines = txt.splitlines()
    is_r = p.suffix.lower() in (".r",)
    name = p.name.lower()
    cfg = name.startswith("config") or "_framework" in str(p).replace("\\", "/").lower()

    # ① 绝对路径 / setwd (逐行;config 与 _framework 自身允许放默认库路径)
    if not cfg:
        for i, ln in enumerate(lines, 1):
            s = ln.split("#")[0]
            if SETWD.search(s):
                out.append(("高危", i, "setwd() 硬编码工作目录;改用相对路径/bio_script_dir"))
            elif ABS_PATH.search(s):
                out.append(("高危", i, "疑似硬编码绝对路径;改用 config 相对路径"))

    # ② 随机却无种子 (文件级)
    rand = RAND_R if is_r else RAND_PY
    seed = SEED_R if is_r else SEED_PY
    if rand.search(txt) and not seed.search(txt):
        out.append(("中危", 0, "用到随机过程但未见固定种子(set.seed/random_state/seed.use)"))

    # ③ 出图疑似非矢量 (用了保存图却全文无 pdf/eps/svg/save_fig)
    saved = GGSAVE.search(txt) or SAVEFIG.search(txt)
    if saved and not (PDF.search(txt) or SAVE_FIG_FW.search(txt)):
        out.append(("中危", 0, "出图疑似只存位图;请用 save_fig() 或导出 PDF/EPS 矢量"))

    # ④ 缺环境快照 (信息级,长脚本才提示)
    if len(lines) > 40 and not ENVSNAP.search(txt) and not cfg:
        out.append(("提示", 0, "未见 sessionInfo/save_session;收尾建议记录依赖版本(铁律6)"))
    return out

def main():
    if len(sys.argv) < 2:
        print("用法: python qc_lint.py <文件或目录>"); return 2
    root = Path(sys.argv[1])
    if root.is_file():
        files = [root]
    else:
        files = [f for f in root.rglob("*") if f.suffix.lower() in (".r", ".py")]
    n_hi = n_file = 0
    for f in sorted(files):
        iss = lint(f)
        if not iss:
            continue
        n_file += 1
        print(f"\n# {f}")
        for lvl, ln, msg in iss:
            loc = f"L{ln}" if ln else "  "
            print(f"  [{lvl}] {loc:>5}  {msg}")
            if lvl == "高危":
                n_hi += 1
    print(f"\n== 扫描完成:{n_file} 个文件有可改进项,其中 [高危] {n_hi} 处 ==")
    return 1 if n_hi else 0

if __name__ == "__main__":
    sys.exit(main())
