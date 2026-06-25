# -*- coding: utf-8 -*-
# =============================================================================
# pubstyle.py  ·  顶刊级绘图共享样式库 (Top-journal figure toolkit, Python)
# -----------------------------------------------------------------------------
# 用途   : 全库 Python 模块统一的 matplotlib 顶刊风格、期刊配色、矢量导出与多
#          panel 合成工具。import 即可获得一致的 Nature/Cell 级图风格。
# 提供   : set_pub_style()   一键设置 rcParams(Arial/Helvetica、矢量字体、刊级字号)
#          PAL               期刊离散配色 (npg/aaas/lancet/nejm/jama/vivid + 色盲安全 okabe_ito)
#          pal(n, name)      取 n 个配色(超长自动插值)
#          CMAP_CONT / CMAP_DIVERGE  连续(viridis)/ 发散(RdBu)色图名
#          save_fig(fig...)  一次导出 矢量 PDF + 300dpi PNG
#          panel_labels(...) 给子图加 A/B/C 角标
#          NATURE_W1 / NATURE_W2  单栏(3.5in≈90mm)/双栏(7.0in≈180mm)宽度常量
# 依赖   : matplotlib, numpy (必需); 其余按模块自备。
# 约定   : 图中文字一律英文(投稿规范)。
# -----------------------------------------------------------------------------
# Nature/Cell 出图规范(默认对齐 nature.com 官方 figure guidelines):
#   尺寸 单栏 90mm(NATURE_W1) / 双栏 180mm(NATURE_W2) / 最大高 170mm;
#   字号 最终图 5–7pt、panel 角标 8pt 粗体;字体 Helvetica/Arial 全文一致;
#   线宽 最细线 ≥1pt(本库 axes.linewidth 已对齐);颜色 RGB,避免色盲红绿对
#   → 离散优先 pal("okabe_ito")、连续 viridis、发散 RdBu(CMAP_DIVERGE);
#   矢量 pdf.fonttype=42 文字可编辑,save_fig() 出 PDF+PNG。
# =============================================================================
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 期刊离散配色板(与 R 端 theme_pub.R 完全一致,跨语言统一观感)
PAL = {
    "npg":    ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
               "#91D1C2", "#DC0000", "#7E6148", "#B09C85"],
    "aaas":   ["#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021",
               "#5F559B", "#A20056", "#808180", "#1B1919"],
    "lancet": ["#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91",
               "#AD002A", "#ADB6B6", "#1B1919"],
    "nejm":   ["#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD",
               "#FFDC91", "#EE4C97"],
    "jama":   ["#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"],
    "vivid":  ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
               "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"],
    # 色盲安全板 Okabe-Ito (Wong 2011 Nature Methods);离散变量首选,满足 Nature
    # 「避免红绿混淆」要求。npg/aaas 等含红绿对,严格色盲场景请改用本板。
    "okabe_ito": ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                  "#CC79A7", "#000000", "#999999"],
}
PAL["cb"] = PAL["okabe_ito"]   # 别名:cb = colorblind-safe

# 发散配色板 RdBu(低=蓝 中=白 高=红),与 R 端 PUB_DIVERGE 一致(ColorBrewer RdBu 反向)
DIVERGE = ["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
           "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]

# Nature 出图尺寸常量(英寸;1in=25.4mm)
NATURE_W1 = 3.5    # 单栏 ≈ 90mm
NATURE_W2 = 7.0    # 双栏 ≈ 180mm
NATURE_HMAX = 6.7  # 最大高 ≈ 170mm


def pal(n: int | None = None, name: str = "npg") -> list[str]:
    """取 n 个期刊配色;n 超过板长时在板内平滑插值扩展。"""
    base = PAL.get(name, PAL["npg"])
    if n is None or n <= len(base):
        return base[: (n if n else len(base))]
    cmap = LinearSegmentedColormap.from_list("_p", base, N=n)
    return [matplotlib.colors.to_hex(cmap(i / (n - 1))) for i in range(n)]


def _pick_font() -> str:
    """优先 Arial/Helvetica(期刊标准无衬线),不可用则降级。"""
    import matplotlib.font_manager as fm
    avail = {f.name for f in fm.fontManager.ttflist}
    for f in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if f in avail:
            return f
    return "DejaVu Sans"


def set_pub_style(base_size: int = 11, grid: bool = False, palette: str = "npg") -> None:
    """一键应用顶刊风 rcParams:无衬线字体、矢量可编辑文字、干净坐标轴。

    base_size : 基础字号(屏幕预览用 11;Nature 单栏最终图建议 7)。
    grid      : 是否保留极淡网格(默认 False,更干净)。
    palette   : 默认配色循环板名(默认 npg;色盲安全场景传 "okabe_ito")。
    """
    font = _pick_font()
    matplotlib.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,          # 矢量字体可编辑(投稿要求)
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": [font, "DejaVu Sans"],
        "font.size": base_size,
        "axes.titlesize": base_size + 1,
        "axes.titleweight": "bold",
        "axes.labelsize": base_size,
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "axes.linewidth": 1.0,        # Nature:最细线 ≥1pt
        "axes.spines.top": False,     # 去顶/右边框 → 期刊常见干净风
        "axes.spines.right": False,
        "axes.grid": grid,
        "grid.color": "#E8E8E8",
        "grid.linewidth": 0.5,
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.labelsize": base_size - 1,
        "ytick.labelsize": base_size - 1,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "legend.frameon": False,
        "legend.fontsize": base_size - 2,
        "legend.title_fontsize": base_size - 1,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.prop_cycle": matplotlib.cycler(color=PAL.get(palette, PAL["npg"])),
    })


def panel_labels(axes: Sequence[plt.Axes], labels: str = "ABCDEFGH",
                 size: int = 16, dx: float = -0.08, dy: float = 1.04) -> None:
    """给一组子图左上角加粗体 A/B/C 角标(多panel合成图规范)。"""
    for ax, lab in zip(np.ravel(axes), labels):
        ax.text(dx, dy, lab, transform=ax.transAxes, fontsize=size,
                fontweight="bold", va="top", ha="right")


def save_fig(fig: plt.Figure, file: str | Path, dpi: int = 300) -> None:
    """一次导出 矢量 PDF + 300dpi PNG(README 预览)。file 可不含扩展名。"""
    stem = str(file)
    for ext in (".pdf", ".png"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    Path(stem).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", dpi=dpi, bbox_inches="tight")


# 连续量统一首选 viridis(色盲友好、灰度印刷稳健)
CMAP_CONT = "viridis"
CMAP_DIVERGE = "RdBu_r"
