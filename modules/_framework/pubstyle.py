# -*- coding: utf-8 -*-
# =============================================================================
# pubstyle.py  ·  顶刊级绘图共享样式库 (Top-journal figure toolkit, Python)
# -----------------------------------------------------------------------------
# 用途   : 全库 Python 模块统一的 matplotlib 顶刊风格、期刊配色、矢量导出与多
#          panel 合成工具。import 即可获得一致的 Nature/Cell 级图风格。
# 提供   : set_pub_style()   一键设置 rcParams(Arial/Helvetica、矢量字体、刊级字号)
#          PAL               期刊离散配色 (npg/aaas/lancet/nejm/jama/vivid)
#          pal(n, name)      取 n 个配色(超长自动插值)
#          save_fig(fig...)  一次导出 矢量 PDF + 300dpi PNG
#          panel_labels(...) 给子图加 A/B/C 角标
# 依赖   : matplotlib, numpy (必需); 其余按模块自备。
# 约定   : 图中文字一律英文(投稿规范)。
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
}


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


def set_pub_style(base_size: int = 11, grid: bool = False) -> None:
    """一键应用顶刊风 rcParams:无衬线字体、矢量可编辑文字、干净坐标轴。"""
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
        "axes.linewidth": 0.8,
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
        "legend.frameon": False,
        "legend.fontsize": base_size - 2,
        "legend.title_fontsize": base_size - 1,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.prop_cycle": matplotlib.cycler(color=PAL["npg"]),
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
