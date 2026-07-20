"""572 · CellVQ — 离散「细胞词表」(vector-quantised codebook)单细胞基础模型。

本模块做两件事:

  1) **可跑基线(始终执行, CPU, 只用本机已有依赖)**:用 PCA + k-means 码本把每个细胞
     压成一个整数 "cell code",量化 VQ 的三件事——码本利用率(perplexity)、量化重构损失、
     以及**离散码相对连续嵌入损失了多少标签信息**(kNN 交叉验证 head-to-head)。
     k-means 码本不是随便选的替代:VQ-VAE 的码本更新本身就是对 encoder 输出做
     k-means / EMA 聚类,所以「PCA 嵌入 + k-means 码本」正是去掉预训练 encoder 之后
     的 VQ 朴素下界。任何"离散词表更好"的主张都必须先赢过它。

  2) **CellVQ 守卫式封装(--cellvq-repo)**:官方 CellVQ 需要 clone 仓库 + 从 ModelScope
     下载预训练权重 + GPU,不可能在本机零依赖跑通。因此这里**不复刻模型**,只做
     环境体检并打印**已核实的真实调用命令**(逐字来自官方 inference.py / run_example.sh)。
     ⚠️ 模型规模上下游说法不一:**论文摘要写 "model parameters totaling 500 million"**,
     **仓库 README 写 "511 million parameters"**;两者都写 "68 million cells" 预训练。
     本模块不在两者间取舍,只如实标注出处。

上游(已核实):
  论文 Wang J, Tan C, Gao Z, Shao S, Liu S, Li SZ. Illuminating cell states by a
    comprehensive and interpretable single cell foundation model. Nat Commun 2026;17:4037.
    doi:10.1038/s41467-026-70071-5 · PMID 41839876 · PMC13139411
  仓库 https://github.com/A4Bio/CellVQ (默认分支 master)
  API 读取来源(本地克隆 + 逐符号定位到源码行,非凭 README 推测):
    inference.py                        CLI 全部 8 个参数与默认值(第 15-22 行)
    model/load.py:124                   def load_model_frommmf(best_ckpt_path, mode='m1',
                                                               params=None, device='cuda')
    model/pretrainmodels/model.py:156   def get_cellcode(self, x, padding_label,
                                            encoder_position_gene_ids, output_attentions=False,
                                            **kwargs) -> (geneemb, cell_code, ...)
    inference.py:74                     torch.save([cellembs, cell_codes], args.save_path)
    run_example.sh / install.sh / README.md
  仓库文件存在性经 GitHub Contents API 核对(examples/cluster_19264.h5ad、
  model/OS_scRNA_gene_index.19264.tsv 均确实在仓库内)。
"""
from __future__ import annotations

import argparse, json, os, sys, warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data", "synthetic_counts.csv")

# 框架统一出图风格
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import set_pub_style, pal, save_fig, CMAP_CONT, NATURE_W1, NATURE_W2  # noqa: E402

SEED = 572


def save_session() -> dict:
    """依赖版本快照(铁律6:可复现)。落进 572_summary.json。"""
    import platform
    snap = {"python": platform.python_version(), "platform": platform.platform()}
    for m in ("numpy", "pandas", "sklearn", "matplotlib", "torch"):
        try:
            snap[m] = __import__(m).__version__
        except Exception:
            snap[m] = "not installed"
    return snap


# ------------------------------------------------------------------ 数据 / 预处理
def load_counts(path: str):
    """读 cells x genes 计数 CSV(首列 cell_id, 第二列 cell_type, 其余为基因)。"""
    import pandas as pd
    df = pd.read_csv(path, index_col=0, comment="#")
    if "cell_type" not in df.columns:
        raise SystemExit(f"{path} 缺少 cell_type 列(基线要用它做标签信息评估)")
    labels = df["cell_type"].astype(str).values
    X = df.drop(columns=["cell_type"]).values.astype("float64")
    return X, labels, list(df.columns[1:])


def embed(X, n_pcs: int = 30):
    """标准 scanpy 式预处理 → PCA 连续嵌入(这是被离散化的对象)。"""
    import numpy as np
    from sklearn.decomposition import PCA
    # CP10K + log1p:与 CellVQ inference.py 对 singlecell 的 pre_normalized=False 分支一致
    depth = X.sum(axis=1, keepdims=True)
    Xn = np.log1p(X / np.clip(depth, 1e-9, None) * 1e4)
    Xn = (Xn - Xn.mean(0)) / np.clip(Xn.std(0), 1e-9, None)
    n_pcs = min(n_pcs, min(Xn.shape) - 1)
    return PCA(n_components=n_pcs, random_state=SEED).fit_transform(Xn)


# ------------------------------------------------------------------ VQ 基线
def quantize(Z, K: int):
    """k-means 码本 = VQ 码本的朴素下界。返回 (码 index, 码本, 量化后嵌入)。"""
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=K, n_init=10, random_state=SEED).fit(Z)
    codes = km.labels_
    return codes, km.cluster_centers_, km.cluster_centers_[codes]


def codebook_perplexity(codes, K: int) -> float:
    """码本利用率的标准度量:exp(H(p))。等于 K 表示码全用满,远小于 K 表示码本坍缩。"""
    import numpy as np
    p = np.bincount(codes, minlength=K) / len(codes)
    p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


def knn_label_accuracy(F, labels, n_neighbors: int = 15) -> float:
    """5 折交叉验证 kNN 分类准确率——衡量该表示里还剩多少细胞类型信息。

    分折在 fit 之外做,避免用同一批细胞既建表示又评估造成的乐观偏倚在折内进一步放大。
    """
    import numpy as np
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    from sklearn.neighbors import KNeighborsClassifier
    k = min(n_neighbors, np.bincount(np.unique(labels, return_inverse=True)[1]).min() - 1)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    sc = cross_val_score(KNeighborsClassifier(n_neighbors=max(k, 1)), F, labels, cv=cv)
    return float(sc.mean())


def run_baseline(Z, labels, k_grid, outdir: str):
    """扫码本大小 K:标签信息保留 / 码本 perplexity / 量化重构 R²,并与连续嵌入对比。"""
    import numpy as np, pandas as pd
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    cont_acc = knn_label_accuracy(Z, labels)          # 连续嵌入的天花板(未离散化)
    tss = float(((Z - Z.mean(0)) ** 2).sum())

    rows, per_k = [], {}
    for K in k_grid:
        codes, book, Zq = quantize(Z, K)
        # 离散码本身作为表示:用 one-hot 码送进 kNN,信息量完全来自"哪个码"
        onehot = np.eye(K)[codes]
        rows.append(dict(
            K=K,
            code_knn_acc=knn_label_accuracy(onehot, labels),
            quant_knn_acc=knn_label_accuracy(Zq, labels),      # 用码本向量代替原嵌入
            perplexity=codebook_perplexity(codes, K),
            recon_r2=1.0 - float(((Z - Zq) ** 2).sum()) / tss,
            ari=adjusted_rand_score(labels, codes),
            nmi=normalized_mutual_info_score(labels, codes),
            n_used=int(len(np.unique(codes))),
        ))
        per_k[K] = codes
    sweep = pd.DataFrame(rows)
    sweep.to_csv(os.path.join(outdir, "572_codebook_sweep.csv"), index=False)
    return sweep, cont_acc, per_k


# ------------------------------------------------------------------ 出图(无条形图)
def fig_sweep(sweep, cont_acc, outdir_assets):
    """码本大小扫描:折线+点(三联),含连续嵌入基线水平虚线。"""
    import matplotlib.pyplot as plt
    c = pal(4)
    fig, axes = plt.subplots(1, 3, figsize=(NATURE_W2, 2.5))

    ax = axes[0]
    ax.plot(sweep.K, sweep.code_knn_acc, "-o", color=c[0], ms=5, label="discrete code (one-hot)")
    ax.plot(sweep.K, sweep.quant_knn_acc, "-s", color=c[1], ms=5, label="quantised embedding")
    ax.axhline(cont_acc, ls="--", lw=1.2, color="0.35", label="continuous PCA (baseline)")
    ax.set_xlabel("Codebook size K"); ax.set_ylabel("kNN label accuracy (5-fold CV)")
    ax.set_title("Information kept after discretisation")
    ax.legend(loc="lower right")

    ax = axes[1]
    ax.plot(sweep.K, sweep.perplexity, "-o", color=c[2], ms=5, label="codebook perplexity")
    ax.plot(sweep.K, sweep.K, ls=":", lw=1.2, color="0.35", label="full usage (= K)")
    ax.set_xlabel("Codebook size K"); ax.set_ylabel("Perplexity  exp(H)")
    ax.set_title("Codebook usage / collapse")
    ax.legend(loc="upper left")

    ax = axes[2]
    ax.plot(sweep.K, sweep.recon_r2, "-o", color=c[3], ms=5)
    ax.set_xlabel("Codebook size K"); ax.set_ylabel(r"Quantisation $R^2$")
    ax.set_title("Reconstruction of the embedding")

    fig.tight_layout()
    save_fig(fig, os.path.join(outdir_assets, "572_codebook_sweep"))
    plt.close(fig)


def fig_code_usage(codes, labels, K, outdir_assets):
    """码使用频次 lollipop(按频次排序,颜色 = 该码的主导细胞类型)。"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    uniq = sorted(set(labels))
    cols = dict(zip(uniq, pal(len(uniq))))
    cnt = np.bincount(codes, minlength=K)
    dom = []
    for k in range(K):
        m = codes == k
        dom.append(max(uniq, key=lambda u: (labels[m] == u).sum()) if m.any() else uniq[0])
    order = np.argsort(-cnt)

    fig, ax = plt.subplots(figsize=(NATURE_W2, 2.8))
    xs = np.arange(K)
    ax.vlines(xs, 0, cnt[order], color="0.75", lw=1.2)
    ax.scatter(xs, cnt[order], s=42, zorder=3,
               color=[cols[dom[i]] for i in order], edgecolor="white", linewidth=0.6)
    ax.set_xticks(xs); ax.set_xticklabels([f"c{i}" for i in order], rotation=90, fontsize=6)
    ax.set_xlabel(f"Cell code (K = {K}, sorted by usage)"); ax.set_ylabel("Cells assigned")
    ax.set_title("Codebook usage, coloured by dominant cell type")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=cols[u], label=u) for u in uniq],
              loc="upper right", ncol=2)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir_assets, "572_code_usage_lollipop"))
    plt.close(fig)


def fig_code_heatmap(codes, labels, K, outdir_assets):
    """码 × 细胞类型列联热图(按码归一化)——这就是「码可解释性」的读出方式。"""
    import numpy as np, matplotlib.pyplot as plt
    uniq = sorted(set(labels))
    M = np.zeros((K, len(uniq)))
    for k in range(K):
        m = codes == k
        if m.any():
            M[k] = [(labels[m] == u).mean() for u in uniq]
    order = np.argsort([-M[k].argmax() * 1e3 - M[k].max() for k in range(K)])[::-1]
    M = M[order]

    fig, ax = plt.subplots(figsize=(NATURE_W1 + 0.6, 3.4))
    im = ax.imshow(M, aspect="auto", cmap=CMAP_CONT, vmin=0, vmax=1)
    ax.set_xticks(range(len(uniq))); ax.set_xticklabels(uniq, rotation=45, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{i}" for i in order], fontsize=6)
    ax.set_xlabel("Cell type"); ax.set_ylabel("Cell code")
    ax.set_title("Code composition")
    fig.colorbar(im, ax=ax, shrink=0.75, label="Fraction of cells in code")
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir_assets, "572_code_celltype_heatmap"))
    plt.close(fig)


def fig_embedding(Z, codes, labels, K, outdir_assets):
    """PC1-PC2 散点:左=真实类型,右=离散码(看词表把连续流形切成了什么)。"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    uniq = sorted(set(labels))
    cols = dict(zip(uniq, pal(len(uniq))))
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, 3.0), sharex=True, sharey=True)

    for u in uniq:
        m = labels == u
        axes[0].scatter(Z[m, 0], Z[m, 1], s=10, color=cols[u], label=u,
                        edgecolor="none", alpha=0.85)
    axes[0].set_title("Ground-truth cell type")
    axes[0].legend(loc="best", markerscale=1.6)

    ccol = pal(K, "vivid")
    axes[1].scatter(Z[:, 0], Z[:, 1], s=10, c=[ccol[i % len(ccol)] for i in codes],
                    edgecolor="none", alpha=0.85)
    axes[1].set_title(f"Discrete cell code (K = {K})")
    for ax in axes:
        ax.set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir_assets, "572_embedding_vs_code"))
    plt.close(fig)


# ------------------------------------------------------------------ CellVQ 守卫式封装
def cellvq_probe(repo: str | None, ckpt: str | None):
    """只做体检 + 打印已核实的官方命令,绝不代替官方实现。

    官方 inference.py 的真实参数(逐字读自 master 分支):
      --input_type {singlecell,bulk}  --pool_type {all,max}  --data_path  --save_path
      --pre_normalized  --mode (默认 m1)  --device  --model_path (默认 ./model/models/models.ckpt)
    官方 Python 入口(model/load.py:124):
      load_model_frommmf(best_ckpt_path, mode='m1', params=None, device='cuda')
        → (model.to(device), config)
    再 (model/pretrainmodels/model.py:156):
      get_cellcode(x, padding_label, encoder_position_gene_ids, output_attentions=False, **kwargs)
        → (geneemb, cell_code, encoder_position_gene_ids[:, indexes2[0]])
      inference.py 以关键字形式调用并解包为 (x, cell_code, _)。
    输出(inference.py:74):torch.save([cellembs, cell_codes], save_path),含两个张量的 .pth。
    上游文档漂移(均已核对源码,以源码为准):
      a) README 说权重放 pretrained_models/(表格里又写 pretrained_model/checkpoint.pt),
         而 inference.py 默认 ./model/models/models.ckpt → 一律用 --model_path 显式指定。
      b) README 的参数表列了 --verbose,inference.py 里**没有这个参数**,传了会报错。
      c) README 的参数表没列 --input_type / --pool_type / --pre_normalized,但 inference.py 有。
    """
    out = {"repo": "https://github.com/A4Bio/CellVQ", "branch": "master"}
    try:
        import torch
        out["torch"] = torch.__version__
        out["cuda_available"] = bool(torch.cuda.is_available())
    except ImportError:
        out["torch"] = "not installed"
        out["cuda_available"] = False

    if not repo:
        out["status"] = "skipped"
        out["reason"] = ("未提供 --cellvq-repo。CellVQ 不是 pip 包,需 git clone "
                         "https://github.com/A4Bio/CellVQ;README 让你 conda create python=3.9.17 "
                         "后跑 ./install.sh(其内容仅两行:pip install torch==2.6.0 / "
                         "pip install scanpy einops cell-gears)")
        return out

    out["repo_path"] = repo
    need = ["inference.py", os.path.join("model", "load.py")]
    missing = [f for f in need if not os.path.exists(os.path.join(repo, f))]
    if missing:
        out["status"] = "skipped"
        out["reason"] = f"{repo} 下缺少 {missing};看起来不是 CellVQ 仓库根目录"
        return out

    ckpt = ckpt or os.path.join(repo, "model", "models", "models.ckpt")
    out["model_path"] = ckpt
    out["checkpoint_present"] = os.path.exists(ckpt)
    if not out["checkpoint_present"]:
        out["reason"] = ("预训练权重需从 ModelScope 下载(仓库 README 称 511M 参数,"
                         "论文摘要称 500M):https://modelscope.cn/models/wj1006/CellVQ/files")
    out["status"] = "ready" if (out["checkpoint_present"] and out["cuda_available"]) else "blocked"
    # 刻意不把本模块的 CSV 填进去:官方模型只吃对齐到 19264 基因面板的 h5ad/npy,
    # 填进去会诱导使用者跑一条注定报错(或更糟,静默错位)的命令
    out["verified_command"] = (
        'python inference.py --data_path <your_19264_aligned.h5ad> '
        f'--save_path <out.pth> --mode m1 --device cuda --model_path {ckpt}'
    )
    out["input_contract"] = ("官方示例为 examples/cluster_19264.h5ad —— 19264 基因固定面板、"
                             "顺序敏感的 h5ad(或 .npy);面板本身是仓库内的 "
                             "model/OS_scRNA_gene_index.19264.tsv;本模块的合成 CSV 不满足该契约,"
                             "跑官方模型前须先按 preprocess/ 把自己的数据对齐到该基因面板")
    out["output_contract"] = "torch.save([cellembs, cell_codes], save_path) → .pth,两个张量"
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description="572 CellVQ 离散细胞词表 · 基线 + 守卫式封装")
    ap.add_argument("--input", default=EXAMPLE, help="cells x genes 计数 CSV(含 cell_type 列)")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--k-grid", default="4,8,12,16,24,32", help="扫描的码本大小,逗号分隔")
    ap.add_argument("--k-show", type=int, default=16, help="用于展示图的码本大小")
    ap.add_argument("--cellvq-repo", default=None, help="已 clone 的 A4Bio/CellVQ 仓库根目录")
    ap.add_argument("--cellvq-ckpt", default=None, help="CellVQ 权重路径(默认 <repo>/model/models/models.ckpt)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True); os.makedirs(ASSETS, exist_ok=True)
    set_pub_style(base_size=9)
    k_grid = [int(x) for x in a.k_grid.split(",") if x.strip()]

    print("[572] Step 1 读入计数并做 CP10K+log1p → PCA 连续嵌入")
    X, labels, genes = load_counts(a.input)
    Z = embed(X, a.n_pcs)
    print(f"       cells={X.shape[0]} genes={X.shape[1]} types={len(set(labels))} PCs={Z.shape[1]}")

    print("[572] Step 2 基线:k-means 码本扫描(VQ 的朴素下界)")
    sweep, cont_acc, per_k = run_baseline(Z, labels, k_grid, a.outdir)
    print(f"       continuous PCA kNN accuracy = {cont_acc:.3f}  ← 离散化要对标的天花板")
    print(sweep.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    K = a.k_show if a.k_show in per_k else k_grid[-1]
    codes = per_k[K]

    print("[572] Step 3 出图(lollipop / heatmap / 折线 / 散点,无条形图)")
    fig_sweep(sweep, cont_acc, ASSETS)
    fig_code_usage(codes, labels, K, ASSETS)
    fig_code_heatmap(codes, labels, K, ASSETS)
    fig_embedding(Z, codes, labels, K, ASSETS)

    print("[572] Step 4 CellVQ 守卫式体检")
    probe = cellvq_probe(a.cellvq_repo, a.cellvq_ckpt)
    for k, v in probe.items():
        print(f"       {k}: {v}")

    import pandas as pd
    pd.DataFrame({"cell_id": range(len(codes)), "cell_type": labels,
                  f"code_K{K}": codes}).to_csv(
        os.path.join(a.outdir, "572_cell_codes.csv"), index=False)
    summary = {
        "n_cells": int(X.shape[0]), "n_genes": int(X.shape[1]),
        "n_pcs": int(Z.shape[1]), "seed": SEED,
        "continuous_pca_knn_acc": round(cont_acc, 4),
        "k_grid": k_grid, "k_shown": K,
        "sweep": sweep.round(4).to_dict(orient="records"),
        "cellvq": probe,
        "session": save_session(),
        "upstream": {
            "paper": "Nat Commun 2026;17:4037",
            "doi": "10.1038/s41467-026-70071-5",
            "pmid": "41839876", "pmc": "PMC13139411",
            "repo": "https://github.com/A4Bio/CellVQ",
        },
    }
    with open(os.path.join(a.outdir, "572_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[572] 完成 → {a.outdir}  图 → {ASSETS}")


if __name__ == "__main__":
    main()
