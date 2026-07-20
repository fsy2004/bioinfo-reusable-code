# -*- coding: utf-8 -*-
"""按「域 / 子类」两级结构重新生成 modules/CATALOG.md 与各域 README.md。

条目元信息(用途、输入输出、依赖、语言、图型、状态)来自:
  ① 上一版 CATALOG.md(git 历史里的 BASE_REV,老条目沿用,不手工重打)
  ② new_modules.json(本次新增模块,由建模块的工作流产出)
两者都没有的条目,标 "—" 并在末尾列成待补清单 —— 不编造。

用法: python modules/_framework/gen_catalog.py [--base <git-rev>]
"""
from __future__ import annotations
import argparse, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODULES = os.path.dirname(HERE)
REPO = os.path.dirname(MODULES)
BASE_REV = "325310ea"          # 重组前最后一次提交,老条目元信息从这里取

DOMAINS = {
    "01_single_cell": ("单细胞分析", "Single-cell analysis", {
        "01_pipeline_qc": ("上游与质控", "Pipeline & QC"),
        "02_integration_batch": ("整合与批次校正", "Integration & batch correction"),
        "03_annotation_typing": ("注释与细胞分型", "Annotation & cell typing"),
        "04_composition_da": ("组成与丰度差异", "Composition / differential abundance"),
        "05_differential_expression": ("差异表达(含 pseudobulk)", "Differential expression"),
        "06_trajectory_velocity": ("轨迹与 RNA 速率", "Trajectory & RNA velocity"),
        "07_cnv_clonality": ("拷贝数与克隆", "CNV & clonality"),
        "08_activity_scoring": ("通路/转录因子活性打分", "Pathway & TF activity"),
        "09_bulk_phenotype_link": ("单细胞↔bulk 表型关联", "Single-cell to bulk phenotype"),
        "10_foundation_models": ("单细胞基础模型", "Foundation models"),
    }),
    "02_spatial_transcriptomics": ("空间转录组", "Spatial transcriptomics", {
        "01_pipeline_segmentation": ("上游与细胞分割", "Pipeline & segmentation"),
        "02_domains_svg_stats": ("空间域、空间可变基因与空间统计", "Domains, SVG & spatial statistics"),
        "03_deconvolution_mapping": ("解卷积与单细胞映射", "Deconvolution & mapping"),
        "04_alignment_3d": ("切片配准与三维重建", "Slice alignment & 3D"),
        "05_cell_communication": ("细胞通讯", "Cell-cell communication"),
        "06_spatial_multiomics": ("空间多组学", "Spatial multi-omics"),
        "07_foundation_models": ("空间基础模型", "Spatial foundation models"),
    }),
    "03_virtual_perturbation": ("虚拟扰动技术", "Virtual perturbation", {
        "01_insilico_knockout": ("虚拟敲除与扰动模拟", "In-silico knockout"),
        "02_grn_inference": ("基因调控网络推断", "GRN inference"),
        "03_causal_perturbation": ("因果表示与反事实", "Causal & counterfactual"),
        "04_drug_perturbation": ("药物扰动与响应", "Drug perturbation"),
        "05_benchmark": ("扰动预测基准", "Perturbation benchmarks"),
    }),
    "04_causal_inference_genetics": ("因果推断与遗传流行病", "Causal inference & genetics", {
        "01_instrument_prep": ("工具变量准备", "Instrument preparation"),
        "02_two_sample_mr": ("两样本孟德尔随机化", "Two-sample MR"),
        "03_cis_mr_drug_target": ("cis-MR 与药靶", "cis-MR & drug targets"),
        "04_mediation_mvmr": ("中介与多变量 MR", "Mediation & MVMR"),
        "05_colocalization": ("共定位", "Colocalization"),
        "06_twas_sceqtl": ("TWAS 与单细胞 eQTL", "TWAS & sc-eQTL"),
        "07_robust_mr_methods": ("稳健 MR 估计量", "Robust MR estimators"),
    }),
    "05_machine_learning": ("机器学习", "Machine learning", {
        "01_feature_selection": ("特征筛选", "Feature selection"),
        "02_classification_models": ("分类模型", "Classification models"),
        "03_survival_ml": ("生存机器学习", "Survival ML"),
        "04_interpretability": ("可解释性", "Interpretability"),
        "05_uncertainty": ("不确定性量化", "Uncertainty quantification"),
        "06_generalization_validation": ("泛化与外部验证", "Generalization & validation"),
    }),
    "06_bulk_omics": ("Bulk 组学", "Bulk omics", {
        "01_differential_expression": ("差异表达", "Differential expression"),
        "02_enrichment": ("富集分析", "Enrichment"),
        "03_coexpression_networks": ("共表达网络(WGCNA 家族)", "Co-expression networks"),
        "04_multiomics_integration": ("多组学整合与分型", "Multi-omics integration"),
        "05_mutation_methylation_proteome": ("突变/甲基化/蛋白/代谢", "Mutation, methylation, proteome, metabolome"),
    }),
    "07_clinical_translational": ("临床与转化", "Clinical & translational", {
        "01_diagnostic_models": ("诊断模型", "Diagnostic models"),
        "02_prognosis_survival": ("预后与生存", "Prognosis & survival"),
        "03_immune_infiltration": ("免疫浸润与解卷积", "Immune infiltration & deconvolution"),
        "04_pharmacovigilance": ("药物警戒", "Pharmacovigilance"),
        "05_epidemiology_burden": ("疾病负担与人群队列", "Disease burden & population cohorts"),
    }),
    "08_structure_drug_design": ("结构与药物设计", "Structure & drug design", {
        "01_docking": ("分子对接", "Molecular docking"),
        "02_md_simulation": ("分子动力学", "Molecular dynamics"),
        "03_virtual_screening": ("虚拟筛选与打分", "Virtual screening & scoring"),
    }),
    "09_network_pharmacology": ("网络药理学", "Network pharmacology", {
        "01_target_databases": ("靶点数据库提取", "Target databases"),
        "02_target_intersection": ("靶点交集与集合图", "Target intersection"),
        "03_druggability": ("成药性评分", "Druggability"),
    }),
    "10_visualization": ("可视化", "Visualization", {
        "01_advanced_plots": ("高级图型", "Advanced plot types"),
        "02_templates_resources": ("模板与外部资源", "Templates & external resources"),
    }),
}

# 没有编号、因而不在历史 CATALOG 里的条目:元信息按实际内容手工登记一次。
EXTRA = {
    "01_single_cell/06_trajectory_velocity/491_sctour_extra_files": dict(
        status="📄", purpose="scTour 官方教程的复现脚本与环境记录(062 的配套材料)",
        io="教程数据 → 复现结果 + 环境说明", deps="Py · sctour", lang="Python/PS", figs="—"),
    "07_clinical_translational/05_epidemiology_burden/comorbidity_paper_template_refs.ris": dict(
        status="🗃️", purpose="共病选题的参考文献(仅本地)", io="—", deps="—", lang="—", figs="—"),
    "07_clinical_translational/05_epidemiology_burden/literature_summary_comorbidity.md": dict(
        status="🗃️", purpose="共病文献综述草稿(仅本地)", io="—", deps="—", lang="—", figs="—"),
    "07_clinical_translational/05_epidemiology_burden/sources_index.csv": dict(
        status="🗃️", purpose="疾病负担数据源索引(仅本地)", io="—", deps="—", lang="—", figs="—"),
    "07_clinical_translational/05_epidemiology_burden/topic_candidates.md": dict(
        status="🗃️", purpose="疾病负担选题候选(仅本地)", io="—", deps="—", lang="—", figs="—"),
    "07_clinical_translational/05_epidemiology_burden/99_external_sources": dict(
        status="🗃️", purpose="GBD/NHANES/CHARLS 上游第三方源码树(git 忽略,仅本地参考)",
        io="—", deps="—", lang="R", figs="—"),
    "10_visualization/02_templates_resources/templates": dict(
        status="📄", purpose="出图模板", io="—", deps="R · ggplot2", lang="R", figs="—"),
    "10_visualization/02_templates_resources/ai_scientific_figures": dict(
        status="🗃️", purpose="AutoFigure-Edit(ICLR'26)本地参考:方法描述 → 可编辑 SVG 示意图",
        io="—", deps="—", lang="Python", figs="schematic"),
    "10_visualization/02_templates_resources/advanced_figure_tools.csv": dict(
        status="📄", purpose="高级图型工具清单", io="—", deps="—", lang="—", figs="—"),
    "10_visualization/02_templates_resources/download_advanced_figure_tools.ps1": dict(
        status="📄", purpose="按清单批量拉取高级图型工具", io="清单 csv → 本地仓库", deps="PowerShell", lang="PS", figs="—"),
    "10_visualization/02_templates_resources/literature_download_links_for_fdm.txt": dict(
        status="📄", purpose="高级图型的文献下载链接", io="—", deps="—", lang="—", figs="—"),
}

HEAD = """# Module catalog

按 **域 → 子类** 两级组织。想做某类分析时,先定位域,再在子类里挑模块;
每行给出用途、输入→输出、依赖、语言与产出图型。新项目脚手架与统一出图样式见
[`_framework/`](_framework/)。

> **复用优先,绝不从头写。** 先在本目录挑模块(或挑一个真实已发表工具)再适配,
> 不要凭记忆手写分析代码 —— 那会带来假 API 和错参数。见
> [`_framework/CONVENTIONS.md` §0](_framework/CONVENTIONS.md)。

## 状态图例

| 标记 | 含义 |
|------|------|
| ✅ | 开箱即跑 —— 用自带合成示例数据本机跑通,零改动 |
| 🟡 | 核心复现或诚实基线本机可跑;完整方法需在分析服务器装包(见 [`_framework/SERVER_DEPENDENCIES.md`](_framework/SERVER_DEPENDENCIES.md)) |
| 🔴 | 重型 / GPU / 外部工具链 —— 守卫式引用封装,不在本机渲染 |
| 📄 | 模板或上游脚本 —— 自带数据 + 自行安装,无捆绑示例 |
| 📦 | Vendored 第三方包 —— 仅保留清单 / 本地 |
| 🗃️ | 仅本地,git 忽略 —— 不在公开仓库中 |

---

"""


def parse_old(rev: str) -> dict:
    out = subprocess.run(["git", "show", f"{rev}:modules/CATALOG.md"],
                         cwd=REPO, capture_output=True, text=True, encoding="utf-8").stdout
    pat = re.compile(r"^\|\s*([^|]*?)\s*\|\s*(\d{3})\s*\|\s*\[([^\]]+)\]\(([^)]+)\)\s*\|"
                     r"\s*([^|]*?)\s*\|\s*([^|]*?)\s*\|\s*([^|]*?)\s*\|\s*([^|]*?)\s*\|\s*([^|]*?)\s*\|\s*$", re.M)
    meta = {}
    for st, i, name, _href, purpose, io, deps, lang, figs in pat.findall(out):
        meta[i] = dict(status=st, purpose=purpose, io=io, deps=deps, lang=lang, figs=figs)
    return meta


def scan() -> dict:
    """{domain: {sub: [(id, entry_name, relpath), ...]}}"""
    tree = {}
    for dom in sorted(os.listdir(MODULES)):
        dp = os.path.join(MODULES, dom)
        if not os.path.isdir(dp) or dom.startswith("_"):
            continue
        tree[dom] = {}
        for sub in sorted(os.listdir(dp)):
            sp = os.path.join(dp, sub)
            if not os.path.isdir(sp):
                continue
            items = []
            for e in sorted(os.listdir(sp)):
                if e in ("README.md", ".gitignore"):
                    continue
                m = re.match(r"^(\d{3})", e)
                ep = os.path.join(sp, e)
                # 无编号的分组目录(如 01_GBD)里若还嵌着编号模块,展开到它们本身,
                # 否则模块会被一个不透明的文件夹名挡住。
                if not m and os.path.isdir(ep):
                    nested = [c for c in sorted(os.listdir(ep)) if re.match(r"^\d{3}", c)]
                    if nested:
                        for c in nested:
                            items.append((c[:3], f"{e}/{c}", f"{dom}/{sub}/{e}/{c}"))
                        continue
                items.append((m.group(1) if m else "", e, f"{dom}/{sub}/{e}"))
            if items:
                tree[dom][sub] = items
    return tree


def row(i, name, rel, meta, new):
    d = new.get(i) or meta.get(i) or EXTRA.get(rel) or {}
    st = d.get("status") or ("📄" if "." in name else "—")
    cells = [st, i or "—", f"[{name}]({rel})",
             d.get("purpose") or "—", d.get("io") or "—",
             d.get("deps") or "—", d.get("lang") or "—", d.get("figs") or "—"]
    return "| " + " | ".join(cells) + " |"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=BASE_REV)
    a = ap.parse_args()

    meta = parse_old(a.base)
    njson = os.path.join(HERE, "new_modules.json")
    new = {}
    if os.path.exists(njson):
        for m in json.load(open(njson, encoding="utf-8")):
            new[m["id"]] = dict(status=m.get("status", "—"), purpose=m.get("summary", "—"),
                                io=m.get("io", "—"), deps=m.get("deps", "—"),
                                lang=m.get("lang", "—"), figs=m.get("figures", "—"))
    tree = scan()

    cat = [HEAD]
    missing, total = [], 0
    for dom, subs in tree.items():
        if dom not in DOMAINS:
            sys.exit(f"未知域 {dom} —— 先在 DOMAINS 里登记")
        cn, en, submap = DOMAINS[dom]
        n = sum(len(v) for v in subs.values())
        cat.append(f"## {dom[:2]} · {cn} — {en}  ({n})\n")
        dom_lines = [f"# {dom[:2]} · {cn} — {en}\n",
                     f"本域共 {n} 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。\n"]
        for sub, items in subs.items():
            scn, sen = submap.get(sub, (sub, ""))
            cat.append(f"### {dom[:2]}.{sub[:2]} · {scn} — {sen}\n")
            cat.append("| St | # | 模块 | 用途 | 输入 → 输出 | 依赖 | 语言 | 图型 |")
            cat.append("|----|---|------|------|------------|------|------|------|")
            dom_lines.append(f"\n## {scn} — {sen}\n")
            for i, name, rel in items:
                cat.append(row(i, name, rel, meta, new))
                d = new.get(i) or meta.get(i) or EXTRA.get(rel) or {}
                dom_lines.append(f"- [{name}]({sub}/{name}) — {d.get('purpose') or '(用途待补)'}")
                if not d:
                    missing.append(rel)
                total += 1
            cat.append("")
        open(os.path.join(MODULES, dom, "README.md"), "w", encoding="utf-8").write("\n".join(dom_lines) + "\n")

    # ---- 图类型 → 模块 反查表(从各模块的图型字段自动生成,不手工维护) ----
    ALIAS = [                      # (归一化名, 匹配的关键词)
        ("Volcano · 火山图", ["volcano"]),
        ("Heatmap · 热图", ["heatmap", "heat map", "tile"]),
        ("ROC / PR", ["roc", "pr curve", "auroc"]),
        ("Calibration / DCA / nomogram · 校准与决策曲线", ["calibration", "dca", "decision curve", "nomogram"]),
        ("Forest · 森林图", ["forest"]),
        ("KM / survival · 生存曲线", ["km", "survival", "kaplan"]),
        ("UMAP / tSNE · 降维嵌入", ["umap", "tsne", "t-sne", "embedding"]),
        ("Violin · 小提琴", ["violin"]),
        ("Raincloud · 云雨图", ["raincloud"]),
        ("Ridgeline · 山脊图", ["ridgeline", "ridge"]),
        ("Box · 箱线图", ["box"]),
        ("Dot / bubble · 点图气泡图", ["dot", "bubble"]),
        ("Lollipop · 棒棒糖图(条形图替代)", ["lollipop"]),
        ("Dumbbell / slopegraph · 哑铃图与斜率图", ["dumbbell", "slope"]),
        ("Venn / UpSet · 集合图", ["venn", "upset"]),
        ("PCA", ["pca"]),
        ("Scatter · 散点", ["scatter"]),
        ("Network · 网络图", ["network", "cnet", "emap", "graph"]),
        ("Chord / circos / alluvial · 弦图环图桑基", ["chord", "circos", "circle", "alluvial", "sankey"]),
        ("Trajectory / vector field · 轨迹与向量场", ["trajectory", "vector", "pseudotime", "velocity"]),
        ("Spatial map · 空间分布图", ["spatial", "domain map", "scatterpie", "niche"]),
        ("Feature map · 基因表达投影", ["feature-map", "feature map", "featureplot"]),
        ("Composite multi-panel · 多面板拼图", ["composite"]),
    ]
    figidx = {k: set() for k, _ in ALIAS}
    for dom, subs in tree.items():
        for sub, items in subs.items():
            for i, name, rel in items:
                d = new.get(i) or meta.get(i) or EXTRA.get(rel) or {}
                f = (d.get("figs") or "").lower()
                if not f or f == "—":
                    continue
                for k, kws in ALIAS:
                    if any(w in f for w in kws):
                        figidx[k].add(i or name)
    cat.append("---\n\n## 图类型 → 模块 反查表\n")
    cat.append("想要某种图,直接查它由哪些模块产出。由各模块的「图型」字段自动生成。\n")
    for k, _ in ALIAS:
        ids = sorted(figidx[k])
        if ids:
            cat.append(f"- **{k}** → {', '.join(ids)}")
    cat.append("")

    if missing:
        cat.append("---\n\n## 元信息待补\n")
        cat.append("以下条目在目录中存在,但没有可靠的用途/输入输出记录 —— 用到时补,不臆造:\n")
        for m in missing:
            cat.append(f"- `{m}`")
        cat.append("")
    open(os.path.join(MODULES, "CATALOG.md"), "w", encoding="utf-8").write("\n".join(cat))
    print(f"CATALOG.md 已生成: {len(tree)} 域 · {sum(len(v) for v in tree.values())} 子类 · {total} 条目"
          f" · 元信息待补 {len(missing)}")


if __name__ == "__main__":
    main()
