#!/usr/bin/env python3
"""Render the curated methods registry and validate all local routing references."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DOCS=ROOT/'docs/research/2026_glucocorticoids_cancer'


def render():
    methods=json.loads((DOCS/'methods.json').read_text(encoding='utf-8'))['methods']
    sources=json.loads((DOCS/'github_sources.json').read_text(encoding='utf-8'))['sources']
    modules={r['id']:r for r in json.loads((ROOT/'catalog/modules.json').read_text(encoding='utf-8'))['modules']}
    papers={p['id'] for f in ['observational_papers.json','mr_gwas_papers.json'] for p in json.loads((DOCS/f).read_text(encoding='utf-8'))['papers']}
    source_ids={s['repository'] for s in sources}
    seen=set()
    for row in methods:
        if row['id'] in seen:raise ValueError('Duplicate method ID: '+row['id'])
        seen.add(row['id'])
        if any(p not in papers for p in row['paper_ids']):raise ValueError('Unknown paper in '+row['id'])
        if any(s['repository'] not in source_ids for s in row['github_sources']):raise ValueError('Unknown source in '+row['id'])
        if any(m not in modules for m in row['related_module_ids']):raise ValueError('Unknown module in '+row['id'])
    extracted={m for f in ['observational_papers.json','mr_gwas_papers.json'] for p in json.loads((DOCS/f).read_text(encoding='utf-8'))['papers'] for m in p.get('methods',[])}
    if extracted-seen:raise ValueError('Unmapped extracted methods: '+', '.join(sorted(extracted-seen)))
    lines=['# 方法 → 论文 → GitHub → 模块', '',
           '范围：纳入论文中逐项提取的方法，加上明确标为扩展的设计环节。每项资源都有实现状态；相关模块只是复用线索，不能据此声称该方法已完整实现。', '',
           '七个新入口 604–610 已完成合成输入行为核查，均需显式输入/配置；其余旧模块保留历史状态，仅参考资源没有安装、运行或被包装成空壳模块。原始文献方法不构成质量背书。', '',
           '机器读取用 [methods.json](methods.json)，来源版本/许可用 [github_sources.json](github_sources.json)。论文 ID 到 DOI、读取级别和参数边界分别见 [观察性文献](observational_papers.json) 与 [MR/GWAS 文献](mr_gwas_papers.json)。', '',
           '状态：`new_adapter_synthetic_checked`=本次新增且合成接口已核查；`reference_only_with_related_entries`=已核来源与相关旧入口，具体方法仍需适配；`reference_only`=仅来源登记；其余状态按行解释。', '']
    for family in sorted({r['family'] for r in methods}):
        lines += ['## '+family,'','| 方法 ID | 论文 ID | GitHub 来源 | 相关模块 | 状态与边界 |','|---|---|---|---|---|']
        for row in methods:
            if row['family']!=family:continue
            github=', '.join(f"[{s['repository']}]({s['url']})" for s in row['github_sources']) or '未确定独立 GitHub；见论文/设计说明'
            module_links=', '.join(f"[{i}](../../../{modules[i]['path']})" for i in row['related_module_ids']) or '—'
            notes=' '.join(row['notes']).replace('|','\\|')
            lines.append('| '+ ' | '.join([f"`{row['id']}`",', '.join(row['paper_ids']) or '设计扩展',github,module_links,f"`{row['implementation_status']}`；{notes}"])+' |')
        lines.append('')
    source_lines=['# GitHub 来源核查清单','','这是代码来源与复用路线清单，不是依赖锁文件。HEAD 未获得的条目明确保留空值；当前源码版本和本机运行包版本分别记录。镜像、维护者 fork、迁移 stub 与有效许可待确认的仓库有单独备注。没有复制第三方 README 全文或上游代码树。','','| 仓库 | 来源 revision | 许可记录 | 状态/备注 |','|---|---|---|---|']
    for source in sorted(sources,key=lambda r:r['repository'].casefold()):
        license_text=source.get('license_declared') or source.get('license_from_pinned_source') or source.get('license') or source.get('license_readme_declaration') or ('已读许可文件标题：'+source['license_file_heading'] if source.get('license_file_heading') else '未确认')
        sha=source.get('head_sha')
        revision=f"[`{sha[:12]}`]({source['url']}/commit/{sha})" if sha else '未冻结；见 JSON 的元数据链接'
        note=source.get('reuse_policy') or source.get('source_note') or source.get('read_level') or source.get('provenance') or source.get('retrieval_note') or '只核查来源元数据；运行状态见对应模块'
        source_lines.append(f"| [{source['repository']}]({source['url']}) | {revision} | {str(license_text).replace('|',' / ')} | {str(note).replace('|',' / ')} |")
    source_lines.append('')
    return {'METHODS.md':'\n'.join(lines),'GITHUB_SOURCES.md':'\n'.join(source_lines)}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args()
    changed=[]
    for name,text in render().items():
        target=DOCS/name
        if not target.exists() or target.read_text(encoding='utf-8')!=text:
            changed.append(name)
            if not args.check:target.write_text(text,encoding='utf-8',newline='\n')
    if changed and args.check:
        print('Stale research inventory: '+', '.join(changed));return 1
    print('Research methods/source inventory valid'+(' and current' if args.check else ''))
    return 0


if __name__=='__main__':raise SystemExit(main())
