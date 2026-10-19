#!/usr/bin/env python3
"""Build a deterministic, local-only module index; never read research data."""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRAMEWORK = ROOT / 'modules' / '_framework'
GENERATED = ('modules.json', 'modules.jsonl', 'domains.json')


def visible_files(root: Path) -> list[str]:
    raw = subprocess.check_output(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=root)
    return sorted(set(p for p in raw.decode('utf-8').split('\0') if p and (root / p).is_file()))


def discover(root: Path, files: list[str]) -> dict[str, tuple[str, list[str]]]:
    groups: dict[str, tuple[str, list[str]]] = {}
    for filename in files:
        parts = filename.split('/')
        if len(parts) < 4 or parts[0] != 'modules' or not re.match(r'^\d{2}_', parts[1]):
            continue
        for i, part in enumerate(parts[3:], 3):
            match = re.match(r'^(\d{3})(?:_|[A-Za-z]|\.)', part)
            if not match:
                continue
            path = '/'.join(parts[:i+1])
            # Results and nested upstream files cannot create extra module identities.
            if path not in groups:
                groups[path] = (match[1], [])
            groups[path][1].append(filename)
            break
    return groups


def load_definitions():
    spec = importlib.util.spec_from_file_location('module_catalog', FRAMEWORK/'gen_catalog.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.DOMAINS


def build(root: Path = ROOT) -> dict[str, str]:
    files = visible_files(root)
    visible = set(files)
    groups = discover(root, files)
    definitions = load_definitions()
    metadata = {m['id']: m for m in json.loads((root/'modules/_framework/module_metadata.json').read_text(encoding='utf-8'))}
    additions = {m['id']: m for m in json.loads((root/'modules/_framework/new_modules.json').read_text(encoding='utf-8'))}
    alias_file = root/'catalog/path_aliases.json'
    aliases = json.loads(alias_file.read_text(encoding='utf-8')) if alias_file.exists() else {}
    records, seen = [], set()
    for path, (identifier, members) in sorted(groups.items(), key=lambda p: (p[1][0], p[0])):
        if identifier in seen:
            raise ValueError(f'duplicate module ID: {identifier}')
        seen.add(identifier)
        target = root/path
        contract_file = target/'module.json' if target.is_dir() else None
        contract = json.loads(contract_file.read_text(encoding='utf-8-sig')) if contract_file and f'{path}/module.json' in visible else {}
        previous, new = metadata.get(identifier, {}), additions.get(identifier, {})
        summary = contract.get('summary') or contract.get('purpose') or contract.get('title') or new.get('summary') or previous.get('purpose')
        entrypoints = [p for p in members if (root/p).suffix.lower() in ('.r','.py','.ps1','.sh') and (not target.is_dir() or (root/p).parent == target)]
        declared = contract.get('entrypoint')
        if declared:
            declared = [declared] if isinstance(declared,str) else declared
            entrypoints = [f'{path}/{entry}' for entry in declared]
            if any(entry not in visible for entry in entrypoints):
                raise ValueError(f'Declared entrypoint is not Git-visible: {path}')
        readme = f'{path}/README.md' if target.is_dir() and f'{path}/README.md' in visible else None
        dom, sub = path.split('/')[1:3]
        dom_labels = definitions.get(dom, (dom, dom, {}))
        sub_labels = dom_labels[2].get(sub, (sub, sub))
        tags = sorted(set([dom, sub] + contract.get('tags', []) + contract.get('task', [])))
        # Query keywords include domain names and the actual module description.
        record = dict(id=identifier, name=target.stem if target.is_file() else target.name,
                      path=path, domain=dom, domain_labels=list(dom_labels[:2]),
                      subcategory=sub, subcategory_labels=list(sub_labels),
                      summary=summary, readme=readme, entrypoints=entrypoints,
                      language=contract.get('language') or new.get('lang') or previous.get('lang'),
                      input_output=contract.get('io') or ({'inputs':contract.get('inputs',[]),'outputs':contract.get('outputs',[])} if contract else None) or new.get('io') or previous.get('io'),
                      dependencies=contract.get('dependencies') or new.get('deps') or previous.get('deps'),
                      tags=tags, legacy_status=new.get('status') or previous.get('status'),
                      verification=contract.get('verification') or contract.get('evidence') or 'historical_metadata_not_revalidated',
                      sources=contract.get('sources', []),
                      citation=new.get('citation'),
                      aliases=sorted(old for old,newpath in aliases.items() if newpath == path),
                      metadata_source=(f'{path}/module.json' if contract else 'modules/_framework/module_metadata.json + new_modules.json'),
                      file_count=len(members))
        if not entrypoints:
            record['entrypoint_status']='reference_or_companion_resource'
        if summary is None:
            record['metadata_gap']='summary_not_recorded'
        records.append(record)
    for old, newpath in aliases.items():
        if not (root/newpath).exists():
            raise ValueError(f'alias target missing: {old} -> {newpath}')
    domains=[]
    for dom in sorted({r['domain'] for r in records}):
        subset=[r for r in records if r['domain']==dom]
        domains.append(dict(id=dom, labels=subset[0]['domain_labels'], modules=len(subset),
                            subcategories=sorted({r['subcategory'] for r in subset})))
    return {'modules.json':json.dumps({'schema_version':1,'scope':'Git-visible numbered module entries; runtime availability is not inferred','modules':records},ensure_ascii=False,indent=2)+'\n',
            'modules.jsonl':''.join(json.dumps(r,ensure_ascii=False,separators=(',',':'))+'\n' for r in records),
            'domains.json':json.dumps(domains,ensure_ascii=False,indent=2)+'\n'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    dest=ROOT/'catalog'
    dest.mkdir(exist_ok=True)
    output=build()
    stale=[]
    for name,text in output.items():
        path=dest/name
        if not path.exists() or path.read_text(encoding='utf-8')!=text:
            stale.append(name)
            if not args.check:
                path.write_text(text,encoding='utf-8',newline='\n')
    if args.check and stale:
        print('Stale AI index: '+', '.join(stale))
        return 1
    count=len(json.loads(output['modules.json'])['modules'])
    print(f'AI index: {count} modules; '+('current' if args.check else 'generated'))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
