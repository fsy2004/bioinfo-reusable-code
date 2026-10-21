#!/usr/bin/env python3
"""Persist portable, dated acquisition evidence after checking the local source cache.

No source code is redistributed or executed. Existing license conflicts are kept.
"""
from datetime import datetime, timezone
import hashlib
import importlib.util
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / 'modules/04_causal_inference_genetics/09_method_library/611_mr_coloc_library/method_library.py'
spec = importlib.util.spec_from_file_location('mr_method_library', MODULE)
lib = importlib.util.module_from_spec(spec); spec.loader.exec_module(lib)


def main():
    data = lib.load_registry()
    reports = []
    for source in data['sources']:
        receipt_path = lib.CACHE / source['id'] / 'receipt.json'
        if not receipt_path.is_file():
            reports.append({'source_id': source['id'], 'status': 'not_acquired'})
            continue
        integrity = lib.verify_source(source)
        if integrity['status'] != 'verified':
            raise ValueError(f"Unverified source cache: {source['id']}: {integrity}")
        receipt = lib.read_json(receipt_path)
        portable = {k:v for k,v in receipt.items() if k != 'previous_archive_attempt'}
        for key in ['source_tree', 'archive_file']:
            if portable.get(key):
                path = Path(portable[key])
                if path.is_absolute():
                    if not path.is_relative_to(ROOT):
                        raise ValueError('Source cache must be within repository for portable receipt')
                    portable[key] = path.relative_to(ROOT).as_posix()
        source['pinned_commit'] = receipt['commit']
        source['acquisition_snapshot'] = portable
        source['acquisition_snapshot_scope'] = 'Dated evidence on acquisition workstation; use verify-sources to check current local cache. Not runtime or algorithm-specific implementation proof.'
        manifest_path = receipt_path.parent / receipt['file_manifest']
        source['source_file_manifest_sha256'] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
        description = Path(receipt['source_tree']) / 'DESCRIPTION'
        if description.is_file():
            fields = {}
            for line in description.read_text(encoding='utf-8', errors='replace').splitlines():
                if line.startswith(('Version:', 'License:')):
                    key, value = line.split(':', 1); fields[key] = value.strip()
            if fields.get('Version'): source.setdefault('declared_version', fields['Version'])
            if fields.get('License'): source.setdefault('declared_license', fields['License'])
        if source.get('license_status') == 'pending_receipt_review':
            source['license_status'] = 'api_or_package_declaration_recorded_not_legal_resolution' if source.get('declared_license') or portable.get('license_api') not in {None, 'NOASSERTION'} else 'license_unresolved_inspect_evidence'
        reports.append(integrity)
    lib.write_json(lib.REGISTRY, data)
    report = {'as_of':datetime.now(timezone.utc).isoformat(), 'sources':len(reports),
              'verified_sources':sum(r['status']=='verified' for r in reports),
              'checked_files':sum(r.get('checked_files',0) for r in reports),
              'scientific_analysis':False,'reports':reports}
    stamp = datetime.now(timezone.utc).strftime('%Y%m%d')
    lib.write_json(ROOT / f'docs/methods/mr_coloc/source_integrity_{stamp}.json', report)
    print(json.dumps({k:v for k,v in report.items() if k!='reports'}, indent=2))


if __name__=='__main__':
    main()
