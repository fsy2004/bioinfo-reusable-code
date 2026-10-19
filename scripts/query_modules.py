#!/usr/bin/env python3
"""Query stable module IDs, domain labels and method keywords without loading code."""
import argparse
import json
import sys
from pathlib import Path
if hasattr(sys.stdout,'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

ROOT=Path(__file__).resolve().parents[1]


def query(records, terms, domain=None):
    for row in records:
        if domain and not row['domain'].startswith(domain):
            continue
        haystack=json.dumps({k:row.get(k) for k in ('id','name','path','summary','tags','domain_labels','subcategory_labels','aliases')},ensure_ascii=False).casefold()
        if all(term.casefold() in haystack for term in terms):
            yield row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('terms',nargs='*')
    p.add_argument('--domain')
    p.add_argument('--json',action='store_true')
    args=p.parse_args()
    records=json.loads((ROOT/'catalog/modules.json').read_text(encoding='utf-8'))['modules']
    matches=list(query(records,args.terms,args.domain))
    if args.json:
        print(json.dumps(matches,ensure_ascii=False,indent=2))
    else:
        for r in matches:
            print(f"{r['id']}\t{r['path']}\t{r['summary'] or 'metadata gap'}")
        print(f'{len(matches)} matches')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
