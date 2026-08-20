#!/usr/bin/env python3
"""Refresh generated README blocks and fail CI when they drift."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"


def replace_block(text: str, name: str, content: str) -> str:
    start = f"<!-- {name}:start -->"
    end = f"<!-- {name}:end -->"
    pattern = re.compile(re.escape(start) + r".*?" + re.escape(end), re.DOTALL)
    if not pattern.search(text):
        raise ValueError(f"README is missing generated block {name}")
    return pattern.sub(f"{start}\n{content.rstrip()}\n{end}", text, count=1)


def local_link_errors(root: Path, text: str) -> list[str]:
    missing = []
    for raw in re.findall(r"\]\(([^)]+)\)", text):
        target = raw.strip().strip("<>")
        if not target or target.startswith(("#", "http://", "https://", "mailto:")):
            continue
        target = target.split("#", 1)[0]
        if target and not (root / target).exists() and target not in missing:
            missing.append(target)
    return missing


def compute_metrics(root: Path) -> dict[str, int]:
    modules = root / "modules"
    domains = [path for path in modules.iterdir() if path.is_dir() and re.match(r"^\d{2}_", path.name)]
    subcategories = [path for domain in domains for path in domain.iterdir() if path.is_dir() and re.match(r"^\d{2}_", path.name)]
    units = [path for subcategory in subcategories for path in subcategory.iterdir() if re.match(r"^\d{3}", path.name)]
    entries = [path for subcategory in subcategories for path in subcategory.iterdir() if path.name not in {"README.md", ".gitignore"}]
    return {"domains": len(domains), "subcategories": len(subcategories), "modules": len(units), "catalog_entries": len(entries)}


def render_metrics(metrics: dict[str, int]) -> str:
    return (
        f"![modules](https://img.shields.io/badge/modules-{metrics['modules']}-2563EB)\n"
        f"![catalog entries](https://img.shields.io/badge/catalog_entries-{metrics['catalog_entries']}-334155)\n"
        f"![domains](https://img.shields.io/badge/domains-{metrics['domains']}-0F766E)\n"
        f"![subcategories](https://img.shields.io/badge/subcategories-{metrics['subcategories']}-7C3AED)\n"
        "![outputs](https://img.shields.io/badge/figures-PDF%20%2B%20300dpi%20PNG-15803D)"
    )


def render_domains(root: Path) -> str:
    rows = ["| Domain | Entries | Entry |", "|---|---:|---|"]
    for domain in sorted(path for path in (root / "modules").iterdir() if path.is_dir() and re.match(r"^\d{2}_", path.name)):
        first = (domain / "README.md").read_text(encoding="utf-8").splitlines()[0]
        label = first.removeprefix("# ").strip()
        count = sum(1 for path in domain.glob("*/*") if path.name not in {"README.md", ".gitignore"})
        rows.append(f"| {label} | {count} | [`modules/{domain.name}`](modules/{domain.name}) |")
    return "\n".join(rows)


def updated_text(root: Path, source: str) -> str:
    source = replace_block(source, "readme-metrics", render_metrics(compute_metrics(root)))
    return replace_block(source, "readme-domains", render_domains(root))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="fail instead of writing when README is stale")
    args = parser.parse_args()
    source = README.read_text(encoding="utf-8")
    expected = updated_text(ROOT, source)
    missing = local_link_errors(ROOT, expected)
    if missing:
        print("README has missing local links: " + ", ".join(missing), file=sys.stderr)
        return 1
    if expected == source:
        print("README generated blocks are current")
        return 0
    if args.check:
        print("README generated blocks are stale; run: python scripts/update_readme.py", file=sys.stderr)
        return 1
    README.write_text(expected, encoding="utf-8", newline="\n")
    print("README generated blocks updated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
