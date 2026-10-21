#!/usr/bin/env python3
"""Query, qualify and acquire pinned official MR/colocalization implementations.

Selection checks declared input compatibility, never IV validity or causality.
No downloaded source is executed; analysis adapters call explicit installed APIs.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from urllib.parse import quote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REGISTRY = HERE / "methods.json"
CACHE = ROOT / "modules/_upstream_sources/mr_coloc"
CATEGORIES = {"mr_estimator", "mr_design", "diagnostic", "colocalization", "fine_mapping", "preparation", "correction_framework", "related_inference"}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def load_registry(path=REGISTRY):
    data = read_json(path)
    validate_registry(data)
    return data


def validate_registry(data):
    ids = [row["id"] for row in data["methods"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate stable method IDs")
    sources = {row["id"] for row in data["sources"]}
    if len(sources) != len(data["sources"]):
        raise ValueError("Duplicate source IDs")
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+", s) for s in sources):
        raise ValueError("Unsafe source ID")
    for row in data["methods"]:
        if row["category"] not in CATEGORIES:
            raise ValueError(f"Unknown method category: {row['id']}")
        if not row.get("name") or not row.get("assumptions"):
            raise ValueError(f"Missing name/assumptions: {row['id']}")
        if not set(row.get("source_ids", [])).issubset(sources):
            raise ValueError(f"Unresolved upstream source: {row['id']}")
        if row.get("adapter") and not row.get("implementations"):
            raise ValueError(f"Adapter without native implementation: {row['id']}")
    return True


def query(data, terms="", category=None):
    tokens = terms.lower().split()
    return [r for r in data["methods"] if (not category or r["category"] == category)
            and all(t in json.dumps(r, ensure_ascii=False).lower() for t in tokens)]


def select_methods(data, context):
    """Return reasons and candidates; a compatible method is not a valid IV model."""
    n = context.get("n_snps")
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError("n_snps must be a positive integer")
    if context.get("ld", "none") not in {"none", "signed_r", "r2"}:
        raise ValueError("ld must be none, signed_r or r2")
    if context.get("data_level", "summary") not in {"summary", "individual"}:
        raise ValueError("data_level must be summary or individual")
    for name in ("independent_instruments", "dense_region", "multiple_signals"):
        if name in context and not isinstance(context[name], bool):
            raise ValueError(name + " must be a JSON boolean")
    traits = context.get("n_traits", 1)
    if isinstance(traits, bool) or not isinstance(traits, int) or traits < 1:
        raise ValueError("n_traits must be a positive integer")
    if context.get("category") and context["category"] not in CATEGORIES:
        raise ValueError("Unknown selection category")
    accepted, blocked = [], []
    for row in data["methods"]:
        if context.get("category") and row["category"] != context["category"]:
            continue
        req = row.get("requirements", {})
        reasons = []
        if not row.get("selectable", False):
            reasons.append("reference_or_unqualified_requirements; read original source")
        if req.get("min_snps") is None:
            reasons.append("minimum_input_not_curated")
        elif n < req["min_snps"]:
            reasons.append(f"requires_at_least_{req['min_snps']}_SNPs")
        if req.get("exact_snps") and n != req["exact_snps"]:
            reasons.append(f"requires_exactly_{req['exact_snps']}_SNPs")
        if req.get("individual_data") and context.get("data_level", "summary") != "individual":
            reasons.append("requires_individual_data")
        if req.get("signed_ld") and context.get("ld", "none") != "signed_r":
            reasons.append("requires_signed_r_not_r2")
        if req.get("independent_iv") and not context.get("independent_instruments", False):
            reasons.append("independent_instruments_not_declared")
        if req.get("dense_region") and not context.get("dense_region", False):
            reasons.append("requires_dense_region_not_clumped_or_significant_only")
        if req.get("min_traits", 1) > context.get("n_traits", 1):
            reasons.append("requires_more_traits")
        if context.get("multiple_signals") and req.get("single_causal_variant"):
            reasons.append("single_causal_variant_assumption_conflicts_with_declared_multiple_signals")
        result = {"id": row["id"], "name": row["name"], "category": row["category"],
                  "tier": row.get("tier", "specialized"), "adapter": row.get("adapter"), "reasons": reasons}
        (blocked if reasons else accepted).append(result)
    ranks = {"baseline": 0, "sensitivity": 1, "specialized": 2, "reference": 3}
    accepted.sort(key=lambda r: (ranks.get(r["tier"], 9), r["id"]))
    return {"context": context, "compatible_for_review": accepted, "blocked": blocked,
            "scientific_validity_checked": False,
            "rule": "Input compatibility is not instrument validity, pleiotropy exclusion, biological relevance or causal proof. Do not run all compatible methods as votes."}


def api(path):
    if shutil.which("gh"):
        proc = subprocess.run(["gh", "api", path], capture_output=True, text=True, encoding="utf-8")
        if proc.returncode:
            raise RuntimeError("GitHub API request failed; status not credentials recorded")
        return json.loads(proc.stdout)
    request = urllib.request.Request("https://api.github.com/" + path, headers={"User-Agent": "mr-coloc-method-library"})
    with urllib.request.urlopen(request, timeout=45) as response:
        return json.load(response)


def safe_extract(archive, target, max_bytes=150_000_000):
    target = Path(target).resolve()
    total = 0
    files = []
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            parts = PurePosixPath(member.name).parts
            if not parts or any(p in {"..", ".git"} or ":" in p or "\\" in p for p in parts) or member.name.startswith("/"):
                raise ValueError("Unsafe archive path")
            if member.issym() or member.islnk():
                continue  # Source code is inspected as regular files; links never followed.
            if not member.isfile():
                continue
            relative = Path(*parts[1:])
            if not relative.parts:
                continue
            dest = (target / relative).resolve()
            if not dest.is_relative_to(target):
                raise ValueError("Archive file escapes source cache")
            total += member.size
            if total > max_bytes or len(files) >= 20000:
                raise ValueError("Expanded source size/file bound exceeded")
            dest.parent.mkdir(parents=True, exist_ok=True)
            stream = tar.extractfile(member)
            if stream is None:
                raise ValueError("Unreadable regular archive member")
            with dest.open("wb") as handle:
                shutil.copyfileobj(stream, handle)
            files.append({"path": relative.as_posix(), "bytes": member.size,
                          "sha256": hashlib.sha256(dest.read_bytes()).hexdigest()})
    return files


def validate_cache_identity(source, receipt, manifest):
    requested = source.get("repository")
    recorded = {str(receipt.get(k, "")).lower() for k in ("repository", "requested_repository")}
    if receipt.get("source_id") and receipt["source_id"] != source["id"]:
        raise ValueError("Cached source ID does not match request")
    if requested and requested.lower() not in recorded:
        raise ValueError("Cached repository identity differs; preserve old version and inspect")
    if source.get("pinned_commit") and source["pinned_commit"] != receipt.get("commit"):
        raise ValueError("Cached commit differs from requested pin; preserve old version and inspect")
    if receipt.get("commit") != manifest.get("commit") or receipt.get("repository") != manifest.get("repository"):
        raise ValueError("Receipt and source manifest identity disagree")


def verify_source(source, cache=CACHE):
    """Verify every regular source file against its acquisition receipt, without execution."""
    source_dir = Path(cache) / source["id"]
    receipt_file = source_dir / "receipt.json"
    if not receipt_file.is_file():
        return {"source_id": source["id"], "status": "not_acquired"}
    receipt = read_json(receipt_file)
    if receipt.get("status") not in {"source_acquired", "source_subset_acquired"}:
        return {"source_id": source["id"], "status": receipt.get("status", "not_acquired")}
    tree = Path(receipt["source_tree"]).resolve()
    if not tree.is_relative_to(source_dir.resolve()):
        raise ValueError("Source tree escapes cache")
    failures = []
    archive = source_dir / "archive.tar.gz"
    if receipt.get("status") == "source_acquired" and (not archive.is_file() or hashlib.sha256(archive.read_bytes()).hexdigest() != receipt["archive_sha256"]):
        failures.append("archive_sha256")
    manifest = read_json(source_dir / receipt["file_manifest"])
    validate_cache_identity(source, receipt, manifest)
    for entry in manifest["files"]:
        path = (tree / entry["path"]).resolve()
        if not path.is_relative_to(tree):
            raise ValueError("Source manifest escapes source tree")
        if not path.is_file() or path.stat().st_size != entry["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            failures.append(entry["path"])
    expected = {entry["path"] for entry in manifest["files"]}
    actual = {p.relative_to(tree).as_posix() for p in tree.rglob("*") if p.is_file() or p.is_symlink()}
    failures.extend("unreceipted_file:" + p for p in sorted(actual - expected))
    return {"source_id": source["id"], "status": "verified" if not failures else "integrity_failed",
            "checked_files": len(manifest["files"]), "failures": failures}


def fetch_source(source, cache=CACHE, max_archive_bytes=30_000_000):
    identifier = source["id"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", identifier):
        raise ValueError("Unsafe source ID")
    repo = source.get("repository")
    if not repo or not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo):
        return {"source_id": identifier, "status": "official_distribution_or_paper_only_not_GitHub"}
    cache = Path(cache).resolve()
    source_dir = cache / identifier
    if not source_dir.resolve().is_relative_to(cache) or source_dir.resolve() == cache:
        raise ValueError("Source directory escapes cache")
    source_dir.mkdir(parents=True, exist_ok=True)
    receipt_file = source_dir / "receipt.json"
    if receipt_file.exists():
        previous = read_json(receipt_file)
        archive = source_dir / "archive.tar.gz"
        if previous.get("status") == "source_acquired":
            integrity = verify_source(source, cache)
            if integrity["status"] != "verified":
                raise ValueError("Existing source integrity failed; preserve and inspect")
            return {**previous, "network_bytes_received": 0, "acquisition_reused": True}
    receipt = {"source_id": identifier, "requested_repository": repo, "checked_on": datetime.now(timezone.utc).date().isoformat(),
               "downloaded_source_executed": False, "redistributed_in_Git": False}
    count = 0
    try:
        metadata = api("repos/" + repo)
        canonical = metadata["full_name"]
        pinned = source.get("pinned_commit") or metadata["default_branch"]
        if source.get("pinned_commit") and not re.fullmatch(r"[0-9a-f]{40}", pinned):
            raise ValueError("Pinned commit must be immutable 40-character SHA")
        head = api("repos/" + canonical + "/commits/" + pinned)
        commit = head["sha"]
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("Invalid immutable upstream commit")
        receipt.update(repository=canonical, default_branch=metadata["default_branch"], commit=commit,
                       official_url=metadata["html_url"], archived=metadata["archived"],
                       license_api=(metadata.get("license") or {}).get("spdx_id"),
                       commit_date=head["commit"]["committer"]["date"])
        url = f"https://codeload.github.com/{canonical}/tar.gz/{commit}"
        archive = source_dir / "archive.tar.gz"
        partial = source_dir / "archive.part"
        request = urllib.request.Request(url, headers={"User-Agent": "mr-coloc-method-library"})
        with urllib.request.urlopen(request, timeout=100) as response:
            with partial.open("wb") as handle:
                while True:
                    chunk = response.read(65536)
                    if not chunk:
                        break
                    count += len(chunk)
                    if count > max_archive_bytes:
                        raise ValueError("Archive exceeds bounded code acquisition size")
                    handle.write(chunk)
        partial.replace(archive)
        receipt.update(archive_url=url, archive_bytes=archive.stat().st_size,
                       archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                       archive_file=archive.relative_to(ROOT).as_posix() if archive.is_relative_to(ROOT) else str(archive))
        extracted = source_dir / commit
        if extracted.exists() and any(extracted.iterdir()):
            raise ValueError("Existing source tree lacks matching receipt; preserve and inspect")
        files = safe_extract(archive, extracted)
        code = [f for f in files if Path(f["path"]).suffix.lower() in {".r", ".py", ".cpp", ".c", ".h", ".hpp", ".sh", ".rs"}]
        evidence = [f["path"] for f in files if Path(f["path"]).name.upper().startswith(("LICENSE", "COPYING")) or Path(f["path"]).name == "DESCRIPTION"]
        receipt.update(status="source_acquired", source_tree=str(extracted), source_files=len(files),
                       code_files=len(code), license_evidence=evidence,
                       file_manifest="source_files.json", source_complete_regular_files_only=True)
        write_json(source_dir / "source_files.json", {"repository": canonical, "commit": commit, "files": files})
    except Exception as exc:
        receipt.update(status="acquisition_failed", error_type=type(exc).__name__, error=str(exc)[:220])
    receipt["network_bytes_received"] = count
    write_json(receipt_file, receipt)
    return receipt


def fetch_code_subset(source, cache=CACHE, max_bytes=10_000_000, max_files=1200):
    """Acquire bounded code/docs only when a full repository includes large data.

    This is explicitly a subset, never an archive or complete repository claim.
    Existing full snapshots are retained. Official source files are not executed.
    """
    identifier = source["id"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", identifier):
        raise ValueError("Unsafe source ID")
    repo = source.get("repository")
    if not repo or not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo):
        return {"source_id": identifier, "status": "not_GitHub"}
    source_dir = Path(cache).resolve() / identifier
    previous = read_json(source_dir / "receipt.json") if (source_dir / "receipt.json").is_file() else {}
    if previous.get("status") in {"source_acquired", "source_subset_acquired"}:
        check = verify_source(source, cache)
        if check["status"] != "verified":
            raise ValueError("Existing source integrity failed; preserve and inspect")
        return previous
    receipt = {"source_id": identifier, "requested_repository": repo,
               "checked_on": datetime.now(timezone.utc).date().isoformat(),
               "downloaded_source_executed": False, "redistributed_in_Git": False,
               "previous_archive_attempt": previous,
               "scope": "bounded code/license/readme/metadata subset; data and large files excluded"}
    try:
        metadata = api("repos/" + repo)
        canonical = metadata["full_name"]
        pinned = source.get("pinned_commit") or previous.get("commit") or metadata["default_branch"]
        commit = api("repos/" + canonical + "/commits/" + pinned)["sha"]
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("Invalid commit")
        tree = api("repos/" + canonical + "/git/trees/" + commit + "?recursive=1")
        if tree.get("truncated"):
            raise ValueError("Tree listing truncated; cannot establish bounded source coverage")
        selected, skipped, total = [], [], 0
        code_ext = {".r", ".py", ".cpp", ".c", ".h", ".hpp", ".sh", ".rs", ".jl", ".f", ".f90", ".m", ".ado", ".do", ".mata", ".stan"}
        for entry in sorted(tree["tree"], key=lambda r: r["path"]):
            path = PurePosixPath(entry["path"])
            size = entry.get("size", 0)
            wanted = path.suffix.lower() in code_ext or path.name.upper().startswith(("README", "LICENSE", "COPYING")) or path.name in {"DESCRIPTION", "NAMESPACE", "pyproject.toml", "setup.cfg", "requirements.txt", "CMakeLists.txt"}
            if entry["type"] != "blob" or entry.get("mode") == "120000" or not wanted:
                continue
            if any(p in {"..", ".git"} or ":" in p or "\\" in p for p in path.parts):
                raise ValueError("Unsafe repository source path")
            if size > 500_000 or total + size > max_bytes or len(selected) >= max_files:
                skipped.append(entry["path"]); continue
            total += size; selected.append(entry)
        target = source_dir / (commit + "-code")
        if target.exists() and any(target.iterdir()):
            raise ValueError("Existing unreceipted subset preserved; inspect before retry")
        target.mkdir(parents=True, exist_ok=True)
        def download(entry):
            url = f"https://raw.githubusercontent.com/{canonical}/{commit}/" + quote(entry["path"])
            req = urllib.request.Request(url, headers={"User-Agent": "mr-coloc-method-library"})
            with urllib.request.urlopen(req, timeout=45) as response:
                raw = response.read(500_001)
            if len(raw) != entry["size"] or len(raw) > 500_000:
                raise ValueError("Source size differs from pinned tree")
            # Git blob SHA independently checks immutable content identity.
            blob_hash = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
            if blob_hash != entry["sha"]:
                raise ValueError("Pinned Git blob hash mismatch")
            dest = target.joinpath(*PurePosixPath(entry["path"]).parts)
            dest.parent.mkdir(parents=True, exist_ok=True); dest.write_bytes(raw)
            return {"path": entry["path"], "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "git_blob_sha": blob_hash}
        with ThreadPoolExecutor(max_workers=6) as pool:
            files = list(pool.map(download, selected))
        code = [f for f in files if Path(f["path"]).suffix.lower() in code_ext]
        if not code:
            raise ValueError("No regular source code found in bounded official subset")
        receipt.update(status="source_subset_acquired", repository=canonical, commit=commit,
                       official_url=metadata["html_url"], license_api=(metadata.get("license") or {}).get("spdx_id"),
                       source_tree=str(target), source_files=len(files), code_files=len(code),
                       subset_bytes=sum(f["bytes"] for f in files), source_complete_regular_files_only=False,
                       license_evidence=[f["path"] for f in files if Path(f["path"]).name.upper().startswith(("LICENSE", "COPYING")) or Path(f["path"]).name == "DESCRIPTION"],
                       file_manifest="source_files.json", eligible_files_skipped_by_bound=skipped)
        write_json(source_dir / "source_files.json", {"repository": canonical, "commit": commit, "scope": receipt["scope"], "files": files})
    except Exception as exc:
        receipt.update(status="acquisition_failed", error_type=type(exc).__name__, error=str(exc)[:220])
    write_json(source_dir / "receipt.json", receipt)
    return receipt


def build_indexes(data, check=False):
    methods = []
    source_records = []
    for source in data["sources"]:
        receipt = CACHE / source["id"] / "receipt.json"
        acquired_record = read_json(receipt) if receipt.exists() else source.get("acquisition_snapshot", {"status": "not_acquired"})
        # Public indexes describe portable, dated receipts rather than machine-specific absolute paths.
        acquired_record = dict(acquired_record)
        acquired_record.pop("previous_archive_attempt", None)
        if acquired_record.get("source_tree"):
            tree = Path(acquired_record["source_tree"])
            if tree.is_absolute() and tree.is_relative_to(ROOT):
                acquired_record["source_tree"] = tree.relative_to(ROOT).as_posix()
        source_records.append({**source, "acquisition": acquired_record})
    acquired = {s["id"] for s in source_records if s["acquisition"]["status"] in {"source_acquired", "source_subset_acquired"} and s["acquisition"].get("code_files", 0) > 0}
    for method in data["methods"]:
        methods.append({**method, "source_repository_code_acquired": bool(set(method.get("source_ids", [])) & acquired),
                        "registry_path": REGISTRY.relative_to(ROOT).as_posix(), "source_collection": "methods"})
    report = {"records_including_overviews": len(methods),
              "detailed_records_including_designs_and_platforms": sum(r.get("record_kind") != "family_overview" for r in methods),
              "family_overviews": sum(r.get("record_kind") == "family_overview" for r in methods),
              "record_kinds": {k:sum(r.get("record_kind") == k for r in methods) for k in sorted({r.get("record_kind", "unspecified") for r in methods})},
              "by_category": {c:sum(r["category"] == c for r in methods) for c in sorted(CATEGORIES)},
              "detailed_by_category": {c:sum(r["category"] == c and r.get("record_kind") != "family_overview" for r in methods) for c in sorted(CATEGORIES)},
              "official_sources":len(source_records), "sources_with_code_acquired":len(acquired),
              "source_archives_acquired":sum(s["acquisition"]["status"] == "source_acquired" for s in source_records),
              "source_subsets_acquired":sum(s["acquisition"]["status"] == "source_subset_acquired" for s in source_records),
              "detailed_records_with_source_repository_code":sum(r.get("record_kind") != "family_overview" and r["source_repository_code_acquired"] for r in methods),
              "count_policy":"Records include variants, designs, diagnostics and software platforms; this is not a count of distinct MR estimators.",
              "runtime_or_scientific_validation_not_inferred_from_source_count":True}
    products = {"mr_coloc_methods.jsonl": "".join(json.dumps(r, ensure_ascii=False) + "\n" for r in methods),
                "mr_coloc_sources.jsonl": "".join(json.dumps(r, ensure_ascii=False) + "\n" for r in source_records),
                "mr_coloc_summary.json":json.dumps(report, ensure_ascii=False, indent=2) + "\n"}
    stale = []
    for name, text in products.items():
        path = ROOT / "catalog" / name
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            stale.append(name)
            if not check:
                path.write_text(text, encoding="utf-8", newline="\n")
    if check and stale:
        raise ValueError("Stale method indexes: " + ", ".join(stale))
    return report


def main(argv=None):
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    q = sub.add_parser("query"); q.add_argument("terms", nargs="*", default=[]); q.add_argument("--category", choices=sorted(CATEGORIES)); q.add_argument("--limit", type=int, default=20)
    s = sub.add_parser("select"); s.add_argument("context", type=Path)
    f = sub.add_parser("fetch"); f.add_argument("--source", action="append"); f.add_argument("--all", action="store_true"); f.add_argument("--max-total-mb", type=int, default=250)
    fc = sub.add_parser("fetch-code"); fc.add_argument("--source", action="append"); fc.add_argument("--all", action="store_true")
    b = sub.add_parser("build"); b.add_argument("--check", action="store_true")
    sub.add_parser("verify-sources")
    run = sub.add_parser("run"); run.add_argument("--config", type=Path, required=True); run.add_argument("--rscript", default="Rscript")
    sub.add_parser("check")
    args = parser.parse_args(argv)
    data = load_registry()
    if args.command == "query":
        if args.limit < 1:
            raise ValueError("query --limit must be positive")
        rows = query(data, " ".join(args.terms), args.category)
        print(json.dumps({"total_matches": len(rows), "returned": min(len(rows), args.limit), "methods": rows[:args.limit]}, ensure_ascii=False, indent=2))
    elif args.command == "select":
        print(json.dumps(select_methods(data, read_json(args.context)), ensure_ascii=False, indent=2))
    elif args.command == "verify-sources":
        reports = [verify_source(s) for s in data["sources"]]
        print(json.dumps(reports, ensure_ascii=False, indent=2))
        return int(any(r["status"] == "integrity_failed" for r in reports))
    elif args.command == "run":
        config = read_json(args.config)
        engine = config.get("engine")
        if engine not in {"TwoSampleMR", "MendelianRandomization", "coloc.abf"}:
            raise ValueError("Supported native engines: TwoSampleMR, MendelianRandomization, coloc.abf")
        script = HERE / ("run_coloc_abf.R" if engine == "coloc.abf" else "run_native_mr.R")
        return subprocess.run([args.rscript, str(script), "--config", str(args.config.resolve())]).returncode
    elif args.command in {"fetch", "fetch-code"}:
        if not args.all and not args.source:
            parser.error("fetch needs explicit --source IDs or --all")
        sources = data["sources"] if args.all else [r for r in data["sources"] if r["id"] in args.source]
        if args.source and set(args.source) - {r["id"] for r in sources}:
            parser.error("Unknown source ID")
        consumed = 0
        for source in sources:
            if args.command == "fetch-code":
                result = fetch_code_subset(source)
                print(json.dumps({k: result.get(k) for k in ("source_id", "status", "commit", "subset_bytes", "code_files", "error")}), flush=True)
                continue
            if consumed >= args.max_total_mb * 1_000_000:
                raise ValueError("Batch acquisition bound reached")
            result = fetch_source(source, max_archive_bytes=min(30_000_000, args.max_total_mb * 1_000_000 - consumed))
            consumed += result.get("network_bytes_received", result.get("archive_bytes", 0))
            print(json.dumps({k: result.get(k) for k in ("source_id", "status", "commit", "archive_bytes", "code_files", "error")}), flush=True)
    else:
        report = build_indexes(data, check=args.command == "check" or (args.command == "build" and args.check))
        if args.command == "check":
            report["registry_valid"] = True
        print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
