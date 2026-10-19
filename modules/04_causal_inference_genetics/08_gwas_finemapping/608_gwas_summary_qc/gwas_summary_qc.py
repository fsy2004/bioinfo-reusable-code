"""Audit dense regional GWAS summary statistics and a signed LD matrix.

This validates supplied inputs. It neither performs GWAS nor clumps variants,
harmonizes alleles, infers ancestry, or changes genome builds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

REQUIRED = ("variant_id", "chr", "pos", "ea", "nea", "beta", "se", "p", "eaf")
META_REQUIRED = ("trait_id", "phenotype", "genome_build", "ancestry", "beta_unit",
                 "effect_type", "source", "region", "dense_region_confirmed")


def read_table(path):
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if not sample.strip():
            raise ValueError("Empty table: " + str(path))
        delimiter = "\t" if "\t" in sample.splitlines()[0] else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames or []
        if not fields or len(set(fields)) != len(fields) or any(not name.strip() for name in fields):
            raise ValueError("Empty or duplicate table header: " + str(path))
        rows = []
        for line_no, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Ragged table row " + str(line_no) + ": " + str(path))
            rows.append(row)
        return rows, fields


def file_record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.resolve()), "sha256": digest.hexdigest()}


def audit(summary_path, metadata_path, ld_path=None, ld_variants_path=None):
    rows, fields = read_table(summary_path)
    with Path(metadata_path).open(encoding="utf-8-sig") as handle:
        meta = json.load(handle)
    findings = []

    def add(level, code, message, row=None):
        findings.append({"level": level, "code": code, "row": row, "message": message})

    for field in REQUIRED:
        if field not in fields:
            add("ERROR", "missing_column", field)
    for field in META_REQUIRED:
        if field not in meta or meta[field] in (None, "", []):
            add("ERROR", "missing_metadata", field)
    if meta.get("genome_build") not in ("GRCh37", "GRCh38"):
        add("ERROR", "genome_build", "genome_build must explicitly be GRCh37 or GRCh38")
    if meta.get("effect_type") not in ("linear_beta", "log_odds"):
        add("ERROR", "effect_type", "Only linear_beta or log_odds is supported; convert OR before use")
    if meta.get("dense_region_confirmed") is not True:
        add("ERROR", "sparse_or_unconfirmed", "Dense unfiltered regional input must be explicitly confirmed")
    def chromosome(value):
        value = re.sub(r"^chr", "", str(value).strip(), flags=re.IGNORECASE).upper()
        if value == "M":
            value = "MT"
        if value in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
            return value
        return None

    region_match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", str(meta.get("region", "")).strip())
    region = None
    if region_match:
        region_chr = chromosome(region_match.group(1))
        region_start, region_end = int(region_match.group(2)), int(region_match.group(3))
        if region_chr and 0 < region_start <= region_end:
            region = (region_chr, region_start, region_end)
    if region is None:
        add("ERROR", "region_format", "region must declare a valid chromosome and 1-based inclusive interval, e.g. chr1:100-200")
    if not rows:
        add("ERROR", "empty_summary", "No variants supplied")
    ids = set()
    positions = set()
    parsed = []
    for index, row in enumerate(rows, 2):
        variant = row.get("variant_id", "").strip()
        if not variant or variant in ids:
            add("ERROR", "variant_id", "Empty or duplicate variant_id: " + variant, index)
        ids.add(variant)
        ea, nea = row.get("ea", "").upper(), row.get("nea", "").upper()
        if ea == nea or not ea or not nea or any(c not in "ACGT" for c in ea + nea):
            add("ERROR", "allele", "Alleles must be distinct explicit DNA sequences; normalize indels upstream", index)
        if {ea, nea} in ({"A", "T"}, {"C", "G"}):
            add("WARNING", "palindromic", "Palindromic allele pair requires reference/frequency-supported orientation", index)
        try:
            pos = int(row.get("pos", ""))
            if pos < 1:
                raise ValueError
            row_chr = chromosome(row.get("chr", ""))
            if row_chr is None:
                add("ERROR", "chromosome", "chr must be 1-22, X, Y or MT (optional chr prefix)", index)
            elif region and (row_chr != region[0] or not region[1] <= pos <= region[2]):
                add("ERROR", "outside_declared_region", "Variant chromosome/position lies outside the declared inclusive region", index)
            key = (row_chr, pos, ea, nea)
            if key in positions:
                add("ERROR", "duplicate_allele_coordinate", str(key), index)
            positions.add(key)
        except ValueError:
            add("ERROR", "position", "pos must be a positive integer", index)
        numbers = {}
        for field in ("beta", "se", "p", "eaf"):
            try:
                value = float(row.get(field, ""))
                if not math.isfinite(value):
                    raise ValueError
                numbers[field] = value
            except ValueError:
                add("ERROR", "numeric", field + " must be finite", index)
        if numbers.get("se", 1) <= 0:
            add("ERROR", "se", "se must be greater than zero", index)
        if "p" in numbers and not 0 <= numbers["p"] <= 1:
            add("ERROR", "p_range", "p must be within [0,1]", index)
        if "eaf" in numbers and not 0 < numbers["eaf"] < 1:
            add("ERROR", "eaf", "Effect allele frequency must be strictly within (0,1)", index)
        if "n" in row:
            try:
                if not math.isfinite(float(row["n"])) or float(row["n"]) <= 0:
                    raise ValueError
            except ValueError:
                add("ERROR", "sample_size", "n must be finite and positive", index)
        if all(k in numbers for k in ("beta", "se", "p")) and numbers["se"] > 0:
            p_wald = math.erfc(abs(numbers["beta"] / numbers["se"]) / math.sqrt(2))
            p_given = numbers["p"]
            if p_given == 0:
                add("WARNING", "p_underflow", "p=0 requires original log-p or sufficiently precise statistics", index)
            elif p_wald > 0 and abs(math.log10(p_wald) - math.log10(p_given)) > 0.5:
                add("WARNING", "p_beta_se", "P differs from two-sided normal Wald approximation; inspect test type/rounding", index)
        parsed.append({"variant_id": variant, "ea": ea, "nea": nea})

    ld_checks = {"provided": ld_path is not None, "psd_test": "not_run"}
    if (ld_path is None) != (ld_variants_path is None):
        add("ERROR", "ld_pair", "Provide both --ld-matrix and --ld-variants")
    elif ld_path is not None:
        for field in ("ld_genome_build", "ld_ancestry", "ld_source", "ld_coded_allele", "ld_metric"):
            if not meta.get(field):
                add("ERROR", "ld_metadata", field)
        if meta.get("ld_genome_build") != meta.get("genome_build"):
            add("ERROR", "ld_build_mismatch", "LD and summary builds differ")
        if meta.get("ld_ancestry") != meta.get("ancestry"):
            add("ERROR", "ld_ancestry_mismatch", "LD ancestry differs; use matched reference before this entry")
        if meta.get("ld_metric") != "signed_r" or meta.get("ld_coded_allele") != "ea":
            add("ERROR", "ld_orientation", "LD must be signed_r coded to ea; r-squared is insufficient")
        ld_rows, ld_fields = read_table(ld_variants_path)
        if not all(k in ld_fields for k in ("variant_id", "ea", "nea")):
            add("ERROR", "ld_variant_columns", "LD variant table needs variant_id, ea, nea")
        if len(ld_rows) != len(parsed):
            add("ERROR", "ld_variant_count", "LD variant count differs from summary")
        else:
            for index, (item, ld_item) in enumerate(zip(parsed, ld_rows), 1):
                if any(item[key] != (ld_item.get(key, "").upper() if key != "variant_id"
                                    else ld_item.get(key)) for key in ("variant_id", "ea", "nea")):
                    add("ERROR", "ld_order_alleles", "SNP order/allele mismatch at LD position " + str(index))
        try:
            with Path(ld_path).open(encoding="utf-8-sig") as handle:
                matrix = [[float(value) for value in line.split()] for line in handle if line.strip()]
            count = len(parsed)
            if len(matrix) != count or any(len(row) != count for row in matrix):
                add("ERROR", "ld_shape", "LD must be a square matrix with one row per summary variant")
            elif any(not math.isfinite(value) for row in matrix for value in row):
                add("ERROR", "ld_finite", "LD contains non-finite values")
            else:
                if any(abs(matrix[i][i] - 1) > 1e-6 for i in range(count)):
                    add("ERROR", "ld_diagonal", "LD diagonal must equal 1")
                if any(abs(value) > 1 + 1e-6 for row in matrix for value in row):
                    add("ERROR", "ld_range", "Signed correlation must be within [-1,1]")
                if any(abs(matrix[i][j] - matrix[j][i]) > 1e-6 for i in range(count) for j in range(i)):
                    add("ERROR", "ld_symmetry", "LD is not symmetric")
                try:
                    import numpy as np
                    eigen_min = float(np.linalg.eigvalsh(np.asarray(matrix)).min())
                    ld_checks.update(psd_test="numpy_eigvalsh", minimum_eigenvalue=eigen_min)
                    if eigen_min < -1e-6:
                        add("ERROR", "ld_psd", "LD is not positive semidefinite")
                except ImportError:
                    add("WARNING", "psd_unchecked", "numpy unavailable; symmetry/range checks do not establish PSD")
        except ValueError:
            add("ERROR", "ld_parse", "LD must be headerless whitespace-separated numbers")
    inputs = [file_record(summary_path), file_record(metadata_path)]
    inputs.extend(file_record(path) for path in (ld_path, ld_variants_path) if path)
    errors = sum(item["level"] == "ERROR" for item in findings)
    warnings = sum(item["level"] == "WARNING" for item in findings)
    return {"module": "608_gwas_summary_qc", "schema_version": "1.0",
            "status": "FAIL" if errors else "PASS_WITH_WARNINGS" if warnings else "PASS",
            "n_variants": len(rows), "errors": errors, "warnings": warnings,
            "metadata": meta, "ld_checks": ld_checks, "findings": findings, "inputs": inputs,
            "scope": "Input audit only; declarations are not independently verified population facts"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--ld-matrix")
    parser.add_argument("--ld-variants")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out = Path(args.out_dir)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error("Output directory must be absent or empty; existing results are preserved")
    try:
        report = audit(args.summary, args.metadata, args.ld_matrix, args.ld_variants)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    out.mkdir(parents=True, exist_ok=True)
    with (out / "audit.json").open("x", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(json.dumps({key: report[key] for key in ("status", "n_variants", "errors", "warnings")}))
    return 2 if report["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
