"""Audit medication records using explicit dictionaries; no abuse classification."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def normalized(value: str) -> str:
    return " ".join((value or "").strip().casefold().split())


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate CSV headers: {path}")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"Ragged CSV rows: {path}")
        return reader.fieldnames, rows


def dictionary(path: Path, required: set[str]) -> dict[str, dict[str, str]]:
    headers, rows = read_csv(path)
    if not required.issubset(headers):
        raise ValueError(f"{path} requires columns {sorted(required)}")
    result = {}
    for row in rows:
        key = normalized(row["alias"])
        if not key:
            raise ValueError(f"Empty dictionary alias: {path}")
        if key in result and result[key] != row:
            raise ValueError(f"Conflicting normalized alias {key!r} in {path}")
        result[key] = row
    return result


def duration(raw: str) -> tuple[float | None, str]:
    if not raw.strip():
        return None, "duration_missing"
    try:
        value = float(raw)
    except ValueError:
        return None, "duration_non_numeric"
    if not math.isfinite(value) or value < 0:
        return None, "duration_invalid"
    return value, ""


def audit(input_path: Path, drug_map_path: Path, route_map_path: Path, output_dir: Path) -> dict:
    headers, rows = read_csv(input_path)
    required = {"person_id", "drug_name", "route", "duration_days"}
    if not required.issubset(headers):
        raise ValueError(f"Input requires columns {sorted(required)}")
    if any(name.startswith("audit_") for name in headers):
        raise ValueError("Input columns cannot use the reserved audit_ prefix")
    if any(not row["person_id"].strip() for row in rows):
        raise ValueError("Every medication record requires a nonempty person_id")
    drug_map = dictionary(drug_map_path, {"alias", "ingredient", "drug_class"})
    route_map = dictionary(route_map_path, {"alias", "route_type"})
    allowed_routes = {"systemic", "inhaled", "topical", "local_other", "unknown"}
    if any(row["route_type"] not in allowed_routes for row in route_map.values()):
        raise ValueError(f"route_type must be one of {sorted(allowed_routes)}")
    if any(not row["ingredient"].strip() or not row["drug_class"].strip() for row in drug_map.values()):
        raise ValueError("Every drug mapping needs ingredient and drug_class")
    annotated, people, flags_count = [], defaultdict(list), Counter()
    seen = Counter()
    for record_index, raw in enumerate(rows, 1):
        flags = []
        drug_key, route_key = normalized(raw["drug_name"]), normalized(raw["route"])
        drug = drug_map.get(drug_key)
        route = route_map.get(route_key)
        if drug is None:
            flags.append("drug_missing" if not drug_key else "drug_unmapped")
        if route is None:
            flags.append("route_missing" if not route_key else "route_unmapped")
        elif route["route_type"] == "unknown":
            flags.append("route_unknown")
        days, issue = duration(raw["duration_days"])
        if issue:
            flags.append(issue)
        # This flags repeat source records; it never deletes them.
        key = tuple(raw[name] for name in headers)
        seen[key] += 1
        if seen[key] > 1:
            flags.append("repeat_identical_record")
        is_gcs = drug is not None and normalized(drug["drug_class"]) == "glucocorticoid"
        result = dict(raw)
        result.update(audit_record_index=record_index, audit_name_normalized=drug_key,
                      audit_ingredient=drug["ingredient"] if drug else "",
                      audit_drug_class=drug["drug_class"] if drug else "",
                      audit_route_type=route["route_type"] if route else "unknown",
                      audit_duration_days=days if days is not None else "",
                      audit_gcs_record=int(is_gcs), audit_flags=";".join(flags))
        annotated.append(result)
        people[raw["person_id"]].append(result)
        flags_count.update(flags)
    person_rows = []
    for person_id, records in people.items():
        gcs = [record for record in records if record["audit_gcs_record"]]
        valid_days = [r["audit_duration_days"] for r in gcs if r["audit_duration_days"] != ""]
        person_rows.append({"person_id": person_id, "observed_records": len(records),
                            "mapped_gcs_records": len(gcs),
                            "any_observed_systemic_gcs": int(any(r["audit_route_type"] == "systemic" for r in gcs)),
                            "any_unmapped_drug": int(any(not r["audit_ingredient"] for r in records)),
                            "any_gcs_unknown_route": int(any(r["audit_route_type"] == "unknown" for r in gcs)),
                            "max_reported_gcs_record_duration_days": max(valid_days) if valid_days else ""})
    destinations = [output_dir / name for name in ("medication_audit.csv", "person_record_summary.csv", "audit.json")]
    if any(path.exists() for path in destinations):
        raise FileExistsError("Output files already exist; select a new output directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_fields = headers + ["audit_record_index", "audit_name_normalized", "audit_ingredient", "audit_drug_class",
                              "audit_route_type", "audit_duration_days", "audit_gcs_record", "audit_flags"]
    person_fields = ["person_id", "observed_records", "mapped_gcs_records", "any_observed_systemic_gcs",
                     "any_unmapped_drug", "any_gcs_unknown_route", "max_reported_gcs_record_duration_days"]
    for path, fields, records in ((destinations[0], output_fields, annotated), (destinations[1], person_fields, person_rows)):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
    report = {"records": len(rows), "observed_people": len(people), "flags": dict(flags_count),
              "input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in (input_path, drug_map_path, route_map_path)},
              "interpretation": "Observed prescription records only. No abuse, cumulative dose, continuous exposure, or unexposed status inferred.",
              "mapping": "Exact normalized alias; raw input columns retained. Duration maxima are not summed or treated as continuous exposure."}
    destinations[2].write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--drug-map", type=Path, required=True)
    parser.add_argument("--route-map", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.input, args.drug_map, args.route_map, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
