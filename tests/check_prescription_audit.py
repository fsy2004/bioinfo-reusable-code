"""Deterministic artificial records exercise the real audit API."""
import csv
import importlib.util
import json
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[1] / "modules" / "17_population_epidemiology" / "05_medication_exposure"
spec = importlib.util.spec_from_file_location("exposure", root / "604_prescription_exposure_audit" / "604_prescription_exposure_audit.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
base = Path(tempfile.mkdtemp(prefix="prescription_interface_")) / "artificial"
base.mkdir(parents=True, exist_ok=False)


def write(name, fields, rows):
    path = base / name
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


drug = write("artificial_drugs.csv", ["alias", "ingredient", "drug_class"], [
    {"alias": "TEST_GCS", "ingredient": "ARTIFICIAL_INGREDIENT", "drug_class": "glucocorticoid"}])
route = write("artificial_routes.csv", ["alias", "route_type"], [{"alias": "TEST_ORAL", "route_type": "systemic"}])
raw = {"person_id": "001", "drug_name": "  TEST_GCS  ", "route": "TEST_ORAL", "duration_days": "10", "source_note": "retain exactly"}
rows = [raw, dict(raw), {"person_id": "002", "drug_name": "UNMAPPED", "route": "", "duration_days": "-2", "source_note": "unknown"},
        {"person_id": "003", "drug_name": "TEST_GCS", "route": "OTHER", "duration_days": "NaN", "source_note": "invalid"}]
inp = write("artificial_input.csv", list(raw), rows)
receipt = module.audit(inp, drug, route, base / "results")
headers, annotated = module.read_csv(base / "results" / "medication_audit.csv")
assert receipt["records"] == 4 and receipt["observed_people"] == 3
assert annotated[0]["drug_name"] == "  TEST_GCS  " and annotated[0]["source_note"] == "retain exactly"
assert "repeat_identical_record" in annotated[1]["audit_flags"]
assert annotated[2]["audit_gcs_record"] == "0" and "drug_unmapped" in annotated[2]["audit_flags"]
assert "route_missing" in annotated[2]["audit_flags"] and "duration_invalid" in annotated[2]["audit_flags"]
assert "duration_invalid" in annotated[3]["audit_flags"] and "route_unmapped" in annotated[3]["audit_flags"]
try:
    module.audit(inp, drug, route, base / "results")
except FileExistsError:
    pass
else:
    raise AssertionError("Existing results must not be overwritten")
conflict = write("conflicting_dictionary.csv", ["alias", "ingredient", "drug_class"], [
    {"alias": "TEST_GCS", "ingredient": "A", "drug_class": "glucocorticoid"},
    {"alias": " test_gcs ", "ingredient": "B", "drug_class": "glucocorticoid"}])
try:
    module.audit(inp, conflict, route, base / "must_not_exist")
except ValueError:
    pass
else:
    raise AssertionError("Conflicting mappings must fail")
assert not (base / "must_not_exist").exists()
report = {"fixture": "artificial records only", "checks": ["raw preservation", "exact normalized mapping", "unmapped names/routes", "invalid durations", "duplicate retention", "overwrite rejection", "conflicting dictionary rejection"], "status": "passed"}
(base / "validation_receipt.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
