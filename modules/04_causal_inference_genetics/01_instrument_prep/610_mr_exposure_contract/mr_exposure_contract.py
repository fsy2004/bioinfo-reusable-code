"""Audit MR exposure meanings and study metadata; no causal estimate is fitted."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

PROXY_KINDS = {"prescribed_use_propensity", "endogenous_cortisol", "gene_expression",
               "protein_abundance", "misuse_or_abuse"}
REQUIRED = ("id", "phenotype", "proxy_kind", "source", "genome_build", "ancestry", "beta_unit", "n")


def audit_contract(contract, instruments):
    issues = []
    def add(level, code, message, exposure_id=None):
        issues.append({"level": level, "code": code, "exposure_id": exposure_id, "message": message})

    exposures = contract.get("exposures", [])
    if not exposures:
        add("ERROR", "no_exposures", "Declare at least one exposure")
    if not contract.get("research_question") or not contract.get("intended_estimand"):
        add("ERROR", "estimand", "Declare research_question and intended_estimand")
    if contract.get("design") not in ("UVMR", "MVMR", "cis_MR"):
        add("ERROR", "design", "design must be UVMR, MVMR or cis_MR")
    outcome = contract.get("outcome", {})
    for key in ("id", "phenotype", "source", "genome_build", "ancestry", "effect_scale", "n"):
        if outcome.get(key) in (None, ""):
            add("ERROR", "outcome_metadata", key)
    try:
        if not math.isfinite(float(outcome["n"])) or float(outcome["n"]) <= 0:
            raise ValueError
    except (KeyError, ValueError, TypeError):
        add("ERROR", "outcome_n", "n must be finite and positive")
    if outcome.get("effect_scale") == "log_odds":
        for key in ("n_cases", "n_controls"):
            try:
                if not math.isfinite(float(outcome[key])) or float(outcome[key]) <= 0:
                    raise ValueError
            except (KeyError, ValueError, TypeError):
                add("ERROR", "outcome_case_control", key + " is required and positive")
    overlap = contract.get("sample_overlap", {})
    if overlap.get("status") not in ("none_documented", "overlap_documented", "unknown"):
        add("ERROR", "overlap_status", "Explicitly declare none_documented, overlap_documented or unknown")
    elif overlap["status"] == "unknown":
        add("WARNING", "overlap_unknown", "Sample overlap remains unknown; absence of documentation is not independence")
    elif not overlap.get("evidence"):
        add("ERROR", "overlap_evidence", "Documented overlap/independence needs source evidence")
    if "fraction" in overlap and (not isinstance(overlap["fraction"], (int, float))
                                  or not 0 <= overlap["fraction"] <= 1):
        add("ERROR", "overlap_fraction", "fraction must be within [0,1], or omit when unknown")

    known_ids = set()
    for exposure in exposures:
        eid = exposure.get("id")
        if eid in known_ids:
            add("ERROR", "duplicate_exposure", str(eid), eid)
        known_ids.add(eid)
        for key in REQUIRED:
            if exposure.get(key) in (None, ""):
                add("ERROR", "exposure_metadata", key, eid)
        kind = exposure.get("proxy_kind")
        if kind not in PROXY_KINDS:
            add("ERROR", "proxy_kind", "Unknown proxy_kind", eid)
        if kind == "endogenous_cortisol":
            add("NOTICE", "cortisol_semantics", "Endogenous cortisol variation does not measure exogenous dose, duration or misuse", eid)
        if kind == "prescribed_use_propensity":
            add("NOTICE", "medication_semantics", "Use propensity includes indication, access and prescribing pathways; it is not a pharmacological dose effect", eid)
        if kind in ("gene_expression", "protein_abundance"):
            for key in ("target_gene", "tissue_or_cell_type", "cis_window", "molecular_qtl_type"):
                if not exposure.get(key):
                    add("ERROR", "molecular_proxy", key + " required", eid)
            add("NOTICE", "target_semantics", "Molecular target perturbation requires direction/on-target validation; expression is not receptor activity or drug misuse", eid)
        if kind == "misuse_or_abuse":
            if not exposure.get("measured_definition") or exposure.get("validated_misuse_phenotype") is not True:
                add("ERROR", "unvalidated_abuse", "Misuse/abuse requires an explicitly measured, validated phenotype and its definition", eid)
        if contract.get("claim") == "glucocorticoid_abuse_causal_effect" and kind != "misuse_or_abuse":
            add("ERROR", "construct_mismatch", "This proxy cannot support the requested glucocorticoid-abuse estimand; state a separate biological question", eid)
        if exposure.get("genome_build") != outcome.get("genome_build"):
            add("ERROR", "build_mismatch", "Explicit harmonization/liftover is required before combining builds", eid)
        if exposure.get("ancestry") != outcome.get("ancestry"):
            add("WARNING", "ancestry_mismatch", "Transportability and ancestry-matched LD require an explicit analysis plan", eid)
        try:
            if not math.isfinite(float(exposure["n"])) or float(exposure["n"]) <= 0:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            add("ERROR", "exposure_n", "n must be finite and positive", eid)

    strengths = {}
    seen = set()
    for row_num, row in enumerate(instruments, 2):
        eid = row.get("exposure_id")
        key = (eid, row.get("variant_id"))
        if eid not in known_ids:
            add("ERROR", "instrument_exposure", "Unknown exposure_id at row " + str(row_num), eid)
            continue
        if not row.get("variant_id") or key in seen:
            add("ERROR", "duplicate_instrument", "Empty/duplicate variant_id at row " + str(row_num), eid)
        seen.add(key)
        try:
            beta, se = float(row["beta"]), float(row["se"])
            if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
                raise ValueError
            f_value = (beta / se) ** 2
            strengths.setdefault(eid, []).append({"variant_id": row["variant_id"], "marginal_F_approx": f_value})
        except (KeyError, ValueError, OverflowError):
            add("ERROR", "instrument_numeric", "Finite beta and positive se required at row " + str(row_num), eid)
    for eid in known_ids:
        vals = strengths.get(eid, [])
        if not vals:
            add("ERROR", "no_instruments", "No valid instruments supplied", eid)
        elif min(item["marginal_F_approx"] for item in vals) < 10:
            add("WARNING", "weak_marginal_iv", "At least one approximate marginal F < 10; this threshold does not establish validity", eid)
    conditional = {"status": "not_applicable"}
    if contract.get("design") == "MVMR":
        evidence = contract.get("conditional_strength", {})
        provided = evidence.get("F_by_exposure", {})
        values_valid = all(eid in provided and isinstance(provided[eid], (int, float))
                           and math.isfinite(provided[eid]) and provided[eid] >= 0 for eid in known_ids)
        if not values_valid or not evidence.get("method") or not evidence.get("covariance_source"):
            conditional = {"status": "UNAVAILABLE", "F_by_exposure": None}
            add("ERROR", "conditional_strength_unavailable", "Supply computed conditional F and exposure-association covariance provenance; marginal F cannot substitute")
        else:
            conditional = {"status": "USER_SUPPLIED_UNVERIFIED", **evidence}
            if any(provided[eid] < 10 for eid in known_ids):
                add("WARNING", "weak_conditional_iv", "At least one supplied conditional F < 10; consider weak-IV robust MVMR")
    errors = sum(item["level"] == "ERROR" for item in issues)
    return {"module": "610_mr_exposure_contract", "status": "FAIL" if errors else "REVIEW_REQUIRED" if any(i["level"] == "WARNING" for i in issues) else "PASS_METADATA_ONLY",
            "errors": errors, "issues": issues, "marginal_instrument_strength": strengths,
            "conditional_strength": conditional,
            "scope": "Metadata and semantic audit only; no causal effect, conditional F or data validity is inferred"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--instruments", required=True, help="TSV with exposure_id, variant_id, beta, se")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    output = Path(args.out)
    if output.exists():
        parser.error("Output already exists; choose a new path to preserve existing results")
    try:
        with Path(args.contract).open(encoding="utf-8-sig") as handle:
            contract = json.load(handle)
        with Path(args.instruments).open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)) or not {"exposure_id", "variant_id", "beta", "se"}.issubset(fields):
                raise ValueError("Instrument table requires unique exposure_id, variant_id, beta, se columns")
            instruments = list(reader)
            if any(None in row or any(value is None for value in row.values()) for row in instruments):
                raise ValueError("Ragged instrument table row")
        report = audit_contract(contract, instruments)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(report["status"])
    return 2 if report["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
