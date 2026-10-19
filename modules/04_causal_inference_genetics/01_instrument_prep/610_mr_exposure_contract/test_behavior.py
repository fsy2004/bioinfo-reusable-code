"""Synthetic semantics and instrument-strength tests; no real GWAS data."""
import unittest
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from mr_exposure_contract import audit_contract


def contract():
    return dict(research_question="synthetic cortisol biology", intended_estimand="endogenous cortisol", design="UVMR",
                exposures=[dict(id="x", phenotype="morning cortisol", proxy_kind="endogenous_cortisol", source="synthetic", genome_build="GRCh38", ancestry="EUR", beta_unit="SD", n=1000)],
                outcome=dict(id="y", phenotype="cancer", source="synthetic", genome_build="GRCh38", ancestry="EUR", effect_scale="log_odds", n=2000, n_cases=200, n_controls=1800),
                sample_overlap=dict(status="unknown"))


class ContractTests(unittest.TestCase):
    def test_approximate_marginal_f(self):
        result = audit_contract(contract(), [dict(exposure_id="x", variant_id="rs1", beta="0.5", se="0.1")])
        self.assertEqual(result["errors"], 0)
        self.assertAlmostEqual(result["marginal_instrument_strength"]["x"][0]["marginal_F_approx"], 25)

    def test_cortisol_cannot_be_abuse_claim(self):
        data = contract(); data["claim"] = "glucocorticoid_abuse_causal_effect"
        result = audit_contract(data, [dict(exposure_id="x", variant_id="rs1", beta="0.5", se="0.1")])
        self.assertIn("construct_mismatch", {item["code"] for item in result["issues"]})
        self.assertEqual(result["status"], "FAIL")

    def test_mvmr_does_not_substitute_marginal_f(self):
        data = contract(); data["design"] = "MVMR"
        result = audit_contract(data, [dict(exposure_id="x", variant_id="rs1", beta="0.5", se="0.1")])
        self.assertIsNone(result["conditional_strength"]["F_by_exposure"])
        self.assertEqual(result["conditional_strength"]["status"], "UNAVAILABLE")

    def test_nonfinite_outcome_sample_sizes_are_rejected(self):
        data = contract(); data["outcome"].update(n=float("nan"), n_cases=float("inf"), n_controls=float("nan"))
        result = audit_contract(data, [dict(exposure_id="x", variant_id="rs1", beta="0.5", se="0.1")])
        codes = {item["code"] for item in result["issues"]}
        self.assertIn("outcome_n", codes)
        self.assertIn("outcome_case_control", codes)

    def test_existing_output_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            contract_path = root / "contract.json"; contract_path.write_text(json.dumps(contract()), encoding="utf-8")
            instruments = root / "instruments.tsv"; instruments.write_text("exposure_id\tvariant_id\tbeta\tse\nx\trs1\t0.5\t0.1\n", encoding="utf-8")
            output = root / "audit.json"; output.write_text("preserve", encoding="utf-8")
            result = subprocess.run([sys.executable, str(Path(__file__).with_name("mr_exposure_contract.py")), "--contract", str(contract_path), "--instruments", str(instruments), "--out", str(output)], capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertEqual(output.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main()
