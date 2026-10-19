"""Synthetic input-contract tests, not a scientific validation dataset."""
import csv
import json
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path
from gwas_summary_qc import audit, read_table


class AuditTests(unittest.TestCase):
    def fixture(self, root):
        summary = root / "summary.tsv"
        summary.write_text("variant_id\tchr\tpos\tea\tnea\tbeta\tse\tp\teaf\nrs1\t1\t100\tA\tC\t0\t0.1\t1\t0.2\nrs2\t1\t101\tG\tT\t0\t0.1\t1\t0.3\n", encoding="utf-8")
        metadata = root / "metadata.json"
        metadata.write_text(json.dumps(dict(trait_id="synthetic", phenotype="synthetic", genome_build="GRCh38", ancestry="EUR", beta_unit="SD", effect_type="linear_beta", source="synthetic test only", region="1:100-101", dense_region_confirmed=True, ld_genome_build="GRCh38", ld_ancestry="EUR", ld_source="synthetic", ld_metric="signed_r", ld_coded_allele="ea")), encoding="utf-8")
        ldv = root / "ld_variants.tsv"
        ldv.write_text("variant_id\tea\tnea\nrs1\tA\tC\nrs2\tG\tT\n", encoding="utf-8")
        ld = root / "ld.txt"
        ld.write_text("1 0.2\n0.2 1\n", encoding="utf-8")
        return summary, metadata, ld, ldv

    def test_matching_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = audit(*self.fixture(Path(tmp)))
            self.assertEqual(result["errors"], 0)
            self.assertEqual(result["n_variants"], 2)

    def test_declared_region_and_chromosome_are_enforced(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self.fixture(Path(tmp))
            metadata = json.loads(paths[1].read_text())
            metadata['region'] = 'chr1:100-101'
            paths[1].write_text(json.dumps(metadata))
            self.assertEqual(audit(*paths)['errors'], 0)
            original = paths[0].read_text()
            paths[0].write_text(original.replace('rs1\t1\t100', 'rs1\t2\t1000'))
            self.assertIn('outside_declared_region', {r['code'] for r in audit(*paths)['findings']})
            paths[0].write_text(original.replace('rs1\t1\t100', 'rs1\tUNKNOWN\t100'))
            self.assertIn('chromosome', {r['code'] for r in audit(*paths)['findings']})
            paths[0].write_text(original)
            metadata['region'] = 'chr1:200-100'
            paths[1].write_text(json.dumps(metadata))
            self.assertIn('region_format', {r['code'] for r in audit(*paths)['findings']})

    def test_ld_allele_order_and_invalid_se_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self.fixture(Path(tmp))
            paths[3].write_text("variant_id\tea\tnea\nrs2\tG\tT\nrs1\tA\tC\n", encoding="utf-8")
            paths[0].write_text(paths[0].read_text().replace("0\t0.1", "0\t0"), encoding="utf-8")
            result = audit(*paths)
            codes = {item["code"] for item in result["findings"]}
            self.assertIn("ld_order_alleles", codes)
            self.assertIn("se", codes)
            self.assertEqual(result["status"], "FAIL")

    def test_empty_duplicate_header_and_ragged_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            file = Path(tmp) / "bad.tsv"
            for content in ("", "x\tx\n1\t2\n", "x\ty\n1\n", "x\ty\n1\t2\t3\n"):
                file.write_text(content, encoding="utf-8")
                with self.assertRaises(ValueError):
                    read_table(file)

    def test_existing_output_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self.fixture(root)
            output = root / "output"; output.mkdir()
            report = output / "audit.json"; report.write_text("preserve", encoding="utf-8")
            result = subprocess.run([sys.executable, str(Path(__file__).with_name("gwas_summary_qc.py")), "--summary", str(paths[0]), "--metadata", str(paths[1]), "--out-dir", str(output)], capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertEqual(report.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main()
