from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures" / "2026_methods"
SCRIPTS = {
    "598": REPO / "modules" / "12_genome_annotation_variants" / "01_somatic_variants" / "598_clairs_longread_somatic_audit" / "598_clairs_longread_somatic_audit.py",
    "599": REPO / "modules" / "12_genome_annotation_variants" / "02_genome_annotation" / "599_eviann_annotation_audit" / "599_eviann_annotation_audit.py",
    "600": REPO / "modules" / "13_comparative_genomics_phylogenetics" / "01_orthology" / "600_orthofinder_v3_audit" / "600_orthofinder_v3_audit.py",
    "601": REPO / "modules" / "14_microbiome_metagenomics" / "01_abundance" / "601_coverm_abundance_audit" / "601_coverm_abundance_audit.py",
    "602": REPO / "modules" / "15_epigenomics_chromatin" / "01_damid" / "602_damidbind_result_audit" / "602_damidbind_result_audit.py",
    "603": REPO / "modules" / "16_proteomics_metabolomics" / "01_dia_proteomics" / "603_diareport_result_audit" / "603_diareport_result_audit.py",
}


class MethodModuleTests(unittest.TestCase):
    maxDiff = None

    def run_module(self, module_id: str, fixture: str, *extra: str):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        outdir = Path(temp.name) / "results"
        process = subprocess.run(
            [sys.executable, str(SCRIPTS[module_id]), "--input", str(FIXTURES / fixture), "--outdir", str(outdir), *extra],
            cwd=REPO,
            text=True,
            capture_output=True,
        )
        return process, outdir

    def read_summary(self, module_id: str, outdir: Path):
        with (outdir / f"{module_id}_summary.json").open(encoding="utf-8") as handle:
            return json.load(handle)

    def assert_figure_pair(self, module_id: str, outdir: Path):
        self.assertTrue((outdir / f"{module_id}_audit.png").is_file())
        self.assertTrue((outdir / f"{module_id}_audit.pdf").is_file())

    def test_598_counts_filter_and_variant_types(self):
        process, outdir = self.run_module("598", "clairs_calls.vcf")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            self.read_summary("598", outdir),
            {"filtered_variants": 1, "indel_variants": 1, "other_variants": 1, "pass_variants": 3, "snv_variants": 2, "total_variants": 4},
        )
        self.assert_figure_pair("598", outdir)

    def test_598_rejects_header_only_vcf(self):
        process, _ = self.run_module("598", "empty.vcf")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("no variant records", process.stderr.lower())

    def test_598_emits_verified_clairs_command(self):
        process = subprocess.run(
            [sys.executable, str(SCRIPTS["598"]), "--print-upstream-command", "--tumor-bam", "tumor.bam", "--normal-bam", "normal.bam", "--reference", "ref.fa", "--platform", "ont_r10_guppy", "--threads", "8", "--upstream-output", "output"],
            cwd=REPO, text=True, capture_output=True,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(process.stdout.strip(), "run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10_guppy")

    def test_599_summarizes_features_and_evidence(self):
        process, outdir = self.run_module("599", "eviann_annotation.gff3")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            self.read_summary("599", outdir),
            {"cds_features": 2, "complete_transcripts": 1, "exon_features": 3, "gene_features": 2, "pseudo_transcripts": 1, "transcript_features": 3, "transcripts_with_evidence": 2},
        )
        self.assert_figure_pair("599", outdir)

    def test_599_rejects_malformed_gff3(self):
        process, _ = self.run_module("599", "malformed.gff3")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("9 tab-separated fields", process.stderr)

    def test_599_emits_verified_eviann_command(self):
        process = subprocess.run(
            [sys.executable, str(SCRIPTS["599"]), "--print-upstream-command", "--genome", "genome.fasta", "--reads-list", "paired.txt", "--proteins", "proteins.faa", "--threads", "24"],
            cwd=REPO, text=True, capture_output=True,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(process.stdout.strip(), "eviann.sh -t 24 -g genome.fasta -r paired.txt -p proteins.faa")

    def test_600_reports_orthogroup_coverage_and_single_copy(self):
        process, outdir = self.run_module("600", "Orthogroups.tsv")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            self.read_summary("600", outdir),
            {"complete_orthogroups": 1, "orthogroups": 3, "single_copy_complete_orthogroups": 1, "species": 3, "total_gene_assignments": 7},
        )
        self.assert_figure_pair("600", outdir)

    def test_600_rejects_table_without_species_columns(self):
        process, _ = self.run_module("600", "orthogroups_no_species.tsv")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("species column", process.stderr.lower())

    def test_600_emits_verified_orthofinder_command(self):
        process = subprocess.run(
            [sys.executable, str(SCRIPTS["600"]), "--print-upstream-command", "--proteomes", "proteomes", "--search-threads", "8", "--analysis-threads", "4"],
            cwd=REPO, text=True, capture_output=True,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(process.stdout.strip(), "orthofinder -f proteomes -t 8 -a 4")

    def test_601_checks_abundance_totals(self):
        process, outdir = self.run_module("601", "coverm_abundance.csv")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            self.read_summary("601", outdir),
            {"genomes": 2, "max_relative_abundance": 0.7, "samples": 2, "samples_summing_to_one": 1, "total_rows": 4},
        )
        self.assert_figure_pair("601", outdir)

    def test_601_rejects_negative_abundance(self):
        process, _ = self.run_module("601", "coverm_negative.csv")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("nonnegative", process.stderr.lower())

    def test_601_emits_verified_coverm_command(self):
        process = subprocess.run(
            [sys.executable, str(SCRIPTS["601"]), "--print-upstream-command", "--bam", "sample.bam", "--genome-definition", "genome-definition.tsv"],
            cwd=REPO, text=True, capture_output=True,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(process.stdout.strip(), "coverm genome --bam-files sample.bam --genome-definition genome-definition.tsv --methods relative_abundance")

    def test_602_classifies_differential_binding(self):
        process, outdir = self.run_module("602", "damidbind_results.csv")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            self.read_summary("602", outdir),
            {"fdr_threshold": 0.05, "gain_loci": 1, "loss_loci": 1, "significant_loci": 2, "tested_loci": 4},
        )
        self.assert_figure_pair("602", outdir)

    def test_602_rejects_missing_padj(self):
        process, _ = self.run_module("602", "damidbind_missing_padj.csv")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("padj", process.stderr.lower())

    def test_603_classifies_differential_proteins(self):
        process, outdir = self.run_module("603", "diareport_results.csv")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            self.read_summary("603", outdir),
            {"downregulated_proteins": 1, "fdr_threshold": 0.05, "significant_proteins": 2, "tested_proteins": 4, "upregulated_proteins": 1},
        )
        self.assert_figure_pair("603", outdir)

    def test_603_rejects_missing_log2fc(self):
        process, _ = self.run_module("603", "diareport_missing_log2fc.csv")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("log2fc", process.stderr.lower())


if __name__ == "__main__":
    unittest.main()
