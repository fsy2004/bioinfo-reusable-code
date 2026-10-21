"""User-facing method retrieval/compatibility and source-archive safety contracts."""
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "modules/04_causal_inference_genetics/09_method_library/611_mr_coloc_library/method_library.py"
spec = importlib.util.spec_from_file_location("method_library", MODULE)
lib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)


def example_registry():
    def row(identifier, category, requirements):
        return {"id": identifier, "name": identifier, "category": category, "family": identifier,
                "source_ids": ["SRC_TEST"], "assumptions": ["Read native method assumptions"],
                "selectable": True, "requirements": requirements}
    return {"sources": [{"id": "SRC_TEST"}], "methods": [
        row("MR-WALD", "mr_estimator", {"min_snps": 1, "exact_snps": 1, "independent_iv": True}),
        row("MR-EGGER", "mr_estimator", {"min_snps": 3, "independent_iv": True}),
        row("COLOC-ABF", "colocalization", {"min_snps": 3, "dense_region": True, "min_traits": 2, "single_causal_variant": True}),
        row("COLOC-SUSIE", "colocalization", {"min_snps": 3, "dense_region": True, "min_traits": 2, "signed_ld": True}),
    ]}


class SelectionTests(unittest.TestCase):
    def test_single_snp_does_not_recommend_egger(self):
        result = lib.select_methods(example_registry(), {"n_snps": 1, "independent_instruments": True})
        self.assertEqual([r["id"] for r in result["compatible_for_review"]], ["MR-WALD"])
        self.assertFalse(result["scientific_validity_checked"])

    def test_r_squared_and_significant_only_are_not_sufficient_for_susie_coloc(self):
        result = lib.select_methods(example_registry(), {"n_snps": 20, "category": "colocalization", "ld": "r2", "dense_region": False, "n_traits": 2})
        self.assertEqual(result["compatible_for_review"], [])
        susie = next(r for r in result["blocked"] if r["id"] == "COLOC-SUSIE")
        self.assertIn("requires_signed_r_not_r2", susie["reasons"])
        self.assertIn("requires_dense_region_not_clumped_or_significant_only", susie["reasons"])

    def test_multiple_signals_redirects_away_from_single_signal_abf(self):
        result = lib.select_methods(example_registry(), {"n_snps": 40, "category": "colocalization", "ld": "signed_r", "dense_region": True, "n_traits": 2, "multiple_signals": True})
        self.assertEqual([r["id"] for r in result["compatible_for_review"]], ["COLOC-SUSIE"])

    def test_query_does_not_confuse_coloc_with_mr(self):
        rows = lib.query(example_registry(), "", category="colocalization")
        self.assertEqual({r["id"] for r in rows}, {"COLOC-ABF", "COLOC-SUSIE"})

    def test_duplicate_or_unresolved_ids_fail(self):
        data = example_registry(); data["methods"].append(data["methods"][0].copy())
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            lib.validate_registry(data)

    def test_string_booleans_cannot_assert_instrument_independence(self):
        with self.assertRaisesRegex(ValueError, "JSON boolean"):
            lib.select_methods(example_registry(), {"n_snps": 3, "independent_instruments": "false"})
        with self.assertRaisesRegex(ValueError, "n_traits"):
            lib.select_methods(example_registry(), {"n_snps": 3, "n_traits": True})

    def test_source_id_cannot_escape_cache(self):
        data = example_registry(); data["sources"][0]["id"] = "../outside"
        with self.assertRaisesRegex(ValueError, "Unsafe source"):
            lib.validate_registry(data)
        data = example_registry(); data["sources"] = []
        with self.assertRaisesRegex(ValueError, "Unresolved"):
            lib.validate_registry(data)


class ArchiveTests(unittest.TestCase):
    def test_fresh_checkout_uses_portable_snapshot_without_local_source_cache(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch.object(lib, "CACHE", Path(temp)):
                report = lib.build_indexes(lib.load_registry(), check=True)
            self.assertGreater(report["official_sources"], 0)
            self.assertTrue(report["runtime_or_scientific_validation_not_inferred_from_source_count"])

    def test_archive_fetch_function_rejects_dot_source_id_before_filesystem_write(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp).resolve() / "cache"
            for identifier in [".", ".."]:
                with self.assertRaisesRegex(ValueError, "Unsafe source ID"):
                    lib.fetch_source({"id": identifier, "repository": "test/repo"}, cache)
            self.assertFalse(cache.exists())

    def test_cache_identity_rejects_changed_repository_pin_and_manifest(self):
        receipt = {"source_id": "SRC_TEST", "repository": "old/repo", "requested_repository": "old/repo", "commit": "a" * 40}
        manifest = {"repository": "old/repo", "commit": "a" * 40}
        for source in [{"id": "SRC_TEST", "repository": "new/repo"}, {"id": "SRC_TEST", "repository": "old/repo", "pinned_commit": "b" * 40}]:
            with self.assertRaisesRegex(ValueError, "Cached"):
                lib.validate_cache_identity(source, receipt, manifest)
        with self.assertRaisesRegex(ValueError, "manifest identity"):
            lib.validate_cache_identity({"id": "SRC_TEST"}, receipt, {**manifest, "commit": "b" * 40})
        lib.validate_cache_identity({"id": "SRC_TEST", "repository": "OLD/REPO", "pinned_commit": "a" * 40}, receipt, manifest)

    def test_windows_drive_or_backslash_path_is_rejected(self):
        for name in ["repo/C:/escape.R", "repo/..\\escape.R"]:
            with tempfile.TemporaryDirectory() as temp:
                root = Path(temp).resolve(); archive = root / "bad.tgz"
                with tarfile.open(archive, "w:gz") as handle:
                    entry = tarfile.TarInfo(name); entry.size = 3
                    handle.addfile(entry, io.BytesIO(b"bad"))
                with self.assertRaisesRegex(ValueError, "Unsafe"):
                    lib.safe_extract(archive, root / "target")

    def test_integrity_check_detects_modified_source_and_unreceipted_files(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp).resolve(); directory = cache / "SRC_TEST"; tree = directory / "commit-code"
            tree.mkdir(parents=True); code = tree / "method.R"; code.write_bytes(b"x<-1")
            lib.write_json(directory / "source_files.json", {"files": [{"path": "method.R", "bytes": 4, "sha256": lib.hashlib.sha256(b"x<-1").hexdigest()}]})
            lib.write_json(directory / "receipt.json", {"status": "source_subset_acquired", "source_tree": str(tree), "file_manifest": "source_files.json"})
            self.assertEqual(lib.verify_source({"id": "SRC_TEST"}, cache)["status"], "verified")
            code.write_bytes(b"x<-2"); (tree / "extra.R").write_bytes(b"unsafe")
            check = lib.verify_source({"id": "SRC_TEST"}, cache)
            self.assertEqual(check["status"], "integrity_failed")
            self.assertIn("method.R", check["failures"])
            self.assertIn("unreceipted_file:extra.R", check["failures"])

    def test_traversal_cannot_write_outside_cache(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve(); target = root / "source"; archive = root / "bad.tgz"
            self.assertTrue(target.resolve().is_relative_to(root))
            with tarfile.open(archive, "w:gz") as handle:
                entry = tarfile.TarInfo("repo/../../escape.txt"); entry.size = 3
                handle.addfile(entry, io.BytesIO(b"bad"))
            with self.assertRaisesRegex(ValueError, "Unsafe"):
                lib.safe_extract(archive, target)
            self.assertFalse((root.parent / "escape.txt").exists())

    def test_regular_code_keeps_identity_but_symlink_is_not_followed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve(); target = root / "source"; archive = root / "good.tgz"
            self.assertTrue(target.resolve().is_relative_to(root))
            with tarfile.open(archive, "w:gz") as handle:
                entry = tarfile.TarInfo("repo/R/method.R"); entry.size = 4
                handle.addfile(entry, io.BytesIO(b"x<-1"))
                link = tarfile.TarInfo("repo/link"); link.type = tarfile.SYMTYPE; link.linkname = "../../outside"
                handle.addfile(link)
            files = lib.safe_extract(archive, target)
            self.assertEqual([r["path"] for r in files], ["R/method.R"])
            self.assertEqual((target / "R/method.R").read_bytes(), b"x<-1")
            self.assertFalse((target / "link").exists())


if __name__ == "__main__":
    unittest.main()
