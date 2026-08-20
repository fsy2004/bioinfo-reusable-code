from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.update_readme import compute_metrics, local_link_errors, replace_block


class ReadmeUpdateTests(unittest.TestCase):
    def test_compute_metrics_matches_catalog_entries_for_numbered_files_and_directories(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "modules" / "01_alpha" / "01_one" / "001_first").mkdir(parents=True)
            (root / "modules" / "01_alpha" / "01_one" / "notes").mkdir()
            (root / "modules" / "02_beta" / "01_two" / "002_second").mkdir(parents=True)
            (root / "modules" / "02_beta" / "01_two" / "003_script.py").write_text("", encoding="utf-8")
            (root / "modules" / "02_beta" / "01_two" / "004legacy.R").write_text("", encoding="utf-8")
            (root / "modules" / "02_beta" / "01_two" / "resource.csv").write_text("", encoding="utf-8")
            (root / "modules" / "_framework").mkdir()
            self.assertEqual(
                compute_metrics(root),
                {"domains": 2, "subcategories": 2, "modules": 4, "catalog_entries": 6},
            )

    def test_replace_block_updates_only_named_generated_region(self):
        source = "before\n<!-- readme-metrics:start -->\nold\n<!-- readme-metrics:end -->\nafter\n"
        self.assertEqual(
            replace_block(source, "readme-metrics", "new"),
            "before\n<!-- readme-metrics:start -->\nnew\n<!-- readme-metrics:end -->\nafter\n",
        )

    def test_local_link_errors_reports_missing_repo_paths_only(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "docs").mkdir()
            (root / "docs" / "ok.md").write_text("ok", encoding="utf-8")
            text = "[ok](docs/ok.md) [missing](docs/missing.md) [web](https://example.org) [anchor](#quick-start)"
            self.assertEqual(local_link_errors(root, text), ["docs/missing.md"])


if __name__ == "__main__":
    unittest.main()
