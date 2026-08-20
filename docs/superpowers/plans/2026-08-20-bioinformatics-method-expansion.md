# 2026 Bioinformatics Method Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use executing-plans to implement this plan task-by-task.

**Goal:** Add six source-audited, turnkey output-audit modules and five missing top-level bioinformatics domains without installing heavy upstream tools or changing remote Git state.

**Architecture:** Each module is an independent Python CLI that validates one stable upstream output format, computes transparent descriptive QC, writes CSV/JSON, and saves PNG/PDF. A frozen method-scan registry records broader 2026 coverage and explains why candidates were selected, deferred, or excluded. Catalog metadata remains generated from `gen_catalog.py` plus `new_modules.json`.

**Tech Stack:** Python 3, standard library, pandas, matplotlib, unittest; existing catalog and QC scripts.

**Spec:** `docs/specs/2026-bioinformatics-method-expansion.md`

## Global Constraints

- Preserve the clean `master` worktree and unrelated files; no commit, push, PR, install, or large download.
- Source every method and command from an official paper/repository.
- Default runs use bundled synthetic files only and never claim full upstream method execution.
- Write behavioral tests before production scripts, observe the expected RED state, then implement minimally.

---

### Task 1: Lock behavior with failing integration tests

**Files:**
- Create: `tests/test_2026_method_modules.py`
- Create: `tests/fixtures/2026_methods/*`

1. Add one success test per module using hand-checked fixture expectations.
2. Add malformed-input tests for missing required columns or invalid values.
3. Add exact command-contract tests for ClairS, EviAnn, OrthoFinder and CoverM.
4. Run `python -m unittest tests.test_2026_method_modules -v`; confirm failure because production scripts do not exist.

### Task 2: Implement shared audit utilities and modules 598–600

**Files:**
- Create: `modules/_framework/output_audit.py`
- Create: `modules/12_genome_annotation_variants/01_somatic_variants/598_clairs_longread_somatic_audit/*`
- Create: `modules/12_genome_annotation_variants/02_genome_annotation/599_eviann_annotation_audit/*`
- Create: `modules/13_comparative_genomics_phylogenetics/01_orthology/600_orthofinder_v3_audit/*`

1. Implement stable CSV/JSON writing and paired PNG/PDF saving.
2. Parse VCF and report PASS/filter and SNV/indel counts; print the official ClairS CLI when requested.
3. Parse GFF3 and report feature hierarchy/evidence completeness; print the official EviAnn CLI when requested.
4. Parse Orthogroups.tsv and report per-group and overall coverage/single-copy metrics; print the official OrthoFinder CLI when requested.
5. Run only tests 598–600 until green.

### Task 3: Implement modules 601–603

**Files:**
- Create: `modules/14_microbiome_metagenomics/01_abundance/601_coverm_abundance_audit/*`
- Create: `modules/15_epigenomics_chromatin/01_damid/602_damidbind_result_audit/*`
- Create: `modules/16_proteomics_metabolomics/01_dia_proteomics/603_diareport_result_audit/*`

1. Validate CoverM long-table abundance, sample totals and nonnegative measures; print official `coverm genome` CLI.
2. Validate damidBind exported differential results and classify significant gains/losses.
3. Validate DiaReport/MSqRob exported differential results and classify up/down proteins.
4. Run tests 601–603 until green, then run the complete test file.

### Task 4: Add documentation, examples, and the frozen method registry

**Files:**
- Create six module `README.md` files and six `example_data` inputs.
- Create: `modules/_framework/METHOD_SCAN_2026_08_20.md`
- Modify: `modules/_framework/TOOL_SELECTION_GUIDE.md`
- Modify: `modules/_framework/SERVER_DEPENDENCIES.md`

1. Document input schema, method boundary, uses, strengths, outputs, upstream commands, dependencies, citation, license, and interpretation limits.
2. Record the 2026 search window, domain coverage, selected methods and deferred candidates with primary URLs/DOIs.
3. Run every module with defaults and copy generated preview PNGs to `assets/`.

### Task 5: Register domains and regenerate catalogs

**Files:**
- Modify: `modules/_framework/gen_catalog.py`
- Modify: `modules/_framework/new_modules.json`
- Regenerate: `modules/CATALOG.md`, domain `README.md` files, repository `README.md`

1. Add domains 12–16 and their subcategories to the hard-coded domain map.
2. Append records 598–603 with honest yellow status and verified citations.
3. Run `python modules/_framework/gen_catalog.py` and the README regeneration command documented in the script/repository.

### Task 6: Final verification

1. Run `python -m unittest tests.test_2026_method_modules -v`.
2. Run each module with its default input in a clean temporary output directory.
3. Run `python modules/_framework/qc_lint.py` against all six module directories.
4. Run catalog consistency checks and inspect `git diff --check`, `git status`, generated summaries and representative PNGs.
5. Report exact changes, verification results, uncommitted status, and any remaining registry-only methods.
