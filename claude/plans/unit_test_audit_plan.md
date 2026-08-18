# Unit test audit — tests that don't earn their keep

Audit against the CLAUDE.md rule:

> A test earns its keep when it covers logic *we* wrote: a branch, a fallback, a calculation,
> a rule that's easy to get wrong later. […] If the test would still pass with our logic
> deleted, or it only restates a field declaration, drop it.

## Scope and method

- 250 project test files, ~32,100 lines, 1,997 test methods (excludes `.venv/` and migrations).
- Read closely: every file added since 2026-01-01 (163 files, ~23,600 lines — the recent,
  largely AI-written cohort) plus every pre-2026 file under ~130 lines.
- Large legacy files (`test_annotation_dispatch.py`, `test_external_annotation.py`,
  `test_annotation_vcf.py`, `test_scheduler.py`, `test_joint_called_vcf.py`, …) were sampled
  and swept with AST/grep rather than read line by line.
- Whole-suite sweeps: assertion counts per test, tests with zero assertions, duplicate test
  method names across files, `@skip` decorators, `self.fail()` in place of `assertRaises`,
  stale `BUG:` / `TODO` / `print()` markers, function-level imports.

Headline: the suite is in better shape than expected. Most of what came out of the recent
test-writing push tests real branches. The problems cluster into five recurring shapes, listed
in priority order below.

---

## Tier 1 — Delete: tests that never run or test nothing

### 1.1 `classification/tests/utils/test_evidence_mixin.py` — whole test class commented out

The only live code is a `BasicEvidence` helper class; the `EvidenceMixinTest` below it is
commented out with `# doesn't work without Transcripts loaded now`. `BasicEvidence` is
referenced nowhere else in the codebase (verified). Delete the file.

### 1.2 `annotation/tests/test_bugs.py:132` — `test_hpo_prefix_slice_length_is_correct`

```python
self.assertEqual(len("Human_Phenotype_Ontology:"), 25)
```

This asserts a Python string literal's length. It never touches `ClinVar._database_terms`, so it
would pass with our logic deleted. The *intent* (catch the hardcoded `name[25:]` slice going
stale) is real — the production code (`annotation/models/models.py:243-244`) is even worse than
the test implies: the `startswith` check uses the 24-char `"Human_Phenotype_Ontology"` (no
colon) while the slice is 25, so the two constants are already inconsistent. The fix is in
production code: derive the offset with `len(PREFIX)` / `removeprefix()` like the neighbouring
Orphanet branch does, at which point the test is moot. `test_hpo_prefix_stripped` (line 128)
already covers the behaviour.

### 1.3 `library/tests/test_utils.py:434` — `TestSha256SumStr`

`sha256sum_str` is a two-line wrapper over `hashlib.sha256`. Both tests hard-code known SHA-256
digests, i.e. they test hashlib. Drop the class. (Keep `TestStringDeterministicHash` — that one
*is* our own accumulator loop — but `test_order_matters` at line 429 asserts only that two
different strings hash differently, which any hash satisfies; fold it away.)

---

## Tier 2 — Tests that pin bugs in place (actively harmful)

These pass today *because* the code is wrong, and will fail when someone fixes it. That inverts
what a test is for: the next person to do the right thing gets a red build.

### 2.1 `annotation/tests/test_bugs.py:85` — `test_v2_old_r_prefix_crashes`

```python
def test_v2_old_r_prefix_crashes(self):
    # Old GRCh37 gnomADv2 VEP output stores "r2.1"; int("r2") → ValueError.
    # This will start passing once the property handles the "r" prefix.
    with self.assertRaises(ValueError):
        self._call("r2.1")
```

The property itself is verified unchanged: `annotation/models/models.py:815` is
`int(self.gnomad.split(".", maxsplit=1)[0])`, with no `r` handling (and a `None` value would
raise `AttributeError`). But the severity is lower than the test's comment implies:
`annotation/migrations/0042_one_off_variant_annotation_version_gnomad.py` already rewrote the
legacy `'r2.1'` values to `'2.1.1'`/`'3.1'`, and `parse_gnomad_version_from_filename`
(`annotation/vep_annotation.py:452`) only produces numeric strings — so migrated deployments
never hit this. Dependents are `gnomad4_or_later` and `has_hemi` (there is no `has_mid`).

Since the crash is unreachable post-migration, lean delete: drop the test, or fix the property
defensively (strip a leading `r`) and flip the assertion to `== 2` in the same change. Either
way a passing test suite should not be the place a "known crash" is documented.

### 2.2 `genes/tests/test_genes_adversarial.py:555` — `test_exon_id_starts_at_1_skips_start_check`

```python
# Currently: no errors (false negative) — document actual behaviour
self.assertEqual(errors, [], "Exon_id=1 silently skips the cdna_start==1 check")
```

This locks a validation false-negative into the test suite: `_validate_cdna_match` misses a
transcript whose first exon starts at cDNA 5 whenever exon numbering is 1-based. Either fix the
`if exon_id == 0:` check and assert the error is raised, or delete the test. Keeping it means a
future correct fix arrives as a test failure.

### 2.3 `annotation/tests/test_bugs.py:149` — `test_99_6_percent_truncates_to_99_not_100`

The load-bearing assertion is `self.assertEqual(int(of * 100), 99)` — the test performs the
truncation itself rather than calling the annotation code that does it, so it would pass with
our logic deleted. Keep the first line (`assertAlmostEqual(of, 0.996)` — that *is*
`overlap_fraction`) and drop the rest, or rewrite it to go through the real
`int(of * 100)` call site so it actually guards the SV-overlap column.

---

## Tier 3 — Stale "BUG:" narratives: the bug was fixed, the story wasn't

A large cluster of tests carry docstrings and names asserting the code is broken. The code has
since been fixed (verified against source in each case below), so the comments now describe
behaviour that does not exist, several test names say the opposite of what they assert, and a
reader can't tell which "BUG:" markers are live. Keep the tests — most cover genuine
falsy-vs-`None` branches — but rewrite them as plain regression tests.

### 3.1 `patients/tests/test_patients_bugs.py` — whole-file header is wrong

File docstring: *"Tests that expose confirmed bugs are expected to FAIL until the bug is fixed."*
All the named bugs are fixed:

| Test docstring claims | Actual source |
|---|---|
| `_age_at_collection_date` 0 treated as missing | `patients/models.py:381` — `if self._age_at_collection_date is not None:` |
| `condition_description` ignores `_deceased=True` | `patients/models.py:217` — `if self._deceased: return "dead"` |
| …plus the `parse_boolean`, `process_record` and `PatientForm` items | all guarded |

Actions: rewrite the module docstring and per-test `BUG:` blocks as "regression: X used to Y";
rename the file to something describing the subject rather than `_bugs`; replace
`assertRaises(Exception)` (line 199) with the specific exception (production now raises
`MultipleObjectsReturned`) — as written it also passes on an `AttributeError`; and move the
function-level `from patients.forms import PatientForm`
(line 316) to the top of the file per the CLAUDE.md import rule.

### 3.2 `pedigree/tests/test_ped_utils.py` — three test names assert their own opposite

```python
def test_lowercase_m_not_recognised(self):        # line 29
    # BUG-4: 'm' is not in SEX_LOOKUP — silently returns None instead of MALE
    self.assertEqual(get_sex('m'), Sex.MALE, "…currently returns None (bug)")
```

`get_sex` now does `SEX_LOOKUP.get(str(sex).upper())` and `get_parent_id` handles NaN
explicitly. Rename `test_lowercase_m_not_recognised` → `test_lowercase_m_is_male` (same for
`_f_`, line 35). `test_nan_is_unknown` (line 52) is accurately named — only its stale `BUG-6`
comment needs to go. `test_unknown_value_returns_none_not_raises` restates `dict.get`'s
default — drop it and keep the asymmetry note (it lives on `TestGetSex`, lines 25–26), which
is the half that's actually a rule we chose.

### 3.3 `analysis/tests/test_node_filter_bugs.py` — class docstring quotes deleted source

`TestPopulationNodeModifiesParents` quotes `population_node.py:60-61` verbatim with
`# BUG: 0 is falsy`. The real line 61 is now `self.gnomad_hom_alt_max is not None`. Quoting
production source with line numbers into a docstring guarantees this rots — replace with one
sentence naming the rule (`gnomad_hom_alt_max=0 is a real filter, not "unset"`).

### 3.4 `library/tests/test_utils.py` — three stale bug annotations

- `test_max_lines_off_by_one` (206) and `test_max_lines_zero` (211) say `# Bug: 'i > max_lines'
  returns max_lines+1 lines`; `library/utils/file_utils.py:49` is `if … i >= max_lines`. Fixed.
- `test_unequal_splits_silent_data_loss` (109) is named for silent data loss, but
  `split_dict_multi_values` raises `ValueError` (`text_utils.py:45`). Rename to
  `test_unequal_splits_raise` and use `assertRaises` instead of the `try/self.fail/except` shape.
- `test_ttl_stale_value_returned_on_expiry` (462) is named for the stale value but asserts the
  function *is* re-called; it also `time.sleep(0.12)`s. Rename to `test_ttl_expiry_recalculates`.

### 3.5 `upload/tests/test_modified_imported_variant.py:33`

Docstring: *"This test is expected to FAIL until the code adds a guard"*. The guard exists at
`upload/models/models.py:770`. Delete the two stale paragraphs; the test itself is good.

### 3.6 `genes/tests/test_genes_adversarial.py` — "Bug candidate" comments that were fixed

`test_three_letter_lowercase_input_normalized` (57, handled by the `.capitalize()` fallback in
`genes/hgvs/phgvs.py:47`), `test_tag_with_space_after_comma_not_matched` (594, name says "not
matched" but asserts it *is*; `.strip()` in `genes/models.py:1116` handles it), and
`test_empty_genome_build_tag_falls_back_to_data` (610 — here nothing was "fixed": the
`or`-fallthrough is the current intended behaviour, just mislabelled as a bug candidate). Same
treatment for all three: rename to the behaviour, drop the speculation.

### 3.7 `analysis/tests/test_explicit_pk_substitution_baseline.py` — no longer a baseline

The docstring frames the file as a *pre-change* snapshot to be re-run "after the optimisation
lands". Issue #546 has landed (`ANALYSIS_NODE_STORE_ID_SIZE_MAX` is live at
`analysis/models/nodes/analysis_node.py:605`). The gate-on-vs-off tests are still valuable —
the gate is a runtime setting, so both paths ship — but the file should be renamed
(`test_explicit_pk_substitution.py`) and the docstring rewritten as "both gate settings must
produce identical result sets".

---

## Tier 4 — Trims inside otherwise-good files

Individually small; together they are the "restates a field declaration" category.

| File:line | Test | Why it doesn't earn its keep |
|---|---|---|
| `snpdb/tests/test_data_archive_mixin.py:38` | `test_mixin_fields_present` | Asserts `hasattr` / `_meta.get_fields()` contains the four mixin field names on VCF/VAV/GCC. Pure field-declaration restatement; every other test in the file exercises those fields for real. Delete. |
| `snpdb/tests/test_data_archive_mixin.py:45` | `test_data_archived_property_reflects_date` | `data_archived` is `bool(data_archived_date)`. Covered incidentally by `test_force_archive_when_uploaded_file_missing`. Delete. |
| `snpdb/tests/test_node_grid_auto_load_setting.py:65` | `test_falls_back_to_settings_constant_when_all_null` | The last three lines re-implement the view's `… or settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS` *inside the test* and assert on that expression. Keep the `assertIsNone(self._resolved())` half (that is the cascade returning nothing); delete the rest, or reach the fallback through the view. |
| `snpdb/tests/test_vcf_export_columns.py:22` | `test_types_are_valid` | Asserts each `c.type` is a `VCFInfoTypes` — guaranteed by the literals in `COLUMN_VCF_INFO`. The sibling uniqueness tests do earn their keep (they preserve constraints dropped with the `ColumnVCFInfo` table). Delete this one. |
| `library/tests/test_preview_request.py:24` | `TestPreviewDataHash` (4 tests) | `PreviewData.__hash__` *is* custom (`library/preview_request.py:312`), so the class should stay — but `test_hash_is_stable` and `test_hash_usable_in_set` are the same identity check twice, and none of the four touch `genome_builds` / `annotation_consortia` / `summary_extra`, the three collection fields the custom hash exists for. Collapse to two tests and cover the collections. |
| `library/tests/test_utils.py:242` | `test_read_exhausted_returns_empty_or_none` | `assertIn(result, ['', None])` accepts two answers, so it pins nothing. The contract is `''` — `IteratorFile.read` (`library/utils/file_utils.py:83-85`) can never return `None` from `read()` — so assert `== ''`. |
| `analysis/tests/test_variant_tags.py` | `test_allele_not_yet_assigned` | Asserts a field we never set is `None`. Delete — the "matched on variant, not allele" rule is proven by the four sibling tests. |
| `snpdb/tests/test_sample_formatter.py` | `test_patient_code_param_blank_when_unset` | **Keep** (corrected on re-verification): the `""` is not Django field behaviour — it comes from our own `patient.patient_code or ""` coalesce at `snpdb/models/models_vcf.py:554`; without it the format string renders `"None"`. Both tests earn their keep. |
| `annotation/tests/test_pathogenicity_predictions.py` | `TestPathogenicityThresholdsTag` | The calibrated/uncalibrated split is our logic (keep `test_excludes_uncalibrated_tools`), but the two band tests hard-code 0.290/0.644 and 22.7/25.3, restating declarations. Also duplicated wholesale by `variantopedia/tests/test_pathogenicity_thresholds_template.py`, which asserts the same three things through a rendered template. Keep one. |
| `classification/tests/utils/test_json_utils.py:24` | `test_json_diff` | Ends with a `JsonDiffs.differences(["A","b"], ["a","b","c"])` call followed by `print(diffs)` and **no assertion** — dead tail. Delete the tail (and the earlier `print`). |
| `ontology/tests/test_ontology.py:13` | `testLoadData` | Whole body is `create_ontology_test_data()`. It's a smoke test that the fixture loaders don't raise — legitimate, but say so in the name (`test_ontology_import_loaders_run`) and drop the no-op `setUpTestData` and `if __name__ == "__main__"` block. |
| `analysis/tests/test_zygosity_nodes.py:41` | `_test_zygosity` | Sets eight zygosity bounds then asserts only `assertIsInstance(count, int)` — it proves the SQL is valid, never that the filters filter. Assert an expected count against the fake trio. |
| `upload/tests/test_uploaded_file.py` | 3 permission tests | `FileUpload.can_view` is ours (superuser-or-owner, not Guardian), so these are defensible, but three tests for a three-line method is heavy. Low priority — fold into one. |

---

## Tier 5 — Structural duplication

### 5.1 `analysis/tests/test_quad_node.py` ↔ `analysis/tests/test_trio_node.py`

`TestQuadNodeInheritance` and `TestTrioNodeInheritance` share **24 identically-named tests**
(`test_all_recessive_*`, `test_denovo_*`, `test_xlinked_*`, `test_zygosity_table_*`, …) across
745 lines, plus several near-duplicates under renamed variants (`require_zygosity` vs
`require_parent_zygosity`), so the real overlap is higher still. The trio file also has ~15
tests the quad file lacks (`test_dominant_matches_*`, `test_only_compound_het_takes_a_parent`,
…). The genuinely quad-specific coverage is the sibling constraint; everything else is
the same assertions over a 4-sample instead of 3-sample zygosity string. Extract a shared
inheritance-mode mixin parameterised on the sample layout, and leave only the sibling tests in
the quad file.

### 5.2 `snpdb/tests/test_library_utils.py::test_sig_digits` ↔ `library/tests/test_utils.py::TestFormatSignificantDigits`

Same function, two files, overlapping cases. `snpdb/tests/test_library_utils.py`'s docstring
("Django unit tests are run in apps so this needs to be here") predates `library/tests/`
existing. Move `test_markdown` and `test_utc_from_timestamp` into `library/tests/`, merge the
sig-digits cases, and delete the snpdb file.

### 5.3 `genes/tests/test_clean_hgvs.py::test_clean_hgvs` ↔ `genes/tests/test_hgvs.py::TestHGVS.test_clean_hgvs`

Same name, same function under test. Consolidate on `test_clean_hgvs.py` (which has the better
framing — it documents which invariants are ours vs cdot's). The `test_hgvs.py` copy is the
weaker one by more than framing: it asserts nothing — it `print()`s each cleaned result and
relies on `create_hgvs_variant` not raising — and carries a dead commented-out
try/`self.fail` block at lines 61–65.

### 5.4 `snpdb/tests/test_cohort.py:30` — `test_not_sub_cohort` duplicates `test_sub_cohort`'s setup

(Corrected on re-verification: an earlier draft claimed this test asserts the opposite of its
name — it doesn't. After the `assertTrue(cohort.is_sub_cohort)` line it goes on to add a
sample from `trio2` and asserts `assertFalse(cohort.is_sub_cohort)`, so it does test what its
name says.) The remaining issue is only that its first half replays `test_sub_cohort`'s body
and re-asserts its result. Trim the duplicated prelude, or fold both into one test.

---

## Hygiene sweep (cheap, mechanical)

- **Function-level imports** — 15 across the suite (list verified complete by an independent
  sweep; `test_annotation_vcf.py:464-465` is two imports), against the CLAUDE.md rule that all
  imports go at the top:
  `analysis/tests/test_scheduler.py:500`, `annotation/tests/test_annotation_dispatch.py:334`,
  `annotation/tests/test_annotation_vcf.py:464-465`, `annotation/tests/test_phenotype_matching.py:93,147`,
  `annotation/tests/test_spliceai.py:89,176,211`, `genes/tests/test_gene_coverage_archive.py:34,40,56`,
  `patients/tests/test_patients_bugs.py:316`, `snpdb/tests/test_partition_archive_model.py:48`,
  `snpdb/tests/test_partition_archive_task.py:106`.
- **`print()` left in tests** — noise in every CI run:
  `analysis/tests/test_clone_nodes.py:86`, `annotation/tests/test_clinvar_xml_parser.py:27,36`,
  `annotation/tests/test_annotation_vcf.py:64,86-88,131-133`,
  `annotation/tests/test_annotation_vcf_cnv.py:55,65-67,96-98`,
  `classification/tests/utils/test_clinvar_export.py:210-216`,
  `classification/tests/utils/test_json_utils.py:29,47`,
  `classification/tests/views/test_classification_view.py:183,267`, `genes/tests/test_hgvs.py:69`.
- **`try / self.fail() / except` instead of `assertRaises`** —
  `upload/tests/vcf/test_vcf_detect_build.py:27,36`, `library/tests/test_utils.py:114`,
  `snpdb/tests/test_data_archive_mixin.py:73` (this one inspects the exception message, so use
  `assertRaises` as a context manager and check `cm.exception`), and
  `classification/tests/models/test_classification_quirks.py:76` (`raise
  self.failureException(...)` — same anti-pattern).
- **`@skip`s** — all seven (count verified by sweep; no `skipIf`/`skipUnless` anywhere), and
  what to do with each:
  - `classification/tests/views/test_classification_view.py:43,188,271,296` — bare `@skip`s,
    added in commit `4b78fa501` ("SACGF/variantgrid_private#3740 - skip failing unit tests").
    No reason in the code (the "awaiting a fix from James" story isn't recorded anywhere in the
    repo). Leave the tests in place; give each `@skip` a reason string citing the private issue
    (`@skip("variantgrid_private#3740")`) so the hold is legible without asking around.
  - `upload/tests/vcf/test_vcf_detect_build.py:40` — open issue #1857, already annotated.
    Legitimate hold, no action.
  - `genes/tests/test_hgvs.py:101` — bare `@skip`; the reason ("needs Ensembl contigs") is in
    the docstring. Mockable, worth un-skipping — or at least move the reason into the decorator.
  - `sync/tests/test_historical_converter.py:38` — bare `@skip` with `# No longer required` as
    the first body line; the only one marked as genuinely dead. Delete the test; the remaining
    `test_latest_keys` keeps the file alive.
- **Stale `TODO` markers inside tests that say the test doesn't work** —
  `analysis/tests/test_clone_nodes.py:139` ("This isn't actually working properly") and `:154`
  ("This doesn't run the test, as it calls `_get_node_q` which MergeNode doesn't implement").
  `test_clone_merge_node` currently asserts nothing. Fix or delete.

---

## Explicitly reviewed and kept (don't re-audit these)

So the next pass doesn't burn effort on them:

- **`*/tests/test_urls.py` smoke tests** across all apps. Zero direct assertions, but they
  delegate to `URLTestCase._test_urls`, which asserts status codes. Cheap protection against
  view-level 500s, and CLAUDE.md documents the pattern.
- **Query-count tests** (`snpdb/tests/test_query_counts.py`, `ontology/tests/test_query_counts.py`,
  `library/tests/test_jqgrid_paginator.py`, `annotation/tests/test_ptc.py`). `assertNumQueries`
  is brittle-looking but these guard specific N+1 regressions and say so.
- **The `str(Q)` substring tests** in `analysis/tests/test_backfill_flag_branches.py` and
  `test_damage_node_v4/v5.py`. Asserting on a `Q` repr is fragile, but they cover real
  version/backfill branching that has no other observable output at that layer. One real defect:
  `test_unbackfilled_required_allow_null_includes_per_field_isnull` (line 58) is named and
  documented for the `allow_null=True` case but passes `splice_allow_null=False` and asserts the
  isnull clauses are *absent* — rename it.
- **`library/tests/test_decorator_audit.py`, `test_signal_receiver_registration.py`** — AST
  audits of the whole codebase. Unusual shape, but each guards a real silent-failure mode
  (Python 3.13 `@classmethod @property` removal; an unused-import autofix unregistering signal
  handlers), and both include a self-check so they can't pass vacuously.
- **The MME suite** (20 files, ~2,200 lines) — consistently good; every test maps to a spec
  requirement or an eligibility/contact-resolution branch.
- **`classification/tests/views/test_classification_view.py`** — 333 lines, all four tests
  `@skip`ped, but this is a deliberate hold (skipped in `4b78fa501` against
  variantgrid_private#3740) rather than dead code.
  Not a deletion candidate. The only ask is a reason string on each `@skip` (see hygiene sweep)
  and dropping the two `print(json.dumps(...))` calls at lines 183 and 267 so they don't fire
  when the tests are re-enabled.
- **`annotation/tests/test_annotation_disk_space.py`, `test_data_archive.py`,
  `test_gene_annotation_command.py`, `test_link_gene_annotation_release.py`,
  `snpdb/tests/test_variant_zygosity_preview_extra.py`, `test_variants_to_vcf.py`** — model
  examples of the rule: each test names a branch and fails for one reason.

---

## Suggested order

1. Tier 1 (deletes) — least judgement required.
2. Tier 2 — decide fix-vs-delete per test. (The gnomAD `r2.1` crash is unreachable on
   migrated deployments — see 2.1 — so it's a cleanup decision, not an urgent fix.)
3. Hygiene sweep — mechanical, mostly `sed`-able, makes the remaining `BUG:` markers meaningful.
4. Tier 3 — rename/rewrite; no behaviour change, so it can land in one commit per file.
5. Tier 4 trims.
6. Tier 5 — the quad/trio mixin is the only item here with real design work in it.
