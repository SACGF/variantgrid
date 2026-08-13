# Issue #1678 — Remove PyHGVS

Agreed at triage 13/08/2026: pyhgvs is no longer called during import or search, the only live caller is the
admin HGVS resolution tool (as a second opinion against biocommons), and everyone is happy for it to go.
Existing `HGVSConverterVersion` rows naming `PYHGVS` stay in the database as plain text — nothing maps that
text back to code.

`hgvs_shim` stays as a dependency. It carries its own `PyHGVSConverter`, but imports `pyhgvs` inside a
`try/except ImportError` and sets `_PYHGVS_AVAILABLE = False` when absent, so the `pyhgvs` package can be
dropped from our requirements while `hgvs_shim` remains. Moving off the shim entirely is a separate,
larger job.

---

## What pyhgvs is in this codebase today

| Layer | Location |
| --- | --- |
| Converter adapter | `genes/hgvs/pyhgvs/hgvs_converter_pyhgvs.py` (122 lines) + empty `__init__.py` |
| Combo checker (only ever biocommons vs pyhgvs) | `genes/hgvs/hgvs_converter_combo.py` (81 lines) |
| Enum members | `HGVSConverterType.PYHGVS`, `HGVSConverterType.COMBO` in `genes/hgvs/hgvs_converter.py:11-19` |
| Factory branches | `genes/hgvs/hgvs_matcher.py:117-151` |
| Transcript data feed | `TranscriptVersion.pyhgvs_data` + `cdot.pyhgvs.pyhgvs_transcript.PyHGVSTranscriptFactory` in `genes/models.py:16,1046-1054` |
| Admin tool | `classification/views/views_hgvs_resolution_tool.py:133-143` |
| One-off command | `classification/management/commands/fix_legacy_imported_allele_info_hgvs_error_provided_reference_length.py:7` |
| Deployment check | `variantgrid/deployment_validation/library_version_checks.py:28` |
| Tests | `genes/tests/test_hgvs.py` (7 pyhgvs test methods + 2 direct pyhgvs uses) |
| Packaging | `requirements.in:43`, `requirements.txt:413`, `environment.yml:98` |
| Prose | `variantgrid/static_files/default_static/page_help/genes/unmatched_hgvs_help.html:3`, `CLAUDE.md:68`, `claude/research/genes.md:127`, assorted comments |

Migrations (`genes/migrations/00{37,48,57,67,69}*`, `classification/migrations/0144`, `0147`) mention pyhgvs
and stay exactly as they are — the genes ones are inert `ManualOperation` stubs whose only live content is a
comment, and the classification ones only touch rows that a fresh database does not have. Two of them import
live code and are worth recognising rather than reacting to: `0144` reads
`HGVSConverterType[settings.HGVS_DEFAULT_METHOD.upper()]` and calls `metadata.version('pyhgvs')` only in the
non-biocommons branch, all behind an `if rvi_qs.exists()` guard; `0147` iterates `for converter_type in
HGVSConverterType`, so a shorter enum simply means fewer no-op `.update()` passes over an empty table.

Nothing in `../variantgrid_sapath` references pyhgvs.

---

## Scope

**In:** everything in the table above except the migrations row.

**Deferred:** `hgvs_shim` itself, and the biocommons-only simplifications that fall out of this
(§7) can be taken or left with the same commit.

---

## §1 — Delete the converter

- Delete `genes/hgvs/pyhgvs/` (both files).
- Delete `genes/hgvs/hgvs_converter_combo.py`. Its `_call_converters` compares results from a converter list
  that only ever held `[BioCommonsHGVSConverter, PyHGVSConverter]`; with one converter left it is an identity
  wrapper.

## §2 — `genes/hgvs/hgvs_converter.py`

- `HGVSConverterType` keeps `BIOCOMMONS_HGVS = 2` and `CLINGEN_ALLELE_REGISTRY = 4` at their current values, so
  any existing `.name` text and any operator habit around the enum stays stable.
- `is_internal_type` becomes `self is HGVSConverterType.BIOCOMMONS_HGVS`. It is read at
  `hgvs_matcher.py:432` and `:570` to decide local-vs-ClinGen, and that distinction survives.
- The abstract `HGVSConverter` class (lines 76-148) loses its last VG subclasses here — `PyHGVSConverter` and
  `ComboCheckerHGVSConverter` were the only two. `BioCommonsHGVSConverter` and `ClinGenHGVSConverter` both
  descend from `hgvs_shim.BioCommonsHGVSConverter` instead. Delete the class and have
  `HGVSConverterFactory.factory`'s return annotation use `hgvs_shim.HGVSConverter`, which is what it has
  actually been returning. Keep `HGVSConverterType`, `HgvsMatchRefAllele` and `HgvsOriginallyNormalized` —
  all three are used widely.

  Check as you go: `_hgvs_string_validation` and `description()` exist on both the VG ABC and the shim base.
  VG's `BioCommonsHGVSConverter` defines its own `description(describe_fallback=True)` at
  `genes/hgvs/biocommons_hgvs/hgvs_converter_biocommons.py:143`, and callers at `hgvs_matcher.py:264,433,484,628`
  pass that kwarg, so that path is already independent of the ABC.

## §3 — `genes/hgvs/hgvs_matcher.py`

- Drop the `ComboCheckerHGVSConverter` and `PyHGVSConverter` imports (lines 29-30).
- In `HGVSConverterFactory.factory`, keep the `BIOCOMMONS_HGVS` and `CLINGEN_ALLELE_REGISTRY` branches and the
  trailing `HGVSImplementationException`. The commented-out `settings.DEBUG → COMBO` block at lines 126-128
  goes with it.
- Update the class docstring at line 119 (the #839 pyhgvs-vs-biocommons decision) and the incidental mentions
  at lines 501 and 587 to name the library actually in play.
- `_gene_symbol_data_provider` (lines 696-704) exists to build a `DjangoTranscriptDataProvider` for converters
  that have no `hdp` — i.e. only PyHGVS. With biocommons and ClinGen both carrying `self.hdp`, this reduces to
  `return self.hgvs_converter.hdp`.

## §4 — `genes/models.py`

- Delete the `PyHGVSTranscriptFactory` import (line 16) and the `pyhgvs_data` cached property (lines 1046-1054);
  its three call sites are all inside the converter deleted in §1. `cdot` itself stays — only its optional
  `cdot.pyhgvs` submodule imports `pyhgvs`, and nothing else touches it.
- Line 888's "We've modified PyHGVS to be able to handle this" comment sits on unrelated code; read it in
  context and reword or drop it.

## §5 — Admin HGVS resolution tool

`classification/views/views_hgvs_resolution_tool.py:133-143`: remove the two `HGVSConverterType.PYHGVS`
matchers. Both branches keep working:

- `display_clingen_separately` on → `[BIOCOMMONS_HGVS (no clingen), CLINGEN_ALLELE_REGISTRY]`, still a
  two-way comparison, which is the case the `all_equal_*` columns were built for.
- off → `[BIOCOMMONS_HGVS]` with ClinGen fallback, plus any "Previously resolved" `ImportedAlleleInfo` rows,
  which is still a useful comparison against what is stored.

`classification/templates/classification/hgvs_resolution_tool.html` needs no change — it renders
`output.matcher_name` from whatever list it is given. Retitle the link label at
`classification/templates/classification/imported_allele_info_detail.html:59` (currently
`"Test pyhgvs/biocommons"`).

## §6 — Remaining single references

- `classification/management/commands/fix_legacy_imported_allele_info_hgvs_error_provided_reference_length.py`:
  drop `from pyhgvs import InvalidHGVSName` and whatever `except` clause uses it — the same file already
  catches `hgvs.exceptions.HGVSInvalidVariantError`. Read the handler before editing; this command reprocesses
  legacy rows and now runs under biocommons, which raises the biocommons exception.
- `variantgrid/deployment_validation/library_version_checks.py`: remove the `"pyhgvs": (0, 12, 4)` entry.
- `variantgrid/settings/components/default_settings.py:427`: the `HGVS_DEFAULT_METHOD` comment lists
  `"pyhgvs"`, `"biocommons_hgvs"`, `"combo"` — leave only `"biocommons_hgvs"` and
  `"clingen_allele_registry"`. All three env files that set it (`shariantcommon.py:153`, `vgaws.py:80`,
  `vgtest2.py:28`) are already `"biocommons_hgvs"`.
- `library/genomics/fasta_wrapper.py`: keep the module — `FastaFileWrapper` is live at
  `snpdb/models/models_genome.py:514`, `snpdb/models/models_variant.py:362,420` and
  `library/genomics/calculate_cancer_mutation_signatures.py:265`. Reword the two docstrings that explain the
  shape in terms of pyhgvs/pygr to describe the pysam contig-slice interface it provides.
- `annotation/external_search_terms.py:33` and `genes/transcripts_utils.py:1-2,72`: comments only. The
  `transcripts_utils` ones are attribution for copied code and are worth keeping as-is; the
  `external_search_terms` one explains why `HGVSName.format_protein()` is not used and can go.
- `variantgrid/static_files/default_static/page_help/genes/unmatched_hgvs_help.html:3`: rewrite the sentence
  to name biocommons hgvs and link https://github.com/biocommons/hgvs. (`variantgrid/sitestatic/` is
  gitignored and regenerated by collectstatic.)
- `CLAUDE.md:68` and `claude/research/genes.md:127`.

## §7 — Packaging

- `requirements.in:43` — remove the `git+https://github.com/SACGF/hgvs#egg=pyhgvs` line, then regenerate with
  `uv pip compile requirements.in -o requirements.txt` (the command in the generated file's header). `pyhgvs`
  is a direct requirement — its `requirements.txt` entry has no `# via` provenance — and the only thing pinned
  through it is `pip==26.1.2`.
- `environment.yml:98` — same line.
- `pysam` stays: `fasta_wrapper` and biocommons both use it directly.

## §8 — Tests

`genes/tests/test_hgvs.py` currently passes (verified 2026-08-13, `--keepdb`).

- Delete the paired pyhgvs test methods, keeping their biocommons twins: `test_hgvs_pyhgvs` (163),
  `test_pyhgvs_reference_diff` (226), `test_pyhgvs_gene_symbol_hgvs` (249),
  `test_pyhgvs_mitochondria_hgvs` (272), `test_pyhgvs_mitochondria_hgvs_on_nuclear_contig_rejected` (287),
  `test_pyhgvs_invalid_trailing_int` (301). The `_test_*` helpers they share with the biocommons tests stay
  as-is, still taking `hgvs_converter_type`.
- `test_hgvs_biocommons` (166) can absorb its `extra_hgvs` inv case inline now that there is no second
  implementation to hold it back — optional tidy-up.
- `test_clean_hgvs` (28-70) builds a `PYHGVS` matcher on GRCh38 purely to reach `clean_hgvs()` and
  `create_hgvs_variant()`. Switch it to `HGVSConverterType.BIOCOMMONS_HGVS`. Biocommons'
  `create_hgvs_variant` is a pure parse (`hgvs_shim/hgvs_converter_biocommons.py:172`) so it needs no
  transcript fixtures, but the converter's `__init__` builds an `AssemblyMapper`/`Babelfish` against the
  Django data provider — run the test to confirm GRCh38 is satisfied by the existing `setUpTestData`, and
  add the fake annotation version for GRCh38 alongside the GRCh37 one if it is not.
- `test_format_hgvs_remove_long_ref` (84-98) constructs `PyHGVSVariant(HGVSName(hgvs_string))` directly.
  Rebuild it on `BioCommonsHGVSVariant(hgvs.parser.Parser().parse(hgvs_string))` — same
  `format(max_ref_length=...)` contract from the shared `HGVSVariant` base. Confirm the expected outputs still
  hold; biocommons formatting of `delins` and long refs is the thing under test, so a differing expectation is
  a real finding, not a test to loosen.
- Module docstring at line 1 ("PyHGVS has its own testing, this is specific to our code") — reword.
- `from pyhgvs import HGVSName` (line 7) and the `genes.hgvs.pyhgvs` import (line 16) go.

## §9 — Verification

1. `grep -ri pyhgvs` over the tree returns only migrations, `.venv`, and `genes/transcripts_utils.py`
   attribution comments.
2. `python3 manage.py test --keepdb genes.tests.test_hgvs classification.tests` — green.
3. `python3 manage.py test --keepdb` — green, with attention to `library.tests.test_decorator_audit`
   (touches `ClinGenHGVSConverter`) and any URL tests covering `hgvs_resolution_tool`.
4. `pip uninstall pyhgvs` in the venv, then re-run 2 and start the server — this is the real check that no
   import path still reaches it and that `hgvs_shim` degrades as advertised.
5. Load `/classification/hgvs_resolution_tool` as a superuser, both with and without
   *display clingen separately*, and confirm both render.
6. Server status settings page (`variantopedia/views.py:322-335`) renders — it constructs an `HGVSMatcher`
   from `settings.HGVS_DEFAULT_METHOD`.
7. `./scripts/linting/run_pylint.sh` for unused imports left behind.

## §10 — Deployment

No migration and no `ManualOperation`: nothing in the database is read back as code. Deployments pick up the
requirements change through the normal upgrade path, and `check_library_versions` stops reporting a pyhgvs row
on the server status page.

The one behaviour change an operator can see is that `HGVS_DEFAULT_METHOD = "pyhgvs"` or `"combo"` in a local
`env_developers/<hostname>.py` now raises `HGVSImplementationException` from the factory at startup of the
first conversion. Worth a line in the release notes; a grep of the env files in this repo finds no such
setting.
