# TSO 500 Phase 2 — per-file loader refinements

Phase 2 of [`tso500_overall_plan.md`](tso500_overall_plan.md): make each of the five files in
`upload/test_data/tso500/` import with the right values in it. Five items, no model dependency on
Phase 1, each testable against one test file.

The stage-by-stage findings this builds on are in [`tso500_ingestion_plan.md`](tso500_ingestion_plan.md);
the files and their gotchas are in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).

Phase 0 landed in [#1709](https://github.com/SACGF/variantgrid/pull/1709), so `_DragenExonCNV.vcf`
reaches the importer and `AllFusions.csv` can be uploaded at all. Neither loads *correctly* yet — that
is this phase.

---

## Where the work lands

The VCF import runs in this order (`upload/import_task_factories/abstract_vcf_import_task_factory.py:137`):

```
  Create Data from VCF Header   ImportCreateVCFModelForGenotypeVCFTask
                                  create_vcf_from_vcf        -> genome build from header
                                  configure_vcf_from_header  -> VCFInfo/Filter/Format, sample fields
                                  handle_vcf_source          -> VCFSourceSettings by ##source
  Preprocess VCF                preprocess_vcf
                                  write_cleaned_vcf_header   -> the header every split file gets
                                  | cat file
                                  | vcf_clean_and_filter     -> contigs, positions, FILTER, skips
                                  | bcftools norm            -> split multi-alt, --check-ref=s
                                  | bcftools view --no-header
                                  | vcf_clean_alts           -> symbolic ALT validation
                                  | split
  Data insertion                BulkGenotypeVCFProcessor.process_entry (per split file)
```

The VCF model — `source`, the sample field columns, anything `VCFSourceSettings` sets — exists before
preprocess runs, so source-driven settings can steer the pipe stages as well as the insert.

## Findings that change the shape of the work

Four things turned up on the way in that the overall plan's one-line descriptions don't account for.

**Undeclared FILTER values are not merely warned about — they kill `bcftools norm`.** Verified against
bcftools 1.20:

```
[W::vcf_parse_filter] FILTER 'LowUniqueAlignments' is not defined in the header
[E::vcf_format] Invalid BCF, the FILTER tag id=4 at chr1:100 not present in the header
```

`bcftools view` passes the record through; `norm`, which is the stage we actually run, does not. So
"let undeclared values through the header stage" is not available as written — the header has to
declare them. Adding a `##FILTER` line for the value makes norm accept it (also verified). That turns
item 3 from a deletion into a design question about where the declaration comes from.

**`cnv.vcf` has no `##source` line, so `VCFSourceSettings` cannot reach it today.** It carries
`##DRAGENVersion` and `##DRAGENCommandLine` instead; `cyvcf2_header_get(reader, "source")` returns
empty and `handle_vcf_source` does nothing. `hard-filtered.vcf` is the same. Only
`SpliceVariants.vcf` (`##source=SpliceGirl 1.0.0.614`) and `_DragenExonCNV.vcf`
(`##source=LrCalculator 1.0.0.11`) are reachable by header alone. Item 1's upload metadata changes
this — see the note under item 4 for why item 4 still doesn't lean on it.

**In the splice VCF, `FORMAT/AD` is already exactly `INFO/ALTDEDUP` and `FORMAT/DP` is exactly
`INFO/REFDEDUP`** — checked, all 17 records. So remapping the two FORMAT fields produces the
`ALTDEDUP`/`REFDEDUP`-derived VAF the plan asks for, with no INFO parsing. Note the header's
descriptions are transposed against their own IDs (`AD` is labelled "SpliceSupportingReads" but its
prose says "do not support deletion", and vice versa for `DP`); the data settles it — `AD` is
alt-supporting, `DP` is reference.

**`SEGID` is already preserved and already queryable.** `CohortGenotype.info` is a real `JSONField`
(`snpdb/models/models_cohort.py:751`), so nothing is lost at import and `info__SEGID=…` works today.
See item 5.

---

## 1. Upload metadata — genome build and source

Two independent facts about a file that its header does not reliably carry. `_DragenExonCNV.vcf` has
no `##contig` lines and a `##reference` of `file:resource_bundle/hg19_decoy/` that
`GenomeBuild.get_name_or_alias` cannot resolve, so the build stays null,
`ImportCreateVCFModelForGenotypeVCFTask` sets `REQUIRES_USER_INPUT` and the pipeline terminates early
(`upload/tasks/vcf/genotype_vcf_tasks.py:63`). `cnv.vcf` and `hard-filtered.vcf` declare no `##source`,
so nothing keyed on the caller can reach them.

**Not by pattern-matching the reference path.** There is both a build *named* `hg19` and `hg19` as
GRCh37's *alias* (`snpdb/migrations/0002_initial_data.py:1425`), so that string is ambiguous by
construction — and this file is GRCh37 (`chrM` 16569), not hg19. The submitter knows the answer;
guessing it from a vendor path is how you import coordinates against the wrong build.

**Not via seqauto either.** The only `genome_build` in the whole seqauto app is a derived property on
`BamFile` (`seqauto/models/models_seqauto.py:725`) that queries linked VCFs and, failing that, logs a
warning and returns `GenomeBuild.legacy_build()`. Backend VCFs do not know their build; the one place
seqauto claims to is a documented guess. Beyond that, the build is a fact the *submitting client*
holds per file, and a pre-registered seqauto object would have to be told it by that same client — an
extra hop that only exists on deployments running seqauto.

**Design: a `FileUpload.metadata` JSONField (default `{}`), supplied by whatever submits the file.**

A blob rather than typed columns because `FileUpload` is shared across every upload type, and the keys
are per-type: a VCF takes `genome_build` and `source`, a BED takes `genome_build`, a gene list takes
neither. `GenomeBuild`'s primary key *is* its name and it is an `ObjectManagerCachingImmutable`, so a
stored name buys nearly the same integrity as an FK — which is what tipped this away from
`FileUpload.genome_build` as a real column. Once there is a blob, one field living outside it is the
inconsistency that makes the next person ask where their key goes.

Surfaces that set it:

- `handle_file_upload(user, django_uploaded_file, path=None, metadata=None)` (`upload/views/views.py:144`).
- `APIFileUploadView` (`upload/views/views_rest.py:42`), beside the existing `path`/`force`. This is
  the route the TSO 500 client uses.
- `manage.py import_vcf`, via typed `--genome-build`/`--source` flags that write into the blob — a
  command line wants named flags, not a JSON string. This is also how you verify this locally.

**Validate at upload time, not at import time.** This is the condition that makes the blob as safe as
a column, and it is easy to skip. `handle_file_upload` already calls `get_uploaded_file_type` before
`process_uploaded_file`, so the file type is known while the client is still connected: reject an
unknown key, or a `genome_build` that will not resolve, with a 400 there. Without it a typo
(`genome_build` vs `genomeBuild`) is silently ignored, the import falls back to header detection, and
the failure surfaces as `REQUIRES_USER_INPUT` with nothing pointing at the cause. Simplest form is a
flat allow-list of known keys; the more correct form has each import task factory declare the keys it
accepts.

How each key is consumed:

- **`genome_build`** — `create_vcf_from_vcf` (`upload/vcf/vcf_import.py:187`) still runs header
  detection first. Where detection finds a build and it **disagrees** with the declared one, fail the
  import rather than picking a winner: disagreeing contig lengths mean the coordinates are not what
  the submitter thinks they are. Where detection finds nothing, the declared build is used.
- **`source`** — populates `vcf.source` in `configure_vcf_from_header`, so `handle_vcf_source` matches
  `VCFSourceSettings.source_regex` against it by exactly the same path as a header-supplied value. A
  client-supplied value wins over a present `##source`; the header text remains in `vcf.header` either
  way.

**Coupling worth naming:** once clients supply `source`, their string becomes part of VG's
configuration contract, because `source_regex` has to match something a client invented. That wants a
stable agreed vocabulary — `DRAGEN TSO500 CNV` and the like — not free text.

Worth adding alongside, cheaply: the "server has exactly one build with annotation" fallback that
`ImportBedFileTask._get_genome_build` (`upload/tasks/import_bedfile_task.py:63`) already applies to
BEDs and VCFs do not have. It doesn't help a two-build deployment so it isn't a substitute.

**Tests.** `upload/tests/vcf/test_vcf_detect_build.py` for the detection-disagrees case;
`upload/tests/vcf/test_vcf_processors.py`'s `_create_fake_upload_step_and_vcf` harness already builds
a `FileUpload` and calls `create_vcf_from_vcf`, so the pickup path tests there. `upload/tests/test_api.py`
for the 400s.

## 2. `VCFSourceSettings` field mapping for SpliceGirl

`configure_vcf_from_header` binds `AD` and `DP` by name (`upload/vcf/vcf_import.py:237`). In this file
`AD` is `Number=1` splice-supporting reads and `DP` is *reference* reads, so:

- `get_ref_alt_allele_depth_function` takes the `DEFAULT_AD_FIELD` branch → `variant.gt_ref_depths` /
  `gt_alt_depths`, which return `-1` because `AD` is not `Number=R`. Allele depth lands empty.
- VAF, computed in `finished_locus` from `locus_ad_sum`, is therefore missing.
- The read depth column shows the reference read count.

The mapping mechanism already exists — `VCF.ref_depth_field`/`alt_depth_field` are columns and
`get_ref_alt_allele_depth_function` already has the explicit branch (it is the FreeBayes AO/RO path).
What is missing is any way to *set* them: FreeBayes and CLCAD2 are hardcoded by name in
`set_allele_depth_format_fields`.

**Where the mapping lives — three candidates, two rejected.**

- *Per-upload metadata* (item 1's blob) — no. Every SpliceGirl VCF has the same mapping. Putting it in
  the upload payload means every client restates it and clients can disagree with each other. The
  clean split is: the client says *what produced the file*, VG says *what that means*.
- *A third hardcode beside FreeBayes/CLCAD2* — cheapest (a few lines, no migration), but it extends a
  pattern that predates `VCFSourceSettings` and that we would rather migrate off than grow.
- *`VCFSourceSettings`* — yes. It is already the "behaviour keyed on source" table. Don't migrate the
  two existing hardcodes as part of this; just stop adding to them.

So: `VCFSourceSettings` grows sample-field overrides, applied in `handle_vcf_source`
(`upload/vcf/vcf_import.py:281`) — which already runs last in `configure_vcf_from_header`, so
overrides land on top of the by-name defaults. For `^SpliceGirl`:

| VCF field | value | why |
|---|---|---|
| `allele_depth_field` | *(null)* | `AD` is not the packed ref,alt array cyvcf2 expects |
| `alt_depth_field` | `AD` | splice-supporting reads |
| `ref_depth_field` | `DP` | reference reads |
| `read_depth_field` | *(null)* | there is no total-depth field; `DP` here is not one |
| `allele_frequency_field` | *(null)* | so VAF is derived rather than read |

With ref and alt set explicitly, `get_ref_alt_allele_depth_function` takes its `vcf.ref_depth_field and
vcf.alt_depth_field` branch (`upload/vcf/bulk_genotype_vcf_processor.py:221`) and `finished_locus`
computes `AD / (AD + DP)` — which is `ALTDEDUP / (ALTDEDUP + REFDEDUP)`, per the finding above.

A data migration creates the row, following `snpdb/migrations/0025_default_vcf_source_settings.py`.

**Model shape — decide first.** Null on an override column cannot distinguish "no override" from
"clear this field", and three of the five values above are clears. Recommendation: a single
`sample_field_overrides` JSONField (default `{}`) keyed on the VCF's own sample-field column names —
one column, "set to null" is expressible, and Phase 6 can add `SM` without another migration. The
alternative is six nullable columns plus an `override_sample_fields` boolean, which is more visible in
the admin and more in keeping with the rest of the model.

**Wrinkle — duplicate positions.** `chr2:47637511` appears twice with different `END`. `process_entry`
keys the locus on `(CHROM, POS, ref)`, so both records join one locus and `locus_ad_sum` totals all
four depths: VAF comes out `1/182` rather than `1/91`. This is inherent to how decomposed multi-alt AD
summing works, and including `svlen` in the locus key would change that shared behaviour for real
multi-allelics. Decide whether to accept it or scope a symbolic-only key.

Leave `sample_variants_type` alone on this row for now: `SOMATIC_ONLY` schedules a Calculate
Mutational Signature step per sample (`upload/tasks/vcf/genotype_vcf_tasks.py:89`), which is not what
you want for 17 splice calls.

**Test.** `upload/tests/vcf/test_vcf_processors.py` runs `create_vcf_from_vcf` →
`configure_vcf_from_header` → `handle_vcf_source` → `process_entry`, which is the whole path.

## 3. Preserve `LowUniqueAlignments`

The splice header declares `LowQ` only, so `vcf_clean_and_filter` strips the undeclared code from 15
of 17 records and reports `Skipped 15 'LowUniqueAlignments' FILTER`
(`upload/management/commands/vcf_clean_and_filter.py:133`). That is the filter separating the two
`PASS` oncology calls from the background — silently deleted, on any caller, not just this one.

Given the bcftools finding above, the value cannot stay in the FILTER column past `bcftools norm`.
Three places it could go — the third is what was built.

**(a) Scan for them when writing the cleaned header** and emit `##FILTER` lines. Costs a full extra
read of the file on *every* import, whether or not it needs one — and on a `.gz` that means
decompressing the whole thing.

**(b) Source-driven declarations** on the same `VCFSourceSettings` row as item 2. Free, but only fixes
callers someone has configured.

**(c) Move the value into an INFO field the pipeline owns, and put it back at insert.** ✔ built.
`vcf_clean_and_filter` already walks every line and already separates declared from undeclared filters
(`upload/management/commands/vcf_clean_and_filter.py:144`) — so the discovery costs nothing new. Instead
of dropping the undeclared values it appends them to `INFO/VG_UNDECLARED_FILTERS`, declared in the
cleaned header beside the `BCFTOOLS_OLD_VARIANT` and `VG_LIFTOVER_SWAPPED_REF_ALT` lines
`_write_cleaned_header` already injects — a fixed ID we own, so nothing has to be discovered ahead of
time. `BulkGenotypeVCFProcessor._restore_undeclared_filters` reads it back and hands it to
`convert_filters`, whose `_add_undeclared_filter_code` already creates the `VCFFilter` rows.

Why (c) wins: no extra pass, no retry, and no cap on header lines — the ceiling is
`_add_undeclared_filter_code`'s existing code assignment, which already handles running out. Verified
against bcftools 1.20 that `norm` copies the INFO onto **both** records of a split multi-allelic and
preserves it through `--check-ref=s` rewriting `REF=N`, so the many→one mapping from source row to
inserted variants comes free rather than needing a variant-keyed sidecar.

Two details it needs: INFO is the last column in a VCF with no samples, so it carries the line's
newline and the value has to be spliced in before it; and the separator is `|`, since the VCF spec
already bars whitespace and `;` from FILTER IDs while `,` would read as a multi-value INFO. The field
is popped in `_get_info_json` so our own plumbing doesn't land in `CohortGenotype.info`.

`create_vcf_filters` reads the *original* header at the header stage, so no `VCFFilter` row is created
for the value there; `_add_undeclared_filter_code` makes one at insert time with description
"Undeclared filter". That path already has a test — `test_genotype_processor_undeclared_filter`.

**Tests.** `TestUndeclaredFiltersMovedToInfo` in `upload/tests/test_vcf_clean_and_filter.py` for the
move, and `test_undeclared_filters_restored_from_info` in `upload/tests/vcf/test_vcf_processors.py` for
the restore.

## 4. Skip copy-neutral `cnv.vcf` rows

9 of 25 records (406 of 500 in a real run) are `ALT=.` with `END`, `REFLEN` and `SEGID`. They import as
reference variants at a single base, with the segment span discarded, and they carry no call.

**Not a blanket `ALT=.` skip.** Reference variants are first class — `Variant.REFERENCE_ALT`
(`snpdb/models/models_variant.py:582`) — and `vcf_clean_alts` allows `.` deliberately ("Can '.' for
reference").

**Rule: skip a record whose ALT is `.` *and* whose INFO carries `END`.** A reference call over a span
is a no-call region, not a variant. That catches `cnv.vcf`'s copy-neutral segments (`END`, no
`SVTYPE`) without touching a single-position reference record, and it is the same posture the importer
already takes to gVCF reference blocks — "we only store variants",
`upload/migrations/0011_one_off_convert_skipped_gvcf_to_standard_vcf_info.py`.

**Why not scope this to a source, now that item 1 makes `cnv.vcf` reachable.** It could be, and that
option only exists because of item 1 — but a general rule is better here. It is not new policy, just
the gVCF reference-block statement applied to the `ALT=.` spelling. Source-scoping means the junk
imports by default and every new CNV caller needs a config row before it stops, and it does nothing
for a human dragging a file into the browser with no metadata attached. If the rule ever drops
something wanted, `VCFSourceSettings` is where the exception goes.

**Where:** `vcf_clean_and_filter`, beside the existing `skip_patterns` machinery, counted into
`--skipped-records-stats-file` so it surfaces as `Skipped 9 '…' records` on the pipeline page like
every other skip. The command does not parse INFO today; `parse_vcf_info_column`
(`library/genomics/vcf_utils.py`) is there if a real parse is wanted — `vcf_clean_alts` already uses it.

**Decide first:** this rule also drops `_DragenExonCNV.vcf`'s `Undetermined` BRCA2 row
(`SVTYPE=CNV;END=32973805`, `ALT=.`), leaving that file with one record. A per-gene "we looked and
could not tell" is arguably a different statement from a copy-neutral segment, and that record exists
in the test data specifically to exercise the `Undetermined` filter. Settle whether it should survive
before writing the rule.

## 5. Resolve gene symbols through `GeneSymbol` — move this to Phase 6

`SEGID=MYCL1` uses an older symbol than the rest of the pipeline (`MYCL`), but nothing is lost at
import: `CohortGenotype.info` is a real `JSONField`, so the value is stored and `info__SEGID=…` is
queryable today. That is why the overall plan's Phase 2 done-when doesn't mention this item.

The actual gap is that the stored value is *the caller's string*. "CNVs in MYCL" will not find the row
that says `MYCL1` unless something resolves aliases, and JSON cannot be joined to
`GeneSymbol`/`GeneSymbolAlias`. So this is a read-time-versus-import-time choice, not a data-loss one:
resolve on display and it is free but unindexable; resolve at import into a real column and grids can
filter and sort on it.

That is precisely Phase 6's question, and it is the same question `SM` and `CN` already pose there
(#1558's "surfacing `CohortGenotype.info["CN"]` … not queryable"). **So this moves to Phase 6, with
those**, not to Phase 5.

The mechanism it will want already exists: `GeneSymbolMatcher.get_gene_symbol_id_and_alias_id`
(`genes/gene_matching.py:45`) resolves through `GeneSymbolAlias` and returns both the resolved symbol
and the alias that got there — the same shape Phase 5's fusion parser needs for `SEPT14` → `SEPTIN14`.
Worth confirming against a real database that `MYCL1` resolves to `MYCL` (NCBI carries it), since that
is the assumption the item rests on.

---

## Decisions — as settled

All five were settled before implementation and built as below (PR #1712).

| Decision | Item | Settled as |
|---|---|---|
| The metadata key vocabulary | 1 | Per-factory declaration — `ImportTaskFactory.get_metadata_keys()`, empty in the base, `{genome_build, source}` on `GenotypeVCFImportFactory`, the one factory reaching `create_vcf_from_vcf`. The `source` strings are documented in the test-data README but nothing is keyed on them yet, so the contract stays uncommitted until something consumes it |
| Sample-field overrides: one JSONField, or six nullable columns plus a boolean | 2 | One `sample_field_overrides` JSONField, keys validated in `clean()` against `OVERRIDABLE_SAMPLE_FIELDS` — three of the five SpliceGirl values are clears, and "key present ⇒ set it, including to null" only works if a typo can't silently vanish |
| Undeclared FILTER: scan when writing the cleaned header, or declare per source | 3 | Neither. Review found a third option that costs nothing: `vcf_clean_and_filter` already separates declared from undeclared filters on the pass it already makes, so it moves the undeclared ones into `INFO/VG_UNDECLARED_FILTERS` and the genotype processor restores them at insert. No extra file read, no retry — see item 3 |
| Should `ALT=. + END` drop `_DragenExonCNV.vcf`'s `Undetermined` BRCA2 row | 4 | Kept — rule narrowed to `ALT=.` + `END` + **no `SVTYPE`**. Verified: all 9 `cnv.vcf` copy-neutral rows carry no `SVTYPE`, the BRCA2 row carries `SVTYPE=CNV`. `SVTYPE` present is the caller declaring a structural-variant call rather than a segmentation interval, which is the distinction the rule draws. gVCF reference blocks are unaffected — they carry `END` with no `SVTYPE`, but their ALT is `<NON_REF>` |
| Duplicate-position splice records: accept halved VAF, or key symbolic loci on svlen | 2 | Accept the join (1/182). Both affected records are `LowQ` background calls, not the PASS oncology calls, and the locus key is shared with every other VCF's multi-allelic AD summing. Recorded rather than left silent: a `ModifiedImportedVariant` with the new `SHARED_LOCUS` operation and an `operation_detail` column carrying the summed denominator and the record's own ref/alt, so the per-record value stays reconstructable |

Settled during planning, recorded so they aren't relitigated: metadata is passed at upload rather than
read off seqauto objects; it is a JSON blob rather than typed columns on `FileUpload`; the sample-field
mapping lives in `VCFSourceSettings` rather than in the upload payload or a third hardcode; item 5
moves to Phase 6.

## What the build found that this plan didn't

- **`hg19` is ambiguous by more than the plan says.** `GenomeBuild.get_name_or_alias("hg19")` raised
  `MultipleObjectsReturned`, not `DoesNotExist`, so "reject a build that will not resolve" didn't cover
  it. That build is now disabled, so it resolves to GRCh37 here, but `GenomeBuild.enabled` is
  per-deployment DB state — the check stays as a 400. The real fix is clients sending a build's own name
  (`GRCh37`), now documented in the API schema, `import_vcf --genome-build` and the test-data README.
- **The undeclared FILTER doesn't need declaring at all.** The plan's two options both put a `##FILTER`
  line in the header, which means knowing the IDs before any record is read — hence a scan. But the
  value only has to get *past* `bcftools norm`, and `_add_undeclared_filter_code` already handles it at
  insert. Carrying it in an INFO field the pipeline owns does that with no extra pass: `norm` copies
  INFO onto both records of a split multi-allelic and preserves it through `--check-ref=s`, so the
  source-row → inserted-variant mapping comes free.
- **`old_multiallelic` could not carry the shared-locus explanation.**
  `get_other_loci_variants_by_multiallelic` excludes same-locus variants ("split into separate loci"),
  which is the opposite of this case — hence the dedicated `operation_detail` column.
- **`create_vcf_from_vcf` never saved the VCF** despite its "Save header ASAP if case something goes
  wrong" comment; the first save was inside `configure_vcf_from_header`. Added, so a build-mismatch
  failure keeps the header.
- **#1711's checklist said `cnv.vcf` keeps 15 `<DUP>`/`<DEL>` rows — it is 16** (8 + 8), and it
  contradicted its own "16 records" line, since 25 − 9 = 16 makes every survivor symbolic. Issue
  corrected.

## Order of work

Nothing here depends on anything else here, so this is grouped by shared surface rather than by
dependency.

1. **Item 1** — `FileUpload.metadata`, `handle_file_upload`, `APIFileUploadView`, `import_vcf`,
   upload-time validation. First because it is the only item that stops a file loading at all, and
   because items 2 and 4 read better once `source` is reachable on every file.
2. **Items 3 and 4** — `vcf_clean_and_filter` and `write_cleaned_vcf_header`, one test file, no model
   changes.
3. **Item 2** — `VCFSourceSettings` override columns, schema migration, `^SpliceGirl` data migration.
4. **Item 5** — moves to Phase 6.

Commits reference [#1506](https://github.com/SACGF/variantgrid/issues/1506); none of these have an
issue of their own, so raise them individually if they are going to different people.

## Done when

All five test files import with the right values in them:

| File | Expected after Phase 2 | Measured |
|---|---|---|
| `hard-filtered.vcf` | unchanged — 93 records | 93 → 93 ✔ |
| `cnv.vcf` | 16 records; the 9 copy-neutral rows skipped and counted on the pipeline page | 25 → 16, 9 skipped ✔ |
| `SpliceVariants.vcf` | 17 records, VAF = `ALTDEDUP/(ALTDEDUP+REFDEDUP)`, `LowUniqueAlignments` on the 15 background calls and visible on the grid | 17 → 17, 15 `LowQ;LowUniqueAlignments` + 2 `PASS` ✔ |
| `_DragenExonCNV.vcf` | imports against a declared GRCh37, no `REQUIRES_USER_INPUT` | 2 → 2 against a declared build ✔ |
| `AllFusions.csv` | unchanged — Phase 5 | unchanged ✔ |

Splice VAFs: `chr7:55087058` (EGFRvIII) 64/65 = 0.985, `chr7:116411708` (MET exon 14) 91/92 = 0.989,
`chr1:120464432` 1/81 = 0.012, `chr2:47637511` 1/182 = 0.0055 for both records at the shared locus.

Measured by driving `write_cleaned_vcf_header` → `vcf_clean_and_filter` → `bcftools norm` against a
GRCh37 database. Still worth running
`python3 manage.py import_vcf --name <n> --user <u> --genome-build GRCh37 <file>` end to end, which
adds the insert and annotation stages these numbers don't cover.
