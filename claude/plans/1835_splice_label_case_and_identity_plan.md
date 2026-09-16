# #1835 — Splice events: label case, and a naming table used as an identity oracle

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-16
Status: in progress

Found importing non-standard classification targets (EGFRvIII, EGFRvIVa, METExon14skipping) under
[#1835](https://github.com/SACGF/variantgrid/issues/1835). Of 4292 `ImportedAlleleInfo` records 14 failed; two are the
case bug in Part 1, two more are the design problem in Part 2. The other ten are the lab's problem (typo'd fusion
partners, multi-allele `c.[a;b]`, bad HGVS, an unresolvable RYBP transcript) except `NM_198253.2(TERT):c.-146C>T`
failing a 5'UTR bounds check, which recurs and is its own bug.

## Part 1 — the alt's case (a bug, fix on its own)

### What happens

A gene-level alt is a `Sequence` row, and everything that inserts a Sequence upper-cases the allele first:
`library/genomics/vcf_utils.py:vcf_get_ref_alt_svlen_and_modification` for VCF records, and
`snpdb/models/models_variant.py:VariantCoordinate` in `from_symbolic_match` and `from_explicit_no_svlen`. So a splice
event that went through the VCF pipeline is stored as `<SPLICE:HGNC:3236:VIII>`. But `genes/gene_splice.py:ResolvedSpliceEvent`
formats its alt with the label in its seeded case - `<SPLICE:HGNC:3236:vIII>` - and that is the coordinate the
classification import hands to `snpdb/variant_pk_lookup.py:VariantPKLookup`. `add` finds no Sequence for the alt, files it
as unknown, and the linking step in `classification/tasks/classification_import_process_variants_task.py` reports
`Variant '... <SPLICE:HGNC:3236:vIII> ...' not inserted!`. Records 1630 (`EGFRvIII`) and 3918 (`METExon14skipping`) on this box.

Only SPLICE is affected. FUSION / GAIN / LOSS alts are kind, namespace and an integer, none of which has case.

The same mismatch bites the records that did import. `genes/gene_splice.py:get_splice_event_variant` parses the label
back off the stored alt (`VIII`, `EX14SKIP`) and `SpliceEventVariant.splice_event` looks the row up by exact label, so
it finds nothing and the display falls back to the canonical string:

| Variant  | alt label | `splice_event` | displays          | should display            |
|----------|-----------|----------------|-------------------|---------------------------|
| 16985914 | EX14SKIP  | None           | MET EX14SKIP      | MET exon 14 skipping      |
| 16983845 | VIII      | None           | EGFR VIII         | EGFRvIII splice variant   |
| 16983799 | V7        | AR-V7          | AR-V7 splice variant | (right, by luck - the seed is upper) |

The tests missed it because `genes/tests/gene_level_test_utils.py:get_sequence` stores the seq exactly as given, so the
test double of the pipeline never upper-cased anything.

### The fix

The alt is a Sequence and a Sequence is upper-case, so upper is the storage form. Canonicalise at the one boundary
that builds and reads the alt, and let the naming table keep its display case:

1. `library/genomics/vcf_enums.py:GeneLevelSymbolicAlt` - `format` upper-cases the label. `parse` returns the label as
   stored (upper after this change; older rows are already upper). The regex `GENE_LEVEL_LABEL` stays permissive so an
   alt written before this change still parses.
2. `snpdb/models/models_variant.py:VariantCoordinate` - `from_gene_level_match` upper-cases the alt like its siblings do,
   so a gene-level coordinate string is the same coordinate whatever case it arrives in.
3. `genes/gene_splice.py:SpliceEventVariant` - `splice_event` matches `label__iexact`; `canonical_str` and `display` use
   the row's label when there is one, so a variant read back from `<SPLICE:HGNC:7029:EX14SKIP>` prints `MET ex14skip`
   / `MET exon 14 skipping`, never `MET EX14SKIP`. The upper label is the storage form only.
   `ResolvedSpliceEvent.label` keeps the seeded case (its `alt` property goes through `format`), so `Matched MET ex14skip`
   messages and `canonical_str` are unchanged.
4. `genes/models/models_splice_event.py:SpliceEvent.label` gets `db_collation='case_insensitive'`, the collation
   `genes/models/models_gene_level.py:GeneLevelId.symbol_str` already uses, so the `(gene_symbol, label, genome_build)`
   uniqueness is case-insensitive too and a curator cannot add `VIII` beside `vIII`. One schema migration.
5. `genes/tests/gene_level_test_utils.py:get_sequence` upper-cases, so the test pipeline stores what the real one does.
   Every existing splice test then exercises the real storage form.

### Data repair

- The three splice Variants are already stored upper-case. Nothing to rewrite in `snpdb_variant` or `snpdb_sequence`.
- Records 1630 and 3918 hold a `variant_coordinate` string with the lower-case alt and status Failed. A re-match derives
  the coordinate again through `resolve_gene_level` and runs a fresh pipeline:
  `manage.py classification_rematch_stuck --status F --gene-level --older-than-hours 0`. Ship it as a
  `ManualOperation` in a classification migration, gated on a Failed gene-level record existing, the way
  `classification/migrations/0184_one_off_revalidate_gene_level_allele_infos.py` is.
- Re-display check on this box after the code change: `vg inspect variant 16985914` and `16983845`, and the
  classification pages for the records on them, print the row's display. No stored evidence to fix - `splice_label` is
  rendered at report time from `classification/report/case_report_context.py`, not written into the evidence.

### Tests

Keep: one test that `resolve_splice_string("EGFRvIII").alt` equals the alt the pipeline stores, one that a Variant
created with the upper-cased alt reads back as its `SpliceEvent` and displays the row's `display`, and one that a
mixed-case gene-level coordinate string parses to the upper-cased coordinate. `scripts/vg tests --explain` after the
edit for the rest.

## Part 2 — a splice string goes through the same stages as an HGVS

### The problem

`genes/models/models_splice_event.py:SpliceEvent` is documented as a naming table: a junction with no row still imports,
labelled with its own coordinates. The TSO 500 path honours that - `genes/gene_splice.py:SpliceEventResolver.resolve`
takes real breakpoints and the row only decides what the junction is called.

The classification path asks the same table "is this string a real splice event?" through
`genes/gene_splice.py:parse_splice_string` and `resolve_splice_string`, and treats absence as fatal. `EGFRvIVa` and
`EGFRvII` fail not because they cannot be represented but because nobody pre-registered the name. Worse, they fail as
if they were HGVS: `resolved_gene_level` returns None, `update_variant_coordinate` runs the HGVS converter on
`EGFRvIVa` and records `No colon (':') provided`, and `_calculate_validation` tags the record
`transcript_type_not_supported` and `cant_resolve_to_variant_coordinate` because `is_gene_level` is only true once
something resolved.

Decision: there is no name registry. The lab's label is the identity, the way a fusion's identity is its two symbols.
A splice string goes through the stages an HGVS does - tidy, canonicalise, validate, coordinate, match - and mints a
Variant when it passes validation, whether or not anything has been observed under that name. Duplicates are prevented
by canonicalisation, not by a table. `SpliceEvent` keeps one job: the TSO 500 importer turns caller breakpoints into the
label a classification arrives under, so a report's `EGFRvIII` and the caller's junction land on one Variant. It is
never consulted on the classification path.

### Canonical labels

The label is our own format, so it is defined for readability in plain text and for surviving the alt's upper-casing:
lower-case tokens joined by underscores, case carrying no information. The alt stores it upper (`<SPLICE:HGNC:3236:V_III>`)
and `GeneLevelSymbolicAlt.parse` lowers it back, so the round trip is lossless and Part 1's boundary stands.

| shape | written forms accepted | canonical label | displays as |
|---|---|---|---|
| numbered variant | `AR V7`, `AR-V7`, `ARV7`, `AR-V7 splice variant` | `v_7` | `AR-V7` |
| roman variant | `EGFRvIII`, `EGFR vIII`, `EGFR-vIVa`, `EGFRvIVa splice variant` | `v_iii`, `v_iva` | `EGFRvIII`, `EGFRvIVa` |
| exon skipping | `MET exon 14 skipping`, `METex14skip`, `MET Exon14Skipping`, `MET ex14 skipping` | `exon_14_skipping` | `MET exon 14 skipping` |
| coordinates | `AR X_66905968_66914514`, `AR chrX:66905968-66914514` | `grch37_x_66905968_66914514` | `AR GRCh37 X:66905968-66914514` |

The coordinate label gains the build. Today `genes/gene_splice.py:coordinate_label` writes `X_66905968_66914514` and the
Variant sits on the build-independent gene-level contig, so the same numbers in GRCh37 and GRCh38 would collide on one
Variant. The build is part of the identity and goes in the label. A lab's coordinate form is canonicalised under the
record's imported build.

One function owns this: `canonical_splice_label(written) -> Optional[str]`, a token parser per shape, and one
`display_splice_label(gene, label) -> str` formats it back (the row's `display` wins where a `SpliceEvent` row exists,
so the TSO 500 events keep the text the panel uses). `splice_key` and `splice_events_by_key` go; canonicalisation is the
key.

### Stages, mirroring the HGVS path

In `classification/models/classification_variant_info_models.py:ImportedAlleleInfo`:

1. **Recognise.** `looks_gene_level(imported)` - the fusion separator, `genes/gene_copy_number.py:COPY_NUMBER_STRING_PATTERN`
   or `genes/gene_splice.py:SPLICE_STRING_PATTERN`, and not `HGVS_UNCLEANED_PATTERN`. `is_gene_level` and
   `imported_as_c_hgvs` read it, so a gene-level-shaped value never reaches the HGVS converter and is never tagged
   `transcript_type_not_supported`, resolved or not. `SPLICE_STRING_PATTERN` widens to the shapes above (`v[IVX]+[a-z]?`
   for `vIVa`, `chrX:start-end` for the coordinate form).
2. **Canonicalise.** `parse_splice_string` returns `(gene name, canonical label)` from `canonical_splice_label`, with no
   table lookup.
3. **Validate.** `resolve_splice_string` resolves the gene with `allow_unknown=False` (as fusions and copy number do)
   and, for the coordinate form, checks the contig exists in the imported build. Each `resolve_*_string` returns a
   result carrying either the resolved identity or a reason (`gene 'ARHGEF' is not a symbol we know`, `'vX' is not a
   splice label shape we accept`), so the message says which stage refused it.
4. **Coordinate.** As now - `ResolvedSpliceEvent.variant_coordinate` through the VCF insert pipeline.
5. **Match.** As now. A failure at stage 3 sets Failed with the reason as the message and a new validation tag
   `gene_level_unresolved` (severity E) in place of `cant_resolve_to_variant_coordinate`, and the
   `classification/views/imported_allele_info_view.py` grid gains a filter on it. Re-matching after a fix is
   `classification_rematch_stuck --status F --gene-level`, whose `--gene-level` widens to values that look gene-level,
   since a record that failed at stage 3 has no coordinate to match the prefix on.

### SpliceEvent's remaining job

- `label` holds canonical labels: the seed rows change to `v_7`, `v_iii`, `exon_14_skipping` in a data migration.
  `display` stays.
- `SpliceEventResolver.resolve` is unchanged except `coordinate_label` now includes the build.
- `SpliceEventVariant.splice_event`, `display` and `canonical_str` read through `display_splice_label`, with the row's
  `display` preferred. `find_splice_events_for_string` (search) canonicalises the same way and stays lookup-only.
- The docstrings in `genes/models/models_splice_event.py` and `genes/gene_splice.py`, and the splice paragraph in
  `genes/CLAUDE.md`, say the new rule: a splice string that validates mints its Variant; `SpliceEvent` names what the
  caller reports and is not consulted on the classification path. `classification/CLAUDE.md` drops the sentence about
  needing a row.

### Data migration

Existing splice Variants carry the old labels (`V7`, `VIII`, `EX14SKIP`) and coordinate labels without a build.
A `ManualOperation` running a new `splice_labels_canonicalise` command re-points each such Variant's `alt` to the
Sequence for its canonical label (same pk, so classifications, VariantAllele and samples follow), taking the build for a
coordinate label from the VCF the Variant was loaded from. A canonical target that already exists is reported and left
for the user rather than merged - on this box the three events have one Variant each, so none is expected. The
`ImportedAlleleInfo.variant_coordinate` strings are rewritten by the same command. It touches `snpdb_variant`, so it
runs with the user, not from a hook.

### Tests to keep

- `canonical_splice_label` on every written form in the table above gives the canonical label, and on a shape we do
  not accept gives None.
- `resolve_splice_string("EGFRvIVa")` mints an alt with no `SpliceEvent` row; `"NOTAGENE vIII"` returns a reason.
- `ImportedAlleleInfo.get_or_create("EGFRvIVa")` matches; `"ARHGEF::TP53"` ends Failed with `gene_level_unresolved`
  and without `transcript_type_not_supported`.
- `SpliceEventResolver.resolve` on a seeded junction gives `v_iii`; on unknown coordinates gives the build-qualified
  coordinate label; `resolve_splice_string` on the report's `EGFRvIII` reaches the same alt.
- `display_splice_label` round-trips each shape.

### Order

1. Part 1 (case) - in progress, unblocks 1630 and 3918.
2. Canonical labels, stages 1-5, the seed data migration and the canonicalise command. One change; the label format
   change and the string path are not separable, since the seed rows must hold what the string path produces.

### Out of scope

`NM_198253.2(TERT):c.-146C>T` failing `c.-146 coordinate is out of bounds` - a recurrent promoter hotspot rejected by a
5'UTR bounds check, its own investigation.
