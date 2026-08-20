# Backfill VariantAnnotation columns from an annotated VCF (#1675)

## Goal

Correct or fill in individual `VariantAnnotation` columns for an existing `VariantAnnotationVersion` without
re-running the whole VEP pipeline over every variant.

The driving case is #1673: COSMIC keeps renaming the INFO field carrying its sample count (`CNT` in v95/97,
`SAMPLE_COUNT` in v99, `GENOME_SCREEN_SAMPLE_COUNT` in v101), so deployments have `cosmic_count` null across
their annotation partitions. Re-annotating from scratch to recover one integer column costs weeks of VEP.

The work splits into three steps, of which VariantGrid owns the first and third:

1. **dump** — write a VCF of the variants already annotated in a version, carrying `variant_id` in INFO
2. **annotate** — the operator annotates that VCF however suits the column (`bcftools annotate` against a
   COSMIC VCF is the fast path; VEP works too), outside VariantGrid
3. **import** — read the annotated VCF back and update the named columns in batches

## Design

### Its own management command

`manage.py annotation_backfill_columns`, with `--dump` / `--import` modes mirroring the shape of
`annotation_external`. Helpers live in a new `annotation/backfill_columns.py`.

It shares the VCF-writing helper (`annotation.annotation_run_files.write_qs_to_vcf`) and the COPY plumbing
(`upload.vcf.sql_copy_files`) with the rest of the annotation code, and owns everything else. What it does is
the mirror image of `annotation_external` at each step: it dumps variants that **are** annotated rather than
ones awaiting annotation, it leaves `AnnotationRun` / range-lock state alone entirely, it UPDATEs named
columns rather than INSERTing whole rows, and its input comes from a **different** annotation source than the
one that produced the version — which is the whole point, and the opposite of the identity that
`vep_check_annotated_file_version_match` exists to enforce on external imports.

### Columns come from the VEPColumnDef registry

`--columns` names `VariantGridColumn` ids that appear in `annotation.vep_columns.VEP_COLUMNS`. The registry
entry supplies:

- the **formatter** — `cosmic_count` needs `format_pick_highest_int` because COSMIC arrives `&`-joined when
  several records overlap a variant; choice columns (`impact`, `sift`) map to their single-char DB values.
  Reusing the registry formatter is what makes a backfilled value identical to a pipeline-written one.
- the **expected INFO field** (`VEPColumnDef.vep_info_field`), as the default source field
- **which columns exist for this version**, via `applies_to(columns_version=…, vep_version=…)`

The INFO field name is overridable per column, since the source file's naming varies by release:

```
--columns cosmic_count,cosmic_legacy_id
--info cosmic_count=GENOME_SCREEN_SAMPLE_COUNT
```

(The operator can equally rename on the way through with
`bcftools annotate -c INFO/COSMIC_SAMPLE_COUNT:=INFO/GENOME_SCREEN_SAMPLE_COUNT`.)

VEP output is supported by reading the value off the PICK'd CSQ entry instead of a top-level INFO field:
`--csq cosmic_count=COSMIC_SAMPLE_COUNT`. Same registry, same formatter, different place to look.

Every named column is carried through one dump, one COPY and one UPDATE, so the row rewrite — the expensive
part — happens once for the set rather than once per column.

### Updates go through a temp table

Per batch (default 1,000,000 records), inside one transaction:

1. COPY the `(variant_id, col…)` rows into a temp table keyed on `variant_id`
2. `ANALYZE` it
3. `UPDATE <partition> a SET col = t.col FROM tmp t WHERE a.variant_id = t.variant_id AND a.col IS DISTINCT FROM t.col`
4. `TRUNCATE` the temp table

The UPDATE targets the version's partition table by name
(`VariantAnnotationVersion.get_partition_table(base_table_name=…)`), so it touches one version's rows and
never scans the other partitions.

Measured on a 2M-row copy of `annotation_variantannotation_version_2` (782MB, production row width), against
Django's `bulk_update` on the same data:

| Approach | 2M rows | Scaled to 24M |
|---|---|---|
| COPY → temp table → `UPDATE … FROM`, ~2% of rows carry a value | 0.6s + 1.7s = **2.3s** | **~30s** |
| Same, every row changes | 0.6s + 25.4s = **26s** | **~5 min** |
| `bulk_update`, batch_size 1000 | **310s** (6,450 rows/s) | **~62 min** |
| `bulk_update`'s prerequisite fetch (`.only()` 3 cols) | **33s** (61,400 rows/s) | **~6.5 min** |

Three things drive the gap:

- `bulk_update`'s cost is Python, not Postgres. Hand-written `UPDATE … CASE id WHEN … END WHERE id IN (…)`
  batches ran at ~80k rows/s against the same table; `bulk_update` managed 6,450 building `Case(When(...))`
  expressions per batch. It also degrades as batches grow — batch_size 5000 dropped to 4,970 rows/s.
- `bulk_update` keys on `VariantAnnotation.pk`, and a VCF gives us `variant_id`, so it needs a full read pass
  over the partition just to build model instances. `UPDATE … FROM` joins on `variant_id` directly and builds
  no Python objects at all.
- `IS DISTINCT FROM` lets Postgres skip rewriting rows whose value is unchanged. For COSMIC only a small
  fraction of variants carry a count, which is what turns 5 minutes into 30 seconds.

### Writing values vs clearing them

Default: records where the source has no value leave the existing column alone, so a partial annotation file
fills in what it knows. `--clear-missing` makes the file authoritative for the dumped variants, writing NULL
where the source is silent — which is what a full re-annotation of that column wants.

## Implementation

### 1. bgzip the VCF dump

`annotation/annotation_run_files.py:write_qs_to_vcf` currently writes `.gz` output through `gzip.open`, which
tabix cannot index — and `bcftools annotate` wants an indexed target. Write bgzip instead, using the
`BGZipWriter` + `io.TextIOWrapper` pattern already in `snpdb/variants_to_vcf.py:vcf_export_to_file`. bgzip
output stays readable everywhere plain gzip was, so the external-annotation dumps gain indexability too.

### 2. `annotation/backfill_columns.py`

**`resolve_backfill_columns(variant_annotation_version, column_ids, info_overrides, csq_overrides)`**
→ list of resolved targets, each carrying: the VariantGrid column id, the source field to read, whether it
comes from INFO or CSQ, the formatter, and the annotation models that hold the column
(`get_model_fields(VariantAnnotation)` / `VariantTranscriptAnnotation` — `cosmic_count` is variant-level, so
only the representative-transcript table). Names that no registry entry claims for this version's
`columns_version` / `vep_version` are reported as an error naming the column, so a typo or a
version-inappropriate column fails before anything is dumped.

**`dump_annotated_variants(variant_annotation_version, output_filename, …)`**
Wraps `get_variants_qs_for_annotation(annotation_version, annotated=True)` — the variants that already hold a
`VariantAnnotation` row for this version — and writes them with `write_qs_to_vcf`, so the dump carries
`variant_id` in INFO exactly as the pipeline's own dumps do. Options:

- `--only-missing` restricts to rows where every target column is null (the common backfill shape, and much
  smaller than a full dump)
- `--min-variant-id` / `--max-variant-id` split a large dump into pieces that can be annotated in parallel

**`import_backfill_vcf(variant_annotation_version, vcf_filename, targets, …)`**
Reads with `cyvcf2`, and per record pulls `variant_id` plus each target's source field. cyvcf2 hands back an
int/float/str for a single value and a tuple for multi-valued INFO, while the registry formatters take the
`&`-joined string VEP produces — so values are normalised to that string form before the formatter runs, and
one code path serves both bcftools and VEP input. Rows accumulate into an in-memory CSV buffer and flush
through `sql_copy_csv_file` (CSV format, where an unquoted empty field is NULL) every `--batch-size` records.
Progress logs rows read / rows changed per batch.

`--dry-run` runs the COPY and a `SELECT count(*)` against the same join predicate, reporting how many rows
would change, and rolls back.

### 3. `annotation/management/commands/annotation_backfill_columns.py`

Thin command wrapper in the style of `annotation_external`: mutually exclusive `--dump` / `--import`,
`--genome-build` (required), `--columns`, `--info`, `--csq`, `--only-missing`, `--min-variant-id`,
`--max-variant-id`, `--batch-size`, `--clear-missing`, `--dry-run`. The version defaults to the build's ACTIVE
`VariantAnnotationVersion`, with `--variant-annotation-version` to name another.

### 4. Tests

- column resolution: a registry column resolves to its `vep_info_field` and formatter; an `--info` override
  wins; a column outside the version's `columns_version` is reported
- value handling: a cyvcf2 tuple normalises and formats the way the pipeline's `&`-joined string does
  (`(12, 7)` → `12` via `format_pick_highest_int`), and an absent value is left alone by default and NULLed
  under `--clear-missing`
- round trip: dump a small annotated version, hand back a VCF carrying a `cosmic_count` INFO field, import,
  and assert the column holds the formatted values while its neighbours are untouched

## Running it for COSMIC (#1673)

```bash
manage.py annotation_backfill_columns --dump --genome-build GRCh38 \
    --columns cosmic_count,cosmic_legacy_id --only-missing --output cosmic_backfill.vcf.gz
tabix -p vcf cosmic_backfill.vcf.gz
bcftools annotate -a Cosmic_GenomeScreensMutant_v101_GRCh38.vcf.gz \
    -c INFO/GENOME_SCREEN_SAMPLE_COUNT,INFO/LEGACY_ID \
    -O z -o cosmic_backfill.annotated.vcf.gz cosmic_backfill.vcf.gz
manage.py annotation_backfill_columns --import --genome-build GRCh38 \
    --columns cosmic_count,cosmic_legacy_id \
    --info cosmic_count=GENOME_SCREEN_SAMPLE_COUNT --info cosmic_legacy_id=LEGACY_ID \
    --vcf cosmic_backfill.annotated.vcf.gz
```

Run it against a version whose annotation has finished: variants annotated after the dump take their values
from the live pipeline config, so those are the ones the config fix below covers.

## Follow-up

The other half of #1673 is the forward fix — teaching the pipeline the v101 field name so newly annotated
variants get a count. That is a `VEPColumnDef` for `GENOME_SCREEN_SAMPLE_COUNT` alongside the existing `CNT`
(≤ v97) and `SAMPLE_COUNT` (v99) entries, gated on the COSMIC release the way `gnomad4_minor_version` gates
gnomAD's, plus pointing `vep_config["cosmic"]` at the v101 VCFs. Independent of this plan, and worth doing
first so the backfill and the pipeline agree on where the number comes from.
