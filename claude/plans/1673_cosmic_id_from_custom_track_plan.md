# COSMIC ID from the custom track as well as VEP's cache (#1673)

## Goal

Populate `VariantAnnotation.cosmic_id` from **both** sources that already carry it, rather than only VEP's
cache:

- `Existing_variation` (`--check_existing`), filtered to `COSV` — what we use today.
- The COSMIC `--custom` track's matched record ID — already present in the CSQ and currently discarded.

Union, never replace. VEP's cache knows COSMIC IDs for variants the GenomeScreensMutant file has no record
of, so replacing would cost coverage; the custom track knows the variants COSMIC matched in this release.
Taking both means a variant with a `cosmic_count` / `cosmic_legacy_id` always has a `cosmic_id` from the
same release, without giving up the cache's breadth.

## Why this is small

**The data is already there.** VEP emits the custom track's matched record ID as the bare `short_name`
column, and we throw it away — `"COSMIC"` sits in `VEP_NOT_COPIED_FIELDS`
(`annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:197`).

Verified on a real annotation run output on `vgaws`
(`/data/annotation_scratch/variantgrid__GRCh38__RefSeq__vep116__cv5__gar104__run33097__standard.vep_annotated.vcf.gz`):

```
CSQ field #171 = COSMIC            populated with COSV IDs
CSQ field #172 = COSMIC_CNT        empty  <- #1673: asked for the v95/97 name against a v99 file
CSQ field #173 = COSMIC_LEGACY_ID  populated
```

Over the first 20,000 records: 4,542 CSQ entries carry a bare `COSMIC` ID, **all 4,542** also carry
`COSMIC_LEGACY_ID`, and **0** carry `COSMIC_CNT`. The custom track always returns ID and legacy ID
together, which is why the backfill runbook measured 0 variants holding a legacy ID but no `cosmic_id`.

No VEP config change is needed. `_get_custom_params_list` (`annotation/vep_annotation.py:42`) builds
`fields` from `cvf.source_field`, so the bare column comes for free from `short_name`; adding a
`VEPColumnDef` with `source_field='COSMIC'` would instead append `COSMIC` to the `fields` list and ask the
COSMIC VCF for an INFO field that does not exist.

## The one subtlety — merge after formatting, not before

`vep_to_db_dict` (`bulk_vep_vcf_annotation_inserter.py:474`) runs a straight-copy loop into `raw_db_data`,
then a cluster of fix-up hooks, then applies per-column formatters into `db_data`:

```python
for source_field, dest_columns in self.source_field_to_columns.items():
    for c in dest_columns:
        if c in model_columns:
            raw_value = vep_transcript_data[source_field]
            if raw_value is not None:
                raw_db_data[c] = raw_value          # <- plain assignment, last write wins
...
self._pick_aloft_values(raw_db_data)
self._pick_promoter_ai_values(raw_db_data)
...
db_data[c] = formatter(value)                        # <- extract_cosmic runs here
```

Two consequences:

1. Stacking a second `VEPColumnDef` on `cosmic_id` does **not** work. The copy loop assigns rather than
   merges, so one source silently wins; and `field_formatters` is keyed by *column*
   (`:325-329`), so `extract_cosmic` — which only makes sense applied to `Existing_variation` — would have
   to serve both sources.
2. At hook time `raw_db_data['cosmic_id']` still holds the **raw** `Existing_variation` string
   (`rs123&COSV456&...`). Only after the formatter loop is it the sorted `&`-joined COSV set. So the merge
   belongs *after* the formatter loop, where both sides are already COSV lists.

## Changes — `annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py`

### 1. Name the CSQ column

Add to `VEPColumns` (`:58`):

```python
COSMIC = "COSMIC"  # --custom short_name column: the matched COSMIC record's COSV ID
```

### 2. Note why it stays in `VEP_NOT_COPIED_FIELDS`

Keep `"COSMIC"` in the list (`:197`) so the straight-copy loop leaves it alone and it does not trip the
unhandled-CSQ-column warning, with a comment pointing at the merge:

```python
"COSMIC",  # merged into cosmic_id by _merge_custom_cosmic_id, not straight-copied
```

### 3. The merge itself

```python
COSMIC_ID_COLUMN = 'cosmic_id'

def _merge_custom_cosmic_id(self, db_data: dict, vep_transcript_data: TranscriptData):
    """ Existing_variation only carries the COSMIC IDs VEP's cache knows; the COSMIC --custom track
        returns the matched record's COSV ID as its bare short_name column. Union them so a variant
        COSMIC matched in this release always has an ID, keeping the cache's broader coverage (#1673). """
    custom_cosmic_id = vep_transcript_data.get(VEPColumns.COSMIC)
    if not custom_cosmic_id:
        return
    cosmic_ids = set(custom_cosmic_id.split(VEP_SEPARATOR))
    if existing := db_data.get(self.COSMIC_ID_COLUMN):
        cosmic_ids.update(existing.split(VEP_SEPARATOR))
    db_data[self.COSMIC_ID_COLUMN] = VEP_SEPARATOR.join(sorted(cosmic_ids))
```

Call it at the end of `vep_to_db_dict`, after the formatter loop and before the return:

```python
if self.COSMIC_ID_COLUMN in model_columns:
    self._merge_custom_cosmic_id(db_data, vep_transcript_data)
return db_data
```

Notes:
- `cosmic_id` is declared on `VariantAnnotation` (`annotation/models/models.py:1762`), not
  `VariantTranscriptAnnotation`, so the `model_columns` guard keeps it out of the per-transcript dict —
  without it we would add a column the transcript insert does not have.
- Empty CSQ fields are already `None` — `empty_to_none` maps `EMPTY_VALUES` (`''`, `'.'`) at
  `:819`/`:967` — so a falsy check is enough.
- `sorted(set(...))` matches what `extract_cosmic` already produces
  (`annotation/vep_field_formatters.py:190`, `get_extract_existing_variation("COSV")`), keeping one
  canonical spelling of the column.
- The custom column can itself be `&`-joined when several COSMIC records hit one position (`num_records=all`,
  `type=exact`), so it is split rather than treated as a scalar.

## Deliberately not doing

- **No `columns_version` bump.** This is a bug fix, it affects one production deployment, and the existing
  data is being backfilled separately.
- **No new `VEPColumnDef`** — see "Why this is small"; it would corrupt the `--custom fields` parameter.
- **No change to `extract_cosmic`** or to `Existing_variation` handling.
- **`cosmic_id` is still never cleared** by this path — the merge only ever adds.

## Testing

Unit-test `_merge_custom_cosmic_id` — it is our merge rule, so it earns its keep:

- custom ID present, `cosmic_id` empty -> custom ID alone
- custom ID present, `cosmic_id` set and disjoint -> both, sorted, `&`-joined
- custom ID already in `cosmic_id` -> unchanged, no duplicate
- custom column `&`-joined -> all values merged
- no custom ID -> `cosmic_id` untouched (including left `None`)
- `cosmic_id` not in `model_columns` -> key never added

**These cannot be run on `vgaws`** — it talks to the production database. Use `python3 -c "import py_compile..."`
here and run the suite on a dev box with `python3 manage.py test --keepdb annotation`.

## Relationship to the backfill

This fixes annotation going forward only. Existing versions pick the same IDs up from §8 of
`/tmp/vg_com_prod_cosmic_backfill.md`, which fills `cosmic_id` gaps from the `INFO/COSMIC_ID` field carried
through the backfill VCF — fill-only, no `--clear-missing`. Expect that gap to be non-zero for the first
time after this backfill, since v101 matches variants v99 never did.
