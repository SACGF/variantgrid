# Plan: Store both raw and masked SpliceAI in one annotation version (#1579)

## Background & decision

SpliceAI ships two precomputed score sets: **raw** and **masked**. The masked files
zero out splice-site changes Illumina flags as typically less pathogenic and are the
version Illumina recommends for interpretation; raw reports everything. Which to trust
is a lab call (research vs diagnostic), so #1579 wants **both** stored per annotation
version, with a setting to show/hide each on the variant details page and the
damage/splicing analysis node.

**Decision (settled in scoping):**

- Store both as **parallel flat column sets** — `spliceai_raw_*` and `spliceai_masked_*`
  — not a JSON/array blob. The individual sub-fields are wired into the flat-column
  machinery (selectable grid columns, classification autopopulate via flat ORM paths,
  export, template), and there have only ever been two flavours in 5+ years, so there is
  no future-proofing payoff to justify a JSON refactor.
- Historically all real data is **raw**, and masked was only ever deployed to a test
  server. So **rename** the existing `spliceai_*` columns to `spliceai_raw_*` (a cheap
  Postgres metadata-only `RENAME COLUMN`, index preserved, prod data correct under its
  true name) and **add** a new `spliceai_masked_*` set.
- This rides a **columns_version bump (4 → 5)** and therefore a new
  `VariantAnnotationVersion` + reannotation. The test-server masked annotation is dropped
  and redone.

## Naming convention

| Current                     | Becomes (raw)                    | New (masked)                        |
|-----------------------------|----------------------------------|-------------------------------------|
| `spliceai_pred_dp_{ag,al,dg,dl}` | `spliceai_raw_pred_dp_*`     | `spliceai_masked_pred_dp_*`         |
| `spliceai_pred_ds_{ag,al,dg,dl}` | `spliceai_raw_pred_ds_*`     | `spliceai_masked_pred_ds_*`         |
| `spliceai_max_ds` (indexed) | `spliceai_raw_max_ds` (indexed)  | `spliceai_masked_max_ds` (indexed)  |
| `spliceai_gene_symbol`      | `spliceai_raw_gene_symbol`       | `spliceai_masked_gene_symbol`       |
| `VariantAnnotationVersion.spliceai` | `…spliceai_raw`          | `…spliceai_masked`                  |

(Keeping the `_pred_` infix matches the existing names and the `SpliceAI_pred_*` VEP
source fields. The two indexed `*_max_ds` columns stay real scalar columns — DamageNode
filters on them.)

## The VEP execution (single pass, both flavours via `--custom`) — DECIDED

The Ensembl **SpliceAI plugin emits fixed field names** (`SpliceAI_pred_DS_AG`, …) with
no prefix option, so it cannot be invoked twice. Instead, **both** flavours are annotated
as VEP `--custom` VCF sources with **distinct `short_name` prefixes** (`SpliceAI_raw`,
`SpliceAI_masked`) — symmetric, no plugin/custom asymmetry, one VEP pass. The
`VEPPlugin.SPLICEAI` path is retired for columns_version ≥ 5 (kept only for historical
column-def records on pre-5 VAVs, which are never re-annotated).

This reuses the existing custom machinery (`_get_custom_params_list`,
`annotation/vep_annotation.py:36`; same pattern as gnomAD/COSMIC, and gnomAD_SV which
already calls `--custom` more than once). The SpliceAI score files are VCFs whose INFO is
`SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`, so each custom
source yields one pipe-delimited string the bulk inserter splits into the eight sub-fields
+ symbol (then computes `*_max_ds`).

Each flavour has **two files** (snv + indel) like the current plugin config, so that's up
to four `--custom` registrations. A variant is either an SNV or an indel, so its snv/indel
files never both match — give a flavour's snv and indel customs the **same `short_name`**
(`SpliceAI_raw` / `SpliceAI_masked`) so each flavour surfaces as one unified field;
confirm VEP accepts the shared short_name during implementation (else use distinct
names and coalesce in the inserter).

Implementation check (not a blocker, just validate): confirm `--custom type=exact`
reproduces the precomputed lookup on a handful of known SNV + indel variants.

## Change list by area

### 1. Model — `annotation/models/models.py`
- **Fields** (`:1457-1466`): rename the 10 `spliceai_*` fields to `spliceai_raw_*`; add the
  10 parallel `spliceai_masked_*` fields (same types; both `*_max_ds` keep `db_index=True`).
- **`SPLICEAI_DS_DP`** (`:1513`): split into `SPLICEAI_RAW_DS_DP` and
  `SPLICEAI_MASKED_DS_DP` (or one dict keyed by flavour). All consumers below iterate this.
- **`has_spliceai()`** (`:1748`): becomes flavour-aware — e.g. `has_spliceai_raw()` /
  `has_spliceai_masked()` (plus a combined `has_spliceai()` for "either").
- **`highest_spliceai()`** (`:1752`): parameterise by flavour.
- **`backfill_spliceai_max_ds()`** (`:1761`): parameterise per flavour (writes
  `spliceai_raw_max_ds` / `spliceai_masked_max_ds`).

### 2. VariantAnnotationVersion — `annotation/models/models.py`
- **`spliceai` field** (`:670`): rename → `spliceai_raw`; add `spliceai_masked`.
- **`uses_raw_spliceai`** (`:732`): both flavours now coexist, so this property's meaning
  changes. Replace with `has_raw_spliceai` / `has_masked_spliceai` presence checks; the
  DamageNode warning logic (below) keys off which flavour the user chose, not "is raw".
- **`backfilled_spliceai_max_ds`** (`:704`): on a fresh columns_version-5 VAV both max_ds
  columns are populated at insert, so this defaults True as before. The historical fallback
  branch in DamageNode still references the old single column for pre-5 VAVs (see §6).

### 3. VEP command — `annotation/vep_annotation.py`
- **Drop the plugin for columns_version ≥ 5** (`:201`): the `SpliceAI,snv=…,indel=…`
  plugin invocation is no longer added. Both flavours go through the `--custom` loop
  (`:217-231`) via their new `VEPCustom` members (§5), emitting prefixes `SpliceAI_raw`
  and `SpliceAI_masked`.
- **`_spliceai_label()`** (`:281`) + **VAV kwargs** (`:403-408`): produce **both**
  `spliceai_raw` and `spliceai_masked` labels from the two file configs.

### 4. VEP column defs — `annotation/vep_columns.py:450-522`
- Set `max_columns_version=4` on the existing 9 plugin-sourced SpliceAI `VEPColumnDef`s
  (they stay valid for pre-5 VAVs, now mapping to the renamed `spliceai_raw_*` grid
  columns). Add **two** new custom-sourced sets gated `min_columns_version=5`, both with
  `source_field_has_custom_prefix=True` and `vep_custom=` their member — one prefixed
  `SpliceAI_raw` → `spliceai_raw_*`, one `SpliceAI_masked` → `spliceai_masked_*`. Each
  parses the pipe-delimited custom field into its sub-fields.

### 5. Enums — `annotation/models/models_enums.py`
- **`VEPCustom`** (`:116`): add `SPLICEAI_RAW` and `SPLICEAI_MASKED` members.
  `VEPPlugin.SPLICEAI` (`:112`) is retained only for historical pre-5 column defs.

### 6. Bulk inserter — `annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py`
- **`_add_spliceai_max_ds`** (`:587`): compute **both** `spliceai_raw_max_ds` and
  `spliceai_masked_max_ds` from their respective DS fields.
- Parse both custom pipe-delimited fields (`SpliceAI_raw`, `SpliceAI_masked`) into their
  `spliceai_{raw,masked}_*` sub-fields.

### 7. Settings — `variantgrid/settings/components/annotation_settings.py`
- `vep_config` (GRCh37 `:123-124`, GRCh38 `:170-171`, T2T): replace the single
  `spliceai_snv`/`spliceai_indel` with four keys —
  `spliceai_raw_snv/indel` + `spliceai_masked_snv/indel` — matching the new `VEPCustom`
  members (lower-cased) so the custom loop picks them up.
- Bump `ANNOTATION_VEP_COLUMNS_VERSION` 4 → 5 (`:60`) and update the legacy
  `pin_annotation_to_columns_version_3()` rollback path accordingly.

### 8. Analysis / DamageNode — `analysis/models/nodes/filters/damage_node.py:246-297`
- Add a **flavour selector** for the splice filter (raw / masked, default **masked** =
  recommended). The `splice_min` query (`:270-297`) targets the chosen flavour's
  `spliceai_*_max_ds` (with the per-DS-field fallback for pre-5 VAVs that only have the
  old single column — keep reading the renamed `spliceai_raw_max_ds` there).
- Rework the warning (`:246-253`): warn only when the user selects the **raw** flavour
  (Illumina recommends masked), instead of warning per-VAV.
- This is the node half of the issue's show/hide requirement.

### 9. Classification autopopulate — `classification/autopopulate_evidence_keys/evidence_from_variant.py:452-496`
- `_get_spliceai_summary` (`:452`) and `va_fields_for_summaries` (`:496`) iterate
  `SPLICEAI_DS_DP`. Decide which flavour feeds the existing SpliceAI evidence key
  (default masked) — or emit both summaries if there's an appetite for a second key.

### 10. Export — `classification/management/commands/classification_export_w_annotations.py:88-99`
- `has_spliceai()` + the four `spliceai_pred_ds_*` cells: point at the chosen flavour
  (default masked) or widen to both.

### 11. Variant details template — `variantopedia/templates/variantopedia/variant_details.html:983-1007`
- Render raw and masked blocks, each gated by its `has_spliceai_*` and by the show/hide
  setting (§12). The JS that fills the `spliceai_pred_ds_*` spans by id needs raw/masked
  id variants.

### 12. Show/hide setting (the core deliverable of #1579)
- Add settings (e.g. `VARIANT_DETAILS_SHOW_SPLICEAI_RAW` / `…_MASKED`, or one enum) in
  `default_settings.py`, defaulting to the lab's preference (masked on, raw on/off).
  Gate the variant-details blocks (§11) and the DamageNode flavour options (§8).
  Check existing show/hide settings conventions for the right pattern before adding.

### 13. Grids — `annotation/grids.py:110`
- The single `"spliceai"` RichColumn → decide raw/masked (default masked) or two columns.

### 14. Misc references
- `VariantGridColumn` rows for the renamed/new fields (see fixture/migration precedents:
  `snpdb/migrations/0023…`, `0086_update_spliceai_columns.py`,
  `annotation/migrations/0027…`). Each sub-field needs its column registration updated/added.

## Migrations & rollout

1. **Schema migration**: `RenameField` ×10 (`spliceai_*` → `spliceai_raw_*`) — metadata-only,
   instant even on the partitioned VariantAnnotation table; `AddField` ×10 for
   `spliceai_masked_*`. Rename `VariantAnnotationVersion.spliceai` → `spliceai_raw`, add
   `spliceai_masked`.
2. **VariantGridColumn** data migration for renamed + new columns.
3. Drop the test-server masked `VariantAnnotationVersion`(s) and reannotate under
   columns_version 5.
4. No prod backfill needed — existing prod rows are correctly raw under the renamed
   columns; masked populates on the next annotation run.
5. Bump `CACHE_VERSION` if any cached column metadata is affected.

## Testing

- `annotation/tests/test_spliceai.py`: extend for both flavours' `*_max_ds` computation.
- `analysis/tests/test_damage_node_spliceai_raw_warning.py`: update — warning now keyed to
  the selected flavour; add masked-default coverage.
- New: bulk-inserter test that a VEP VCF with plugin(masked) + custom(raw) populates both
  column sets correctly (Approach A parsing).
- Autopopulate/export tests pointing at the default flavour.
- If Approach A: a validation check comparing custom-raw allele matching vs plugin output
  on known variants (gate the A-vs-B decision).

## Open questions / decisions to confirm

1. **Default flavour** for the single-valued consumers (evidence key, export, grid): masked
   recommended. Confirm with the labs.
2. **Show/hide granularity** — per-flavour booleans vs a single setting; and whether the
   default hides raw by default (diagnostic) or shows both.

**Decided:** both flavours via `--custom` with distinct `short_name` prefixes
(`SpliceAI_raw` / `SpliceAI_masked`), single VEP pass, `VEPPlugin.SPLICEAI` retired for
columns_version ≥ 5.
