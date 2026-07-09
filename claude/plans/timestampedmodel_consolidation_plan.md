# Plan: Consolidate on `model_utils.TimeStampedModel`, remove `django_extensions`

Issue: sacgf/variantgrid_private#3437 — VariantGrid uses two `TimeStampedModel` implementations.

## Decision

Standardise every model on `model_utils.models.TimeStampedModel` and eliminate
`django_extensions.db.models.TimeStampedModel`. `django_extensions` is used for
**nothing else** in the codebase (all references are the timestamp mixin or its
`CreationDateTimeField` / `ModificationDateTimeField`, the latter only surviving
inside historical migrations), so migrating the mixin makes the whole package
removable. `model_utils` stays regardless — it is load-bearing via
`InheritanceManager` (8+ model files).

The one behaviour `django_extensions` provides that `model_utils` does not — the
`update_modified=False` opt-out — is replaced by queryset (`.update(...)`) writes,
which bypass `save()`/`pre_save` entirely and so never touch `modified` under
either base. This is the sanctioned "don't bump modified" path going forward.

## Key behavioural change to be aware of

`model_utils.TimeStampedModel.save()` force-adds `modified` to any
`save(update_fields=[...])` call. `django_extensions` does not. So after migration,
**every partial `save(update_fields=[...])` on a formerly-django_extensions model
would also persist `modified = now()`**, where previously it left `modified`
untouched (the attribute was set in `pre_save` but never written because it was
absent from `update_fields`).

**Requirement: preserve current behaviour — a save that does not bump `modified`
today must not bump it after migration.** `modified` is surfaced/used
(`classification_stats.py:105` orders by `-modified`, admin "Modified" columns,
export last-modified basis), so a spurious bump re-orders dashboards and can
re-export a record.

The faithful preservation mechanism is a **queryset `.update(...)`**: it compiles
to raw SQL, never calls `save()`/`pre_save`, and so never touches `modified` under
either base — exactly the workaround James identified. Every no-bump partial save
on a migrating model must therefore be converted to `.update(...)` **before** the
base class changes (Step 1). Every such site across all migrating apps has been
audited and enumerated in Step 1 (classification in 1a–1c, everything else in 1d);
there is no open audit remaining.

## Step 1 — Convert every no-bump save on migrating models to queryset updates

Do this **before** the base-class swap. The `update_modified=False` sites otherwise
raise `TypeError` under `model_utils`; the partial `update_fields` sites otherwise
silently start bumping `modified`.

### Classification-app analysis (verified)

Findings that scope the work:

- **`ClassificationModification` is `models.Model`, NOT `TimeStampedModel`**
  (`classification/models/classification.py:2226`). It has no `modified` field, so
  **all `ClassificationModification.save(update_fields=[...])` sites need no
  change** — leave them as `save()`. This covers: `classification.py:1491`,
  `classification.py:2493`, `evidence_key_rename.py:233`, `change_owner.py:55`,
  `classification_db_refs.py:58`, `classification_history_censor.py:148`,
  `classification_populate_sort_orders.py:15`, `evidence_key_to_unit.py:55`.
- **`qs.update()` is a faithful substitute for `Classification` and grouping
  partial saves** — nothing non-create depends on `save()` there:
  - The only `post_save` receiver on `Classification`
    (`classification/signals/classification_hooks_share_flags.py:23`) is guarded by
    `if created:` — a no-op on updates of existing rows.
  - `Classification.save()`'s only side effect (`fix_permissions`,
    `classification.py:1058`) is create-only (`self.id is None`).
  - `ClassificationGrouping` / `AlleleOriginGrouping` have no `save()` override and
    no `post_save` receiver.
  - There are no sender-less `post_save`/`pre_save` receivers in the codebase.

### 1a. The two `update_modified=False` sites (must convert — would `TypeError`)

`classification/models/classification_inserter.py:359` — only `last_import_run` /
`last_source_id` are mutated on `record` (set at lines 236–239):
```python
# replace  record.update_modified = False / record.save()
Classification.objects.filter(pk=record.pk).update(
    last_import_run=record.last_import_run,
    last_source_id=record.last_source_id,
)
```

`classification/management/commands/fix_variant_matching.py:154,157` — one-off
command; both saves persist the `allele_info` linkage set by `ensure_allele_info()`
/ `update_allele_info_from_classification()`:
```python
# replace each  c.save(update_modified=False)
Classification.objects.filter(pk=c.pk).update(allele_info=c.allele_info)
```
`allele_info` is the complete field set: `ensure_allele_info_with_created`
(`classification.py:843`) mutates only `self.allele_info` on the row;
`update_allele_info_from_classification` saves the `ImportedAlleleInfo` object
separately and changes no other `Classification` column.

### 1b. `Classification` partial `save(update_fields=[...])` sites (must convert — would bump `modified`)

Each becomes `Classification.objects.filter(pk=<obj>.pk).update(<field>=<obj>.<field>, …)`:

| Location | Current | Fields to `.update()` |
|---|---|---|
| `classification/models/classification.py:1494` (`patch_history`) | `self.save(update_fields=['evidence'])` | `evidence=self.evidence` |
| `classification/models/condition_text_matching.py:1170` | `classification.save(update_fields=['condition_resolution'])` | `condition_resolution=classification.condition_resolution` |
| `classification/models/clinical_context_utils.py:99` | `classification.save(update_fields=['clinical_context'])` | `clinical_context=classification.clinical_context` |
| `classification/signals/classification_hooks_pending_flags.py:13` | `c.save(update_fields=["summary"])` | `summary=c.summary` |
| `classification/evidence_key_rename.py:227` | `vc.save(update_fields=["evidence"])` | `evidence=vc.evidence` |

Management commands (one-off; convert so re-runs stay behaviour-identical):
`change_owner.py:59` (`evidence`, `user`), `classification_db_refs.py:52`
(`evidence`), `classification_history_censor.py:135` (`evidence`),
`classification_set_allele_context.py:18` (`allele_origin_bucket`),
`classification_groupings.py:28` (`summary`), `fix_legacy_labs.py:37` (`lab`).

### 1c. `ClassificationGrouping` / `AlleleOriginGrouping` `dirty` saves (must convert)

In `classification/models/classification_grouping.py` — lines 246, 248, 294, 372,
all `save(update_fields=["dirty"])`:
```python
# self.save(update_fields=["dirty"])  ->
type(self).objects.filter(pk=self.pk).update(dirty=True)
# self.allele_origin_grouping.save(update_fields=["dirty"])  ->
AlleleOriginGrouping.objects.filter(pk=self.allele_origin_grouping.pk).update(dirty=True)
```

### 1d. Non-classification apps (audited — only two sites need conversion)

Every other `save(update_fields=[...])` site was audited against the same criteria.
Only **two** sit on a migrating `TimeStampedModel`; the rest target models whose
base chain terminates at plain `models.Model` (via `DataArchiveMixin`,
`SubVersionPartition`/`RelatedModelsPartitionModel`, `AbstractVariantAnnotation`) and
so have **no `modified` column** — the base swap cannot affect them and they stay
as `save(update_fields=...)`.

**Convert:**

`snpdb/common_variants.py:28` — `CohortGenotypeCommonFilterVersion` (django_extensions
TimeStampedModel; no `save()` override, no signal receiver):
```python
CohortGenotypeCommonFilterVersion.objects.filter(pk=common_filter.pk).update(
    additional_gnomad_versions=additional_gnomad_versions,
)
```

`flags/models/models.py:611` — inside `FlagsMixin.flag_collection_safe`, a shared
mixin method. Its concrete subclasses are `Allele` (plain `models.Model`),
`ClinicalContext` and `Classification` (both **django_extensions** TimeStampedModel).
No subclass is a model_utils TimeStampedModel, so none bumps `modified` here today;
the polymorphic form preserves that for all of them (and is a no-op behaviour change
for `Allele`). The property only runs on persisted instances, so `pk` is set:
```python
type(self).objects.filter(pk=self.pk).update(flag_collection=flag_collection)
```

**No change (verified — target model has no `modified` field):**
`annotation/models/models.py:818` (`VariantAnnotationVersion`),
`annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:176`,
`annotation/management/commands/fix_columns_version4_damage_counts.py:47`,
`fix_annotation_sv_overlaps.py:41`, `fix_columns_version2_damage_counts.py:94`,
`fix_historical_spliceai_max_ds.py:40`, `snpdb/archive.py:69`,
`snpdb/partition_archive.py:94,102,117`,
`snpdb/tasks/partition_archive_tasks.py:35,49,59,64,122`,
`snpdb/management/commands/migrate_common_filter_gnomad_versions.py:93`.

Sites that already include `"modified"` in `update_fields` (e.g.
`genes/panel_app.py:46,163`) intend to bump and need no change. No sender-less
`post_save`/`pre_save` receivers exist in the repo, so no migrating model has a
hidden update-time signal side effect.

## Step 2 — Migrate the 32 model modules to `model_utils`

In each file below, change
`from django_extensions.db.models import TimeStampedModel`
to
`from model_utils.models import TimeStampedModel`
(leaving the `class Foo(TimeStampedModel)` declarations untouched):

```
analysis/models/models_analysis.py
analysis/models/models_karyomapping.py
analysis/models/models_variant_tag.py
analysis/models/nodes/analysis_node.py
annotation/models/models.py
annotation/models/models_citations.py
annotation/models/models_cohort_stats.py
annotation/models/models_sample_stats.py
classification/models/classification.py
classification/models/classification_grouping.py
classification/models/clinical_context_models.py
classification/models/discordance_models.py
classification/models/evidence_key.py
flags/models/models.py
genes/models.py
pathtests/models.py
patients/models.py
seqauto/models/models_seqauto.py
snpdb/models/models.py
snpdb/models/models_clingen_allele.py
snpdb/models/models_cohort.py
snpdb/models/models_cohort_stats.py
snpdb/models/models_dbsnp.py
snpdb/models/models_jobs_control.py
snpdb/models/models_somalier.py
snpdb/models/models_user_settings.py
snpdb/models/models_variant.py
snpdb/models/models_vcf.py
snpdb/models/models_zygosity_counts.py
sync/models/models.py
sync/models/models_classification_sync.py
upload/models/models.py
```

Some of these files already `import django_extensions` elsewhere only for the
mixin — after the swap, remove any now-unused `django_extensions` imports the
linter flags.

## Step 3 — Generate migrations

Each migrated model gets an `AlterField` on `created` and `modified`
(`CreationDateTimeField` → `AutoCreatedField`, `ModificationDateTimeField` →
`AutoLastModifiedField`). These are **Python-level** field defaults, so the
generated SQL is a no-op / state-only change — no data or column change.

```bash
python3 manage.py makemigrations
python3 manage.py migrate --keepdb   # verify it applies clean
```

Review the generated migrations to confirm they are field-type-only and touch no
DB columns unexpectedly.

## Step 4 — Add the lint guard banning future `django_extensions` mixin imports

`config/pylint3.rc` has an empty `deprecated-modules=` (line 350). Pylint's
`deprecated-module` (W0402) fires on the module named in a `from … import …`, so:

```ini
deprecated-modules=django_extensions.db.models
```

This flags any new `from django_extensions.db.models import TimeStampedModel`
without a custom checker. Target `django_extensions.db.models` **only** — do not
ban `django_extensions.db.fields` yet, since historical migrations still reference
`CreationDateTimeField` / `ModificationDateTimeField` until the squash (Step 5).
Migrations import from `db.fields`, not `db.models`, so this rule leaves the
migration tree clean.

Confirm `deprecated-module` (W0402) is not in the `disable=` block (line 407+); if
it is, enable it. Then run `./scripts/linting/run_pylint.sh` and confirm the
migrated tree is clean (should be zero hits once Step 2 is complete).

## Step 5 — Follow-up after the migration squash (note, not this PR)

`django_extensions` cannot be dropped from `requirements` or removed from
`INSTALLED_APPS` (`variantgrid/settings/components/default_settings.py:690`) while
historical migrations still deconstruct `django_extensions.db.fields.*` (115
references). Once the planned migration squash collapses those references:

1. Remove `'django_extensions'` from `INSTALLED_APPS`.
2. Remove `django_extensions` from the project requirements.
3. Extend the lint ban to `django_extensions` as a whole (or `django_extensions.db.fields`).

Add a `ManualOperation` / note in the squash PR so this cleanup is not forgotten.

## Testing

- `python3 manage.py makemigrations --check` — no further model/state drift.
- `python3 manage.py test --keepdb` — full suite; pay attention to classification
  import, sync, and export tests (the `modified`-on-partial-save behaviour change).
- Targeted: exercise the classification patch path that previously used
  `update_modified=False` (Step 1a) and confirm `modified` is unchanged after an
  import-run-only update.
- `./scripts/linting/run_pylint.sh` — clean, and the new rule trips on a
  deliberately-added test import.
