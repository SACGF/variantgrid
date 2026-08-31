# Clean up `import_processing` (#928)

Written by Claude Fable 5 (claude-fable-5), 2026-08-31

## What the issue asks

1. "A lot of processes don't clean up import_processing" – find out which, and why.
2. Consolidate `VCF_IMPORT_DELETE_TEMP_FILES_ON_SUCCESS` / `VARIANT_ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS`
   into one setting.
3. (comment) `UploadPipelineFinishedTask` and `pipeline_success_task` look like duplicates – roll into one?
   "One of the reasons for no cleanup is that UploadPipelineFinishedTask is setting status to SUCCESS while
   pipeline success task only does anything if status = PROCESSING."

Item 2 was done in 7f7a62d67 (both settings became `IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS`,
now `True` by default in `variantgrid/settings/components/default_settings.py:449`). Item 3 is the main
live bug and the comment's diagnosis is exactly right – details below.

## Diagnosis

### What is on disk (this dev box, `data/import_processing`, 5,276 entries, ~30 GB)

| prefix | dirs | size | owner state | why it is still there |
|---|---|---|---|---|
| `pipeline_*` | 258 | 7.1 GB | 133 SUCCESS, 124 row deleted, 1 ERROR | §A (success path never runs), §B (nothing on delete) |
| `annotation_run_*` | 971 | 23 GB | all FINISHED, all pk ≤ 4472, none after 2026-07-06 | historical: setting was `not DEBUG` until f44fd2db3. Leak has stopped; needs a one-off sweep only |
| `clingen_allele_registry_<uuid>` | 3,856 | 16 MB | – | all **empty**: §C |
| `liftover_*` / `manual_variants_*` / `classification_import_*` | 114 / 16 / 45 | ~2 MB | pipelines SUCCESS | §D (generated input VCF has no lifecycle hook) |
| `gene_annotation_*` | 3 | 9 MB | – | §E management commands leave their COPY CSV |
| `test/<uuid>`, plus `pipeline_4919` etc. with no DB row | 26 + | small | – | §F tests write into the real dir |

### A. `UploadPipelineFinishedTask` defeats `pipeline_success_task` (the comment's hunch, confirmed)

`schedule_pipeline_stage_steps` (`upload/tasks/vcf/import_vcf_step_task.py:207-213`) builds the FINISH
stage as a **chain**: `[<FINISH-dependent steps in sort order>..., pipeline_success_task]`. For every
factory that inherits `AbstractVCFImportTaskFactory.get_finish_task_classes()` (`[UploadPipelineFinishedTask]`,
`abstract_vcf_import_task_factory.py:50`) the chain is therefore:

1. `UploadPipelineFinishedTask.process_items` (`import_vcf_tasks.py:152`) – sets `status = SUCCESS`, saves.
2. `pipeline_success_task` (`import_vcf_step_task.py:233`) – guarded by `if status == PROCESSING`; it is
   now SUCCESS, so it returns without calling `UploadPipeline.success()`.

`UploadPipeline.success()` (`upload/models/models.py:226`) is the only place that (a) calls
`remove_processing_files()`, (b) records `processing_seconds_*`, (c) fires the `import_*_success` event.
None of that happens for those pipelines.

DB proof (local): of 178 SUCCESS pipelines, 156 have `processing_seconds_wall_time IS NULL`. Every one of
those has `UploadPipelineFinishedTask` among its FINISH steps (Liftover, Manual Variant Entry,
Insert-variants-only, Variant Tags, TSO500 …). The 22 with timings are precisely the factories whose FINISH
list lacks it: ClinVar (`[ImportClinVarSuccessTask]`), Patient Records (no FINISH steps), and genotype VCF
(`settings.FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = []`). Those three prove `pipeline_success_task` on its own
closes a pipeline correctly – FINISH is scheduled by `check_pipeline_stage` whether or not any
FINISH-dependent steps exist.

So the answer to the comment is yes: `UploadPipelineFinishedTask` is redundant, and its presence is the
bug. The `status == PROCESSING` guard in `pipeline_success_task` must stay – FINISH can legitimately be
scheduled more than once (two DATA_INSERTION steps finishing together each call `check_pipeline_stage`;
`schedule_pipeline_stage_steps` de-duplicates the *steps* via `start_date` but appends
`pipeline_success_task` unconditionally), and the guard is what makes the second call a no-op.

### B. Deleting an `UploadPipeline` row leaves its directory

124 `pipeline_*` dirs have no row. There is no `post_delete` for `UploadPipeline` (only
`pre_delete_uploaded_vcf`). `annotation/signals/annotation_run_cleanup.py` already establishes the pattern
for AnnotationRun; UploadPipeline needs the same.

### C. `get_import_processing_dir()` creates the directory on every call

`library/django_utils/django_file_utils.py:8` – `mk_path` inside the getter. `ClinGenAlleleRegistryAPI.__init__`
(`snpdb/clingen_allele_api.py:73`) computes `api_failure_output_filename` eagerly for every instance, so each
of the ~7 `.instance()` call sites mints an empty `clingen_allele_registry_<uuid>` dir per call. The file is
only ever written in the `except` branch of `_get_or_post` (line 122). Same getter also means every
`remove_processing_files()` recreates the dir before deleting it (harmless but backwards).

### D. Generated input VCFs

`snpdb/liftover.py:95`, `annotation/manual_variant_entry.py:80`, `classification/classification_import.py:122`
write a VCF into their own `<prefix>_<pk>` dir, point `FileUpload.path` at it and run a pipeline. Nothing
removes that dir. (Variant Tags and TSO500 fusions write theirs *inside* `pipeline_<pk>` via a
`pre_vcf_task`, so they are regenerated on retry and cleaned with the pipeline – the right shape, but
moving liftover/manual/classification onto that shape is more than this issue needs.)

Trade-off to accept: once removed, "Retry import" on an already-successful liftover / manual entry /
classification import is unavailable. `view_upload_pipeline` already handles a missing file (warning
"File does not exist on disk, cannot reload", `allow_retry_import=False`), and ERROR pipelines keep their
input, which is the case retry exists for.

### E. Management commands

`annotation/management/commands/gene_annotation.py:539` and `human_protein_atlas_import.py:103` COPY from a
CSV under `gene_annotation_<pk>` / `human_protein_atlas_<pk>` and leave it. `genes/models.py:2231` (gene
coverage) cleans up but gates on `not settings.DEBUG` instead of the setting.

### F. Tests

`annotation/fake_annotation.py:42` puts fake-annotation scratch under the *real*
`PRIVATE_DATA_ROOT/import_processing/test/<uuid>` and never removes it; tests that build an `UploadPipeline`
without `override_settings(IMPORT_PROCESSING_DIR=…)` write `pipeline_<test-db pk>` straight into the real
dir (`pipeline_4919`, mtime today, no row in the dev DB).

### Settings

`IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS` (import_processing scratch) and
`ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS` (#1670, `ANNOTATION_VCF_DUMP_DIR` – VEP input/output, AnnotSV
dir) govern different trees with different retention reasoning (a failed VEP run keeps its dump for
investigation and retry-upload-only). Keep both; this plan uses the import_processing one everywhere under
`IMPORT_PROCESSING_DIR`.

## Plan

### 1. One success path for VCF pipelines

- Delete `UploadPipelineFinishedTask` (`upload/tasks/vcf/import_vcf_tasks.py:152-158` and the
  `register_task` at line 282). Remove the import in `abstract_vcf_import_task_factory.py`.
- `AbstractVCFImportTaskFactory.get_finish_task_classes()` returns `[]`.
- `LiftoverImportFactory.get_finish_task_classes()` (`import_task_factories.py:403`) returns
  `[LiftoverCompleteTask]`.
- `pipeline_success_task` stays as the single closer; keep its `PROCESSING` guard. Add a short docstring
  saying it is the only thing that takes a VCF pipeline out of PROCESSING and why the guard exists.
  The comment already on `GeneLevel…get_finish_task_classes` (`import_task_factories.py:168`) remains true.
- `UploadStep` rows on existing deployments that still name
  `upload.tasks.vcf.import_vcf_tasks.UploadPipelineFinishedTask` in `script` belong to finished pipelines
  and are never re-launched (`schedule_pipeline_stage_steps` filters `start_date__isnull=True`), so no data
  migration. A retry of an old pipeline re-creates its steps from the factory.

### 2. Path helpers that only create on write

`library/django_utils/django_file_utils.py`:

```python
def import_processing_dir_path(pk, prefix='pipeline') -> str:      # pure – no mkdir
def get_import_processing_dir(pk, prefix='pipeline') -> str:       # path + mk_path (existing callers)
def get_import_processing_filename(pk, base_filename, prefix='pipeline') -> str   # unchanged
def remove_import_processing_dir(pk, prefix='pipeline'):           # rmtree(ignore_errors=True) of the pure path, logs
```

Use `remove_import_processing_dir` from:
- `UploadPipeline.remove_processing_files()` (`upload/models/models.py:210`) – also drops the bare
  `rmtree` that raises if the dir is gone.
- `BulkVEPVCFAnnotationInserter.remove_processing_files()` (`bulk_vep_vcf_annotation_inserter.py:952`).
- `AnnotationRun.reset_for_retry` (`annotation/models/models.py:1376`) – replace the hard-coded
  `f"annotation_run_{self.pk}"` with the helper. Put the prefix constant
  `ANNOTATION_RUN_IMPORT_PROCESSING_PREFIX = "annotation_run"` in `annotation/annotation_run_files.py`
  (imports only library/snpdb, so both models and the bulk inserter can import it) and have
  `BulkVEPVCFAnnotationInserter.PREFIX` read it.

### 3. Pipeline directory lifecycle

`upload/models/models.py`:
- `UploadPipeline.success()` keeps calling `remove_processing_files()` under the setting (now reachable
  thanks to §1) and additionally calls a new `remove_generated_input_file()` (see §4).
- New `post_delete` receiver for `UploadPipeline` next to `pre_delete_uploaded_vcf` (`models.py:528`):
  `remove_processing_files()` + `remove_generated_input_file()`, unconditional – with the row gone nothing
  can name the directory (same reasoning as `annotation_run_post_delete_handler`). Registering a receiver
  also disables Django's fast-delete path for queryset deletes, which is what makes it fire for
  `FileUpload`/`VCF` cascades.
- ERROR pipelines keep everything on disk (matches #1670's decision for AnnotationRun); retry already
  wipes the pipeline dir before re-running.

### 4. Generated input VCFs (liftover / manual variants / classification import)

- `FileUpload.is_import_processing_scratch` property: `path` is set, `import_source != WEB_UPLOAD`, and
  `os.path.realpath(path)` is under `settings.IMPORT_PROCESSING_DIR`.
- `UploadPipeline.remove_generated_input_file()`: when the file is scratch **and** it is outside the
  pipeline's own dir, `rmtree(dirname(path), ignore_errors=True)`. (Files inside `pipeline_<pk>` – variant
  tags, TSO500 – are already covered by `remove_processing_files`.)
- Called from `success()` (gated on the setting) and from the `post_delete` receiver (unconditional).
  `retry_upload_pipeline` (`upload/uploaded_file_type.py:98`) keeps calling only `remove_processing_files()`
  so a retry still has its input.
- `LiftoverRun.source_vcf` stays as a record of what was written, like AnnotationRun's filename fields.

### 5. ClinGen failure dump

`snpdb/clingen_allele_api.py`: drop the constructor's eager `get_import_processing_filename`. Keep the
`api_failure_output_filename` kwarg for callers/tests that pass one; when it is `None`, compute the path
inside the `except` branch of `_get_or_post` immediately before writing, as
`get_import_processing_filename("failures", f"{uuid4()}.json", prefix="clingen_allele_registry")` – one
shared `clingen_allele_registry_failures/` dir, created only when a failure is actually dumped. Check
`snpdb/tests/utils/mock_clingen_api.py` and `snpdb/tests/test_clingen_allele.py` for anything relying on
the attribute being pre-set.

### 6. Management commands and gene coverage

- `gene_annotation.py::_write_records` and `human_protein_atlas_import.py`: after `sql_copy_csv`, if
  `settings.IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS`, `remove_import_processing_dir(pk, prefix)`.
- `genes/models.py:2231`: gate on `settings.IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS` instead of
  `not settings.DEBUG`, and go through the pure-path helper (`gene_coverage/gene_coverage_collection_<pk>` is
  a nested layout – keep it, just replace the inline `rmtree`).
- `snpdb/tasks/somalier_tasks.py` already cleans up; leave as is.

### 7. Tests write to a throwaway directory

`variantgrid/settings/components/default_settings.py`, inside an existing `if UNIT_TEST:` block:

```python
IMPORT_PROCESSING_DIR = os.path.join(tempfile.gettempdir(), "variantgrid_unit_test", "import_processing")
```

and `VariantGridTestRunner.teardown_test_environment` removes it. `annotation/fake_annotation.py:42` then
builds its per-test dir as `os.path.join(settings.IMPORT_PROCESSING_DIR, "test", str(uuid4()))` so it lands
in the same tree. Existing `override_settings(IMPORT_PROCESSING_DIR=tempdir)` in individual tests keep
working.

### 8. Sweep for existing deployments

New management command `upload/management/commands/import_processing_cleanup.py` (`--dry-run` flag,
prints per-prefix counts and bytes). It walks `settings.IMPORT_PROCESSING_DIR` top level and removes:

| entry | remove when |
|---|---|
| `pipeline_<pk>` | no `UploadPipeline` row, or `status == SUCCESS` |
| `annotation_run_<pk>` | no `AnnotationRun` row, or `get_status() == FINISHED` |
| `clingen_allele_registry_<uuid>` | directory is empty |
| `liftover_<pk>` | `UploadedLiftover(liftover_id=pk)` missing, or its pipeline SUCCESS/missing |
| `manual_variants_<pk>` | `UploadedManualVariantEntryCollection(collection_id=pk)` – same rule |
| `classification_import_<pk>` | `UploadedClassificationImport(classification_import_id=pk)` – same rule |
| `somalier_vcf_extract_<pk>` / `somalier_relate_<pk>` | owning row missing or status SUCCESS |
| `gene_annotation_<pk>` / `human_protein_atlas_<pk>` | owning version row exists |
| `test/` | always |

Anything else (unknown prefix, ERROR/PROCESSING owners) is listed and left alone. Add
`upload/migrations/00XX_one_off_import_processing_cleanup.py` with
`ManualOperation.task_id_manage(["import_processing_cleanup"])` and a `test=` that returns True only when
`IMPORT_PROCESSING_DIR` exists and is non-empty, so the upgrade script surfaces it where there is something
to reclaim. The command is worth keeping (developers running with the setting off can run it by hand).

### 9. Tests

- `upload/tests/test_pipeline_success.py`: pipeline in PROCESSING with a FINISH-dependent step; run
  `pipeline_success_task` → status SUCCESS, `processing_seconds_wall_time` set, pipeline dir removed with the
  setting True and kept with it False; a second call is a no-op.
- `post_delete` on `UploadPipeline` removes the pipeline dir and a scratch input dir; a WEB_UPLOAD
  `FileUpload` (path under upload storage) is left alone.
- `retry_upload_pipeline` removes the pipeline dir and keeps the generated input.
- `import_processing_dir_path` creates nothing; `ClinGenAlleleRegistryAPI()` creates nothing.
- `import_processing_cleanup --dry-run` / real run against a tempdir seeded with one entry per row of the
  §8 table (owner present-and-successful, owner present-and-error, owner missing).

Audit at the end per CLAUDE.md – drop any that only restate framework behaviour.

## Files

- `library/django_utils/django_file_utils.py`
- `upload/models/models.py` (success, remove_*, post_delete receiver, `FileUpload.is_import_processing_scratch`)
- `upload/tasks/vcf/import_vcf_tasks.py`, `upload/tasks/vcf/import_vcf_step_task.py`
- `upload/import_task_factories/abstract_vcf_import_task_factory.py`, `import_task_factories.py`
- `upload/management/commands/import_processing_cleanup.py` + `upload/migrations/00XX_…py`
- `annotation/annotation_run_files.py`, `annotation/models/models.py`,
  `annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py`
- `annotation/management/commands/gene_annotation.py`, `human_protein_atlas_import.py`
- `annotation/fake_annotation.py`, `variantgrid/test_runner.py`,
  `variantgrid/settings/components/default_settings.py`
- `genes/models.py`, `snpdb/clingen_allele_api.py`
- `variantgrid/templates/default_templates/changelog.html` – entry for #928
- tests as in §9
