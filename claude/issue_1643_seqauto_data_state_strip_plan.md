# Issue #1643 — Remove `data_state` from `SeqAutoRecord`

> The filesystem-scan / cluster-job deletion (Phases 1–4b of the original plan) is done and verified.
> This file now covers only the deferred `data_state` strip.

## Background

`data_state` meant "expected file that should exist on disk but doesn't" — a concept from the old
filesystem scan. Under the REST-API model the pipeline POSTs records for files that exist, so both
serializers just hard-write `DataState.COMPLETE` on every create. The field no longer carries
meaning and should be dropped from the `SeqAutoRecord` base.

## The two facts that shape this work

1. **`SeqAutoRecord` is the shared polymorphic base of every sequencing/QC model.** `data_state`
   (plus `hash`, `file_last_modified`, `is_valid`) is a column on that base, inherited by all 14
   API-populated subclasses (`SequencingRun`, `SampleSheet`, `IlluminaFlowcellQC`, `Fastq`, `FastQC`,
   `BamFile`, `Flagstats`, `SingleSampleVCF`, `JointCalledVCF`, `QC`, `QCGeneList`, `QCExecSummary`,
   `QCGeneCoverage`). So this is a **field-strip across the base + all readers/writers + a migration**,
   not a model delete.

2. **The `DataState` enum (`snpdb/models/models_enums.py`) must NOT be deleted.** It is used
   independently by the `genes`/`upload` gene-coverage feature on a *separate* model field,
   `GeneCoverageCollection.data_state` (`genes/models.py:2212`), which has nothing to do with
   `SeqAutoRecord`. Confirmed external users to leave untouched:
   - `genes/tasks/gene_coverage_tasks.py`, `genes/archive.py`, `genes/models.py`,
     `genes/serializers.py`, `genes/admin.py`, `genes/tests/*`
   - `upload/tasks/import_gene_coverage_task.py`

## Why this is its own phase (not risky, but wide)

Failures surface at **request/render time**, not at `manage.py check` or migration time:
- Both REST serializers currently **write** `DataState.COMPLETE` on create. Drop the column but miss
  a writer → `AttributeError` on every API ingest.
- Miss a template badge or grid filter → a QC / sequencing detail page 500s.

The 7-test seqauto suite won't reliably catch these, so exercise the touched pages + endpoints.

## Where `data_state` is read/written (seqauto only — all must be cleared before the migration)

- **Serializers (writers):** `serializers/seqauto_qc_serializers.py`, `serializers/sequencing_serializers.py`
  — stop setting `data_state = DataState.COMPLETE`.
- **Grids:** `grids/qc_data_grids.py` (filters), `grids/sequencing_data_grids.py` (~10 refs).
- **Forms:** `forms.py`.
- **Views:** `views.py`.
- **Model/helpers:** `models/models_seqauto.py`, `qc/sequencing_run_utils.py`,
  `tasks/gold_summary_tasks.py`, `management/commands/reload_qc_gene_coverage.py`.
- **Template tag:** `templatetags/seqauto_record_tags.py` (`record_data_state_helper`).
- **Templates:** `templates/seqauto/tags/record_data_state_helper.html`,
  `tabs/view_qc_exec_summary_tab.html`, `tabs/view_sequencing_run_stats_tab.html`,
  `view_bam_file.html`, `view_joint_called_vcf.html`, `view_qc.html`, `view_sequencing_run.html`,
  `view_single_sample_vcf.html`, `view_unaligned_reads.html`.
- **Tests:** `tests/test_joint_called_vcf.py`, `tests/test_urls.py`.

## Steps

1. Rework/remove the `record_data_state_helper` tag and the templates that render data-state badges.
2. Stop both serializers writing `data_state`; drop the field from grids and forms.
3. Clear the remaining readers in views, `models_seqauto.py`, `qc/sequencing_run_utils.py`,
   `gold_summary_tasks.py`, `reload_qc_gene_coverage.py`.
4. Adjust `tests/test_joint_called_vcf.py` and `tests/test_urls.py`.
5. Migration dropping `data_state` from `SeqAutoRecord`. Evaluate `hash`, `file_last_modified`,
   `is_valid` at the same time — drop any that are also dead; if dropped, also remove
   `SeqAutoRecord.get_file_last_modified`, `last_modified_datetime`, the `validate()` no-op, and
   simplify `SeqAutoRecord.save()`.
6. Keep the `DataState` enum and every `genes`/`upload` reference to it.

Optional tidy-up while here: the scan-only `get_path_from_*` helpers left after the earlier phases —
low payoff (harmless dead-ish code) and the `get_params`/path chain the API needs is intertwined, so
only if it stays clean.

## Verification

- `python3 manage.py test --keepdb seqauto` and `python3 manage.py makemigrations --check`.
- `grep -rn "data_state" seqauto/` returns nothing (the enum import stays only in `snpdb`).
- Confirm `genes`/`upload` gene-coverage still reference `DataState` and pass their tests.
- Smoke-test the REST bulk-create endpoints (sequencing_files, qc_exec_summary, qc_gene_coverage,
  qc_gene_list) and the sequencing-run / QC / BAM / VCF detail pages — the writers and badge
  templates live on those paths.
