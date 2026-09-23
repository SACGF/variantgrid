# CombinedVariantOutput: single-shot import, and a page for the TSO500 pair

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-23
Status: draft

Follow-up to #1904 (the `DragenTSO500CombinedVariantOutput` row, landed ddf07d123) and #1903 step 3 (the CVO
stopped being a variant source). Two things are left over from those:

1. **The import still writes a VCF.** `upload/import_task_factories/import_task_factories.py:DragenTSO500CombinedVariantOutputImportTaskFactory`
   is an `AbstractVCFImportTaskFactory` whose pre-VCF step writes a header-only VCF naming the RNA arm as its one
   sample, so the normal header step makes a Sample for the patient chain to hang off
   (`upload/tasks/import_dragen_tso500_combined_variant_output_task.py`). Nothing needs that Sample any more:
   the patient / specimen / extraction chain is made by `resolve_pair` without it,
   `link_samples_to_extractions` links every Sample carrying the arm's `vcf_sample_name`, and the real arm VCFs get
   their seqauto link and extraction by filename (cbf5a94d2). What the tester sees is a 0-variant VCF on the data
   page and the run page, an empty Sample that can be put in an analysis, a `genome_build` the upload has to
   declare for a file with no coordinates, and the CVO row's `rna_sample` pointing at that empty Sample.
2. **Neither DRAGEN record has a page.** The analysis row is listed on the specimen page and the run page's Data
   tab (`seqauto/templates/seqauto/combined_variant_output_table.html`) with half its columns; the MetricsOutput's
   `LibraryQC` rows are listed on the specimen page only, with their metrics flattened into one cell, and the run
   page does not show them at all. The MetricsOutput import itself is already the single-shot shape
   (`upload/import_task_factories/import_task_factories.py:DragenTSO500MetricsOutputImportTaskFactory`), so it
   needs only the page side of this plan.

## Data

No new column on either row. `DragenTSO500CombinedVariantOutput` and `LibraryQC` gain `PreviewModelMixin` and
`get_absolute_url` (behaviour, no migration).

```python
class UploadedDragenTSO500CombinedVariantOutput(UploadData):   # upload/models/models_uploaded_files.py
    """ What makes the upload 'processed' (the API's sha256 de-duplication keys on it), as
        UploadedDragenTSO500MetricsOutput does for MetricsOutput.tsv. The analysis row is keyed on (run, pair)
        and replaced by a re-analysis, so this points at nothing: get_data() is the row whose file_upload is ours """
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
```

Data migration (new file on top of the frozen `snpdb/migrations/0264_vcf_source_settings_combined_variant_output_reads.py`
and `snpdb/migrations/0266_vcf_source_settings_combined_variant_output_genome_build.py`): delete the
`VCFSourceSettings` row `^DRAGEN TSO500 CombinedVariantOutput`. Nothing writes a VCF with that source any more.

Existing deployments keep the empty VCFs already imported. A `ManualOperation` in the same migration runs a one-off
command that soft-deletes VCFs whose `source` starts with `DRAGEN TSO500 CombinedVariantOutput` through
`snpdb/tasks/soft_delete_tasks.py:soft_delete_vcfs`, so they leave the data page the way a user-deleted VCF does
and the CVO rows' `rna_sample` (SET_NULL) is refilled by the arm linking below.

## Import

`DragenTSO500CombinedVariantOutputImportTaskFactory` becomes an `ImportTaskFactory` shaped exactly on
`DragenTSO500MetricsOutputImportTaskFactory`: `get_metadata_keys` is `{SEQUENCING_RUN}` only, `get_data_classes` is
the new UploadData, `create_import_task` returns `ImportDragenTSO500CombinedVariantOutputTask.si(upload_pipeline.pk)`.

`upload/tasks/import_dragen_tso500_combined_variant_output_task.py` keeps one task,
`ImportDragenTSO500CombinedVariantOutputTask(ImportTask)`, whose `process_items(file_upload)` is today's insert step
minus the VCF:

- parse, `parse_pair_identifiers`, `resolve_pair`, `link_samples_to_extractions` as now;
- the run: `sequencing_run` metadata, else `sequencing_run_for_sample_ids`. Neither is an import failure today, it
  is a message saying nothing was recorded. With no UploadStep to hang a `SimpleVCFImportInfo` on, a file that
  would record nothing raises `ValueError` and fails the pipeline saying so, as the MetricsOutput task does;
- a chain that cannot be made is not a failure: `write_combined_variant_output` already parks the claim with the
  error on the row (`specimen_match_error`), which the page below shows;
- `get_or_create` the UploadData, `reconcile_pending_extractions.delay()`, return 1.

Gone with the VCF: `DragenTSO500CombinedVariantOutputCreateVCFTask`, `_write_sample_vcf`, `_sample_name`,
`ALT_READS_FORMAT` / `REF_READS_FORMAT`, `source_from_analysis_details` and `SOURCE_PREFIX`,
`upload/tso500/dragen_combined_variant_output_records.py:link_to_sequencing_run`, and the gene-level contig imports.
The module docstring loses its VCF paragraphs. `seqauto/views.py` names the CVO as an example of an uploaded VCF
linked after the fact (line 99): the MetricsOutput or the SpliceGirl VCF is the example now.

### Arm samples

`seqauto/models/models_seqauto.py:DragenTSO500CombinedVariantOutput.link_arm_sample` matches a landing arm VCF's
Sample on `vcf_sample_name`. The SpliceGirl VCF's sample column is `SAMPLE`, so the RNA arm never matches by name,
and the empty VCF was the only Sample that did. Both linkers (`link_arm_sample`, and the reconcile pass in
`patients/tasks/extraction_matching_tasks.py:link_combined_variant_outputs`) accept an arm Sample either way: its
`vcf_sample_name` is the arm's sample ID, or its `SampleFromSequencingSample` is the arm's linked
`SequencingSample`. The case report already walks both routes (`classification/report/case_report_context.py:case_combined_variant_output`).

## Page

Both records are keyed on (run name, pair ID) - the analysis once, the QC six times - and that pair is the thing a
tester wants to see whole: the analysis and its QC together, with the two arms it was run on. So the page is the
pair on the run, not either row.

- `seqauto/urls.py`: `view_tso500_pair/<str:sequencing_run_name>/<str:pair_id>` → `views.view_tso500_pair`.
  `get_absolute_url` on both rows builds it. Permission: readable by whoever can read the claimed specimen's
  patient, else the run (the rule the table's links already imply); a pair with neither is admin-only.
- A new template, seqauto/templates/seqauto/view_tso500_pair.html:
  - the arms: a two-row table (DNA, RNA) linking each arm's `SequencingSample` and `Sample` where linked and
    saying which is still waiting; the specimen (linked, or the parked reference with its match error); the run;
  - the analysis (`DragenTSO500CombinedVariantOutput`, absent where only the MetricsOutput has landed): every field
    grouped as the file has them (`[Analysis Details]`, `[TMB]`, `[MSI]`, `[GIS]`), the computed TMB and MSI
    calls with the band they came off (`tmb_call` / `msi_call`), and the original upload (`file_upload`);
  - the QC (`LibraryQC`, absent where only the CVO has landed): one table per arm, a row per category with
    pass / fail and completion, expanding to the metrics with value, unit, LSL / USL, the guideline source and
    each metric's own verdict - the `metrics` JSON the specimen page currently prints in one cell - and the
    MetricsOutput upload it came off.
- `PreviewModelMixin` on both models: categories "TSO500 analysis" and "Library QC", `preview_enabled` on
  `settings.VARIANT_GENE_LEVEL_ENABLED` as the factories are, so search and hover previews work. A `LibraryQC`
  preview is its pair's page too.
- `combined_variant_output_table.html`: the pair cell links to the page. The specimen page's Library QC table
  (`patients/templates/patients/view_specimen.html`) collapses to one row per (run, pair, arm) with the categories'
  verdicts, linking to the page for the metrics.
- `seqauto/templates/seqauto/view_sequencing_run.html` Data tab: the same QC summary for the run's pairs under the
  analyses table, so a run's QC is visible from the run.
- `patients/templates/patients/view_patient_specimens.html`: under the specimen links, the analyses table for
  the patient's specimens (`show_specimen=True`), so a tester lands on the pair from the patient without opening
  each specimen. The patient page's "related data" stays sample- and cohort-centric. The data page gets no tab: it
  lists variant-bearing files, and each row's provenance is the upload page (`UploadData.get_data` gives the
  uploaded file's "view data" link the row).

## Tests

`upload/tests/test_import_dragen_tso500_combined_variant_output.py`: the VCF-shaped tests
(`test_names_the_rna_sample_and_the_module_and_has_no_records`, `test_the_rna_arm_links_the_vcf_to_its_sequencing_run`,
`test_a_sample_sheet_without_the_arm_leaves_no_link_rows`, `test_the_seqauto_rows_are_written_for_the_rna_arm`) go;
the end-to-end pipeline tests run the single-shot task and assert the row, the UploadData and no VCF. New:
a SpliceGirl-shaped Sample (column `SAMPLE`, linked to the RNA SequencingSample) fills `rna_sample`; a pair no run
names fails the pipeline with the message; the pair page renders for a pair with both records, with the analysis
only, with the QC only, and for a parked one (`vg page` on each). `scripts/vg tests --explain` names the rest.

## Docs

The CVO bullets in `upload/AGENTS.md` (the "record-less VCF" paragraph and its promised follow-up) and
`seqauto/AGENTS.md` (both the `LibraryQC` and the analysis bullets) are rewritten to describe the single-shot import,
the two arm-sample routes and the pair page.
