# seqauto — agent notes
Owns: SequencingRun, SampleSheet/SequencingSample, the sequencing file records (UnalignedReads, AlignmentFile, SingleSampleVCF,
JointCalledVCF), QC models, EnrichmentKit and gold coverage, the seqauto REST API, and the link from what upload
imports to those records (BackendVCF, `seqauto/upload_links.py`) with the DRAGEN TSO 500 run-level importers.
- seqauto sits above upload: upload never imports it (`.importlinter`). seqauto plugs in through upload's signals
  (`upload/signals/signals.py`, received in `seqauto/signals/upload_signal_handlers.py`) and import factories
  (`seqauto/import_task_factories.py`, listed in settings.IMPORT_TASK_FACTORY_IMPORTS). Its models name upload's by
  string and subclass `upload/upload_data.py:UploadData`, as upload's models import apps that import seqauto's.
- `LibraryQC` is per-library caller QC that arrives by upload (DRAGEN TSO500 MetricsOutput today), keyed on
  (run name, pair, category) with a nullable SequencingRun and a claimed Specimen - not a SeqAutoRecord, no path.
  Each row is one arm (`nucleic_acid`) and links that arm's `SequencingSample` through the sheet's `Pair_ID` /
  `Sample_Type` columns, which the pipeline posts as `SequencingSampleData`
  (`seqauto/models/models_seqauto.py:sequencing_sample_for_pair`); a sheet posted without them leaves the link null
  and `patients/tasks/extraction_matching_tasks.py:link_library_qc_to_sequencing_samples` fills it later. Reports
  read it by the case's pair, falling back to the specimen (`classification/report/case_report_context.py:case_library_qc`).
- `DragenTSO500CombinedVariantOutput` is one DRAGEN analysis of one pair (#1904) - TMB, MSI, GIS, the caller's tumour
  fraction and ploidy - keyed on (run name, pair) like `LibraryQC`, sharing its specimen claim
  (`seqauto/models/models_seqauto.py:SpecimenClaimMixin`). The lab's MSI / TMB calls are properties over the
  `TSO500_*_CALL_BANDS` settings, not columns. The case report picks the row for its samples
  (`classification/report/case_report_context.py:case_combined_variant_output`).
- An arm of a CVO links to its `SequencingSample` and to a `Sample` two ways: the `Sample` whose `vcf_sample_name` is
  the arm's sample name, or - for a single-sample VCF that names its one sample something of its own, as SpliceGirl's
  'SAMPLE' - the `Sample` seqauto matched to that arm's `SequencingSample` by filename. Both
  `seqauto/models/models_seqauto.py:DragenTSO500CombinedVariantOutput.link_samples` and `link_arm_sample` take either,
  each filled later where it lands after the file
  (`patients/tasks/extraction_matching_tasks.py:link_combined_variant_outputs`, and `link_arm_sample` from
  `seqauto/upload_links.py:link_samples_and_vcfs_to_sequencing`).
- Both rows are keyed on the same (run, pair), and that pair is the page: `seqauto/views.py:view_tso500_pair` shows the
  analysis, the QC of each arm with its metrics, and the arms' sequencing samples and Samples, which is where both
  models' `get_absolute_url` and preview point. It is readable by whoever can read the claimed specimen's patient, else
  by anyone who can read the run, else admin-only. The specimen page and the run's Data tab list a line per
  (run, pair, arm) through `seqauto/qc/library_qc_summary.py:summarise_library_qc` and link here for the metrics.

- A TSO 500 pair's CombinedVariantOutput tsv is no variant source (#1903) - its `[Splice Variants]` are the VCF's PASS
  calls on EGFR, MET and AR - and has no coordinates, so it is a single-shot import rather than a VCF pipeline
  (seqauto/import_task_factories.py:DragenTSO500CombinedVariantOutputImportTaskFactory,
  seqauto/tasks/import_dragen_tso500_combined_variant_output_task.py:ImportDragenTSO500CombinedVariantOutputTask, with an
  `UploadedDragenTSO500CombinedVariantOutput` UploadData one-to-one with the file). It writes no VCF and takes no
  `genome_build`; `sequencing_run` is its only metadata key.
- What that file is for is the pair's identity (seqauto/tso500/dragen_combined_variant_output_records.py).
  `[Analysis Details]` names the Patient (the code
  `settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX` reads out of `Pair ID` - the lab writes that either as the whole pair
  sample name, whose leading sequencing sample ID changes when the patient is re-sequenced, or as the bare patient
  code, and the default regex reads both), the Specimen
  (the ten-digit accession inside each sample ID) and the two Extractions (its container suffix), created when absent;
  the DNA/RNA sample IDs are exact `Sample.vcf_sample_name` and `SequencingSample.sample_name`, which links both arms'
  samples to their extraction without seqauto's filename matching. The
  analysis itself - `[Analysis Details]` and the `[TMB]`, `[MSI]`, `[GIS]` scalars - is one
  seqauto/models/models_seqauto.py:DragenTSO500CombinedVariantOutput per (run, pair) (#1904), keyed on the upload's
  `sequencing_run` metadata as the MetricsOutput is (the file names its run 'NA'); a CVO sent without it takes the run
  whose current sheet names one of its sample IDs. The run is half the row's key, so a pair no registered run names,
  or metadata naming a run seqauto has not registered, fails the import saying so, as does one naming no `Pair ID`. A chain that cannot be made does not fail it: the row
  is still written with its specimen claim parked saying why, which the pair's page shows
  (seqauto/views.py:view_tso500_pair).
- The run's MetricsOutput.tsv is a separate, single-shot import (seqauto/tasks/import_dragen_tso500_metrics_output_task.py,
  seqauto/tso500/dragen_metrics_output_parser.py + _records.py): it has no variants and no coordinates, and recognises itself
  by a banner line ending 'Metrics Output' (with the module version, as the CVO's has it). One file covers a whole run
  - the pipeline sends Results/MetricsOutput_orig.tsv, the copy the lab's run wrapper leaves untouched - and a 2.6.2
  column is a *pair* (the CVO's Pair ID) carrying both arms, not a sample. `sequencing_run` upload metadata is
  required, must name a run seqauto has registered, and is part of the key, since a pair column alone does not
  identify a pair across runs. It writes one
  seqauto/models/models_seqauto.py:LibraryQC per (run, pair, QC category), naming the arm each category is about and judging every
  metric of a known section by the file's own LSL/USL guideline - settings.TSO500_LIBRARY_QC_GUIDELINES overrides that
  per (section, metric) where a lab quotes its own number, and each metric records which it was judged by. The column
  names a Specimen by its trailing ten-digit accession and never creates one (that is the CVO's job). Each row links its
  arm's SequencingSample where the run's current sheet carries the TSO500 Pair_ID / Sample_Type columns as
  SequencingSampleData (seqauto/models/models_seqauto.py:sequencing_sample_for_pair). A bare-patient-code Pair ID has no
  accession, so the run supplies it: off the linked arm's name, or on a sheet without Pair_ID data off the SequencingSample
  whose name carries the code - the sheet's names carry the code and the accession together. A column nothing on the run is named for
  (a control, another assay on the flowcell) keeps its rows with the claim parked saying so, and an unresolvable
  accession parks like a VCF sample's, with reconcile_pending_extractions fired after the writes.
API:
- Reads are open to any logged-in user; writes need a superuser or a member of `SEQAUTO_API_WRITE_GROUP` (the
  pipeline's token user). A new write endpoint subclasses `seqauto/views_rest.py:SeqAutoModelViewSet` or
  `seqauto/views_rest.py:SeqAutoWriteAPIView`; a POST that only reads (batch lookups) stays a plain `APIView`.
  An upload links to sequencing data under the same rule (`seqauto/seqauto_permissions.py:can_write_seqauto`, through
  `seqauto/upload_links.py:check_may_link_to_sequencing_data`) - by `path` (a VCF to its SingleSampleVCF /
  JointCalledVCF, gene coverage to its QCGeneCoverage, `seqauto/upload_links.py:get_sequencing_record_for_path`)
  or by run name (TSO500 MetricsOutput / CombinedVariantOutput). A link refused fails the upload, checked before the
  import writes anything of its own.
- An endpoint for records of a run subclasses `seqauto/views_rest.py:SequencingRunModelViewSet` or
  `seqauto/views_rest.py:SequencingRunWriteAPIView` (calling `log_sequencing_runs` with what it saved), so its calls
  appear on the run page's API Log tab.
- The API's `alignment_files` (an old client's single `bam_file` is `[bam_file]`) are `AlignmentFile`s: BAM or CRAM,
  by `file_type` or a `.cram` extension, several per sequencing sample. `SingleSampleVCF` and `QC` link all of theirs
  (M2M, never one picked out) and carry `sequencing_sample` directly - query through that, not the alignment files.
  A VCF is keyed on (sequencing sample, variant caller): a re-analysis re-posting that pair replaces the path and
  the alignment files, whatever their order. A QC is keyed on (sequencing sample, VCF); a re-post naming alignment
  files replaces them (`seqauto/serializers/sequencing_serializers.py:SingleSampleVCFSerializer`,
  `seqauto/serializers/seqauto_qc_serializers.py:QCSerializer`). They reach IGV as `SampleFilePath` rows, copied at
  VCF import (`seqauto/upload_links.py:link_samples_and_vcfs_to_sequencing`) and, for a file posted after that, by
  `seqauto/models/models_seqauto.py:AlignmentFile.link_to_samples`.
- Serializers write through `library/django_utils/change_tracking.py:update_or_create_if_changed` (and
  `get_or_create_recorded` where an existing row is left alone), never `objects.update_or_create`, so a replayed post
  is a true no-op: `modified` means a field changed, and the API Log row shows the call as unchanged (#1413).
- A client-visible API change needs a name in `variantgrid/views_rest.py:API_FEATURES` (see
  claude/guides/operations.md#authentication-surface).
