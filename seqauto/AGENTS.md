# seqauto — agent notes
Owns: SequencingRun, SampleSheet/SequencingSample, the sequencing file records (UnalignedReads, BamFile, SingleSampleVCF,
JointCalledVCF), QC models, EnrichmentKit and gold coverage, and the seqauto REST API.
- `LibraryQC` is per-library caller QC that arrives by upload (DRAGEN TSO500 MetricsOutput today), keyed on
  (run name, pair, category) with a nullable SequencingRun and a claimed Specimen - not a SeqAutoRecord, no path.
  Each row is one arm (`nucleic_acid`) and links that arm's `SequencingSample` through the sheet's `Pair_ID` /
  `Sample_Type` columns, which the pipeline posts as `SequencingSampleData`
  (`seqauto/models/models_seqauto.py:sequencing_sample_for_pair`); a sheet posted without them leaves the link null
  and `patients/tasks/extraction_matching_tasks.py:link_library_qc_to_sequencing_samples` fills it later. Reports
  read it by specimen (`classification/report/case_report_context.py:specimen_library_qc`).
API:
- A client-visible API change needs a name in `variantgrid/views_rest.py:API_FEATURES` (see
  claude/guides/operations.md#authentication-surface).
