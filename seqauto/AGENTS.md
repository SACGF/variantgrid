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
- `DragenTSO500CombinedVariantOutput` is one DRAGEN analysis of one pair (#1904) - TMB, MSI, GIS, the caller's tumour
  fraction and ploidy - keyed on (run name, pair) like `LibraryQC`, sharing its specimen claim
  (`seqauto/models/models_seqauto.py:SpecimenClaimMixin`). Its arms link to their `SequencingSample` and to the `Sample`
  whose `vcf_sample_name` is the arm's sample name, each filled later where it lands after the file
  (`patients/tasks/extraction_matching_tasks.py:link_combined_variant_outputs`, and `link_arm_sample` from
  `upload/vcf/vcf_import.py:link_samples_and_vcfs_to_sequencing`). The lab's MSI / TMB calls are properties over the
  `TSO500_*_CALL_BANDS` settings, not columns. The case report picks the row for its samples
  (`classification/report/case_report_context.py:case_combined_variant_output`).
API:
- A client-visible API change needs a name in `variantgrid/views_rest.py:API_FEATURES` (see
  claude/guides/operations.md#authentication-surface).
