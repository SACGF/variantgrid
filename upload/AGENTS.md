# upload — agent notes
Owns: FileUpload, UploadPipeline/UploadStep, UploadedVCF + Uploaded* satellites, the celery VCF import pipeline (preprocess
→ unknown-variant insert → parallel genotype insert → annotation wait → finish), ModifiedImportedVariant, VCFImportInfo.
Start with:
- models/models.py — FileUpload (the model; there is no UploadedFile), UploadPipeline, UploadStep, UploadedVCF,
  UploadedVCFPipelineMaxVariant, UploadedVCFPendingAnnotation, ModifiedImportedVariant. models/models_enums.py has
  VCFPipelineStage and UploadedFileTypes.
- file_type_icons.py — the icon each UploadedFileTypes value wears on the upload pages (FILE_TYPE_ICONS); a new
  type needs an entry there, test_file_type_icons checks. Drawn ones are file-icon-* in uicore's svg_icon_sprite.html.
- upload_processing.py — process_uploaded_file / process_upload_pipeline (retry) / process_vcf_file: entry points
  that pick a factory and fire the celery chain. A file's type comes off its contents on every entry point
  (import_task_factories/import_task_factory.py:get_import_task_factory_from_extension - the upload page and API through
  uploaded_file_type.py:get_uploaded_file_type, `import_vcf` and seqauto through
  upload_processing.py:get_vcf_file_type), so a SpliceGirl or gene-level CNV VCF takes its own path however it comes
  in. Only a VCF we wrote ourselves (ClinVar) passes an explicit file_type.
- import_task_factories/abstract_vcf_import_task_factory.py — AbstractVCFImportTaskFactory.create_import_task builds
  the step graph; import_task_factories.py has the per-type factories (GenotypeVCFImportFactory = plain VCF).
- tasks/vcf/import_vcf_step_task.py — ImportVCFStepTask (base of every VCF step) and schedule_pipeline_stage_steps.
- vcf/vcf_preprocess.py (bcftools pipe + split) · tasks/vcf/unknown_variants_task.py (variant creation) ·
  vcf/bulk_genotype_vcf_processor.py (genotype insert) · vcf/vcf_import.py (VCF from header, import loop)
Patterns here:
- Two task bases. Single-shot file types (BED, gene list, PED, patient records) subclass
  upload/tasks/import_task.py:ImportTask and return an item count from process_items(file_upload). Every VCF pipeline
  step subclasses upload/tasks/vcf/import_vcf_step_task.py:ImportVCFStepTask, takes an upload_step_id and returns
  items_processed (None skips the items_to_process == items_processed check). Register each task instance at module
  bottom with app.register_task() and list the module in variantgrid/settings/components/celery_settings.py:CELERY_IMPORTS,
  or the step dies with AttributeError when scheduled.
- Steps are rows. An UploadStep with pipeline_stage runs now; one with pipeline_stage_dependency is created inert and
  launched by upload/tasks/vcf/import_vcf_step_task.py:schedule_pipeline_stage_steps once the last step of that stage
  ends (ImportVCFStepTask.check_pipeline_stage). Stages: PRE_DATA_INSERTION → DATA_INSERTION → ANNOTATION_COMPLETE →
  FINISH (upload/models/models_enums.py:VCFPipelineStage). Add a step by adding its class to the factory's
  get_post_vcf_header_classes / get_post_data_insertion_classes / get_finish_task_classes
  (upload/import_task_factories/abstract_vcf_import_task_factory.py:AbstractVCFImportTaskFactory), or create the
  UploadStep inside a running step and call upload/models/models.py:UploadStep.launch_task.
- Preprocess is one shell pipe, not vt. upload/vcf/vcf_preprocess.py:_build_pipe_commands chains
  manage.py vcf_clean_and_filter | bcftools norm --multiallelics=- --old-rec-tag | vcf_clean_alts | split, retrying
  with bcftools sort when the unsorted marker file appears. Normalisation history rides in the INFO tag
  upload/models/models.py:ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG. Change record filtering in
  upload/management/commands/vcf_clean_and_filter.py:Command, not in the processors.
- The split stage needs GNU split, bash and the bgzip binary (htslib, apt package tabix) - a missing one only fails
  inside split's --filter at import time. upload/vcf/vcf_preprocess.py:get_split_vcf_command is the one place the
  command lives; manage.py deployment_check runs it on a tiny VCF ("VCF import split pipe") so a new box fails early.
- Variants are created in bulk by hash. upload/tasks/vcf/unknown_variants_task.py:SeparateUnknownVariantsTask runs per
  split file, batches coordinates through snpdb/variant_pk_lookup.py:VariantPKLookup and writes CSVs of unknowns;
  upload/tasks/vcf/unknown_variants_task.py:InsertUnknownVariantsTask re-checks and inserts them under a cache lock on
  the variant_id_single_worker queue. Create Variant/Locus rows only from that single worker, never from a parallel step.
- Genotypes go in via SQL COPY. upload/vcf/bulk_genotype_vcf_processor.py:BulkGenotypeVCFProcessor.process_entry
  accumulates rows, resolves variant ids per batch (settings.SQL_BATCH_INSERT_SIZE), splits common/rare by gnomAD AF,
  writes a CSV with upload/vcf/sql_copy_files.py:write_sql_copy_csv and spawns an
  upload/tasks/vcf/import_sql_copy_task.py:ImportCohortGenotypeSQLCopyTask step. Variants-only imports use
  upload/vcf/bulk_minimal_vcf_processor.py:BulkMinimalVCFProcessor (override batch_handle_variant_ids to act on ids).
- Every processor calls AbstractBulkVCFProcessor.set_max_variant; upload/vcf/vcf_import.py:import_vcf_file then upserts
  upload/models/models.py:UploadedVCFPipelineMaxVariant (one row per annotation pipeline type, only ever raised) —
  that is how the pipeline knows annotation is finished.
- Annotation wait is signal driven. upload/tasks/vcf/genotype_vcf_tasks.py:VCFCheckAnnotationTask creates
  UploadedVCFPendingAnnotation; upload/signals/signal_handlers.py:annotation_run_complete_signal_handler (wired in
  upload/apps.py:UploadConfig) re-checks is_fully_annotated after each annotation run and, when true, sends
  schedule_pipeline_stage_steps for ANNOTATION_COMPLETE. FINISH runs as a chain ending in pipeline_success_task;
  upload/tasks/vcf/genotype_vcf_tasks.py:ImportGenotypeVCFSuccessTask sets ImportStatus.SUCCESS and sends
  upload/signals/signals.py:vcf_import_success_signal (analysis/signals/signal_handlers.py:handle_vcf_import_success
  auto-creates analyses off it; connect further consumers in an AppConfig.ready).
Gotchas:
- A killed worker child (SIGTERM/SIGKILL, `revoke(terminate=True)`) never reaches `ImportVCFStepTask.run`'s except
  blocks; the master's `task_failure`/`task_revoked` receivers in `upload/tasks/vcf/import_vcf_step_task.py` fail the
  step and its pipeline instead. A whole-worker SIGKILL fires neither and still leaves the pipeline PROCESSING.
- A CNV VCF whose header declares a segment field naming a gene (settings.VCF_GENE_LEVEL_SEGMENT_FIELDS, DRAGEN
  TSO 500's `SEGID`) is claimed by import_task_factories/import_task_factories.py:GeneLevelCNVImportTaskFactory and
  rewritten onto the gene-level contig - the caller's segment is the panel's target window, not the event, so it is
  never stored as a Variant (tasks/import_gene_level_cnv_task.py, @see snpdb.gene_level_variants). A file naming a gene
  on partial calls (DragenExonCNV's `GENE=`) is not a segment field and keeps importing as coordinate SVs.
- A file type gated by a setting overrides `import_task_factories/import_task_factory.py:ImportTaskFactory.enabled`;
  a disabled factory is left out of `get_import_task_factories`, so it is neither picked for an upload nor listed by the
  capabilities endpoint. The five gene-level factories return `settings.VARIANT_GENE_LEVEL_ENABLED`, and with it off a
  `SEGID` CNV VCF or a SpliceGirl VCF imports as an ordinary VCF on its written coordinates.
- A SpliceGirl VCF (the TSO 500 RNA arm's SpliceVariants.vcf, recognised by its `##source=SpliceGirl` header line -
  the pipeline sends it as a plain VCF) is claimed by
  import_task_factories/import_task_factories.py:SpliceGirlImportTaskFactory and each `<DEL>` junction rewritten as a
  gene-level splice Variant, every record with its FILTER (tasks/import_splicegirl_vcf_task.py). The record names no
  gene: it comes from the SpliceEvent row at the breakpoints, else the one gene a transcript puts at the donor
  (genes/gene_splice.py:SpliceEventResolver.resolve_junction, which also takes a stored coordinate Variant). A junction in
  no gene, or ambiguously in two, is skipped and counted on the import page. The header keeps `##source`, so the
  `^SpliceGirl` VCFSourceSettings row still binds AD/DP as alt/ref depth.
- A TSO 500 pair's CombinedVariantOutput tsv is no variant source (#1903) - its `[Splice Variants]` are the VCF's PASS
  calls on EGFR, MET and AR (import_task_factories/import_task_factories.py:DragenTSO500CombinedVariantOutputImportTaskFactory,
  tasks/import_dragen_tso500_combined_variant_output_task.py). It writes a VCF of no records, whose one sample is the RNA
  arm, and the file declares no genome build, so one is declared at upload or comes off the
  `^DRAGEN TSO500 CombinedVariantOutput` VCFSourceSettings row.
- What that file is for is the pair's identity, written once the header step has made the Sample
  (get_post_vcf_header_classes - a VCF with no records skips every DATA_INSERTION-dependent step) by
  tasks/import_dragen_tso500_combined_variant_output_task.py:DragenTSO500CombinedVariantOutputInsertTask
  (tso500/dragen_combined_variant_output_records.py). `[Analysis Details]` names the Patient (the code
  `settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX` reads out of `Pair ID` - the lab writes that either as the whole pair
  sample name, whose leading sequencing sample ID changes when the patient is re-sequenced, or as the bare patient
  code, and the default regex reads both), the Specimen
  (the ten-digit accession inside each sample ID) and the two Extractions (its container suffix), created when absent;
  the DNA/RNA sample IDs are exact `Sample.vcf_sample_name` and `SequencingSample.sample_name`, which links both arms'
  samples to their extraction and the CVO's VCF to its sequencing run without seqauto's filename matching. The
  analysis itself - `[Analysis Details]` and the `[TMB]`, `[MSI]`, `[GIS]` scalars - is one
  seqauto/models/models_seqauto.py:DragenTSO500CombinedVariantOutput per (run, pair) (#1904), keyed on the upload's
  `sequencing_run` metadata as the MetricsOutput is (the file names its run 'NA'); a CVO sent without it takes the run
  whose current sheet names one of its sample IDs, and one no registered run names is not recorded. None of it fails
  the import - a chain that cannot be made parks the row's specimen claim and is a SimpleVCFImportInfo message.
  The record-less VCF exists only to make the RNA arm's Sample for the patient chain: once the chain no longer needs it
  (#1903 step 3), the factory becomes a single-shot ImportTask like the MetricsOutput's, with an
  `UploadedDragenTSO500CombinedVariantOutput` UploadData one-to-one with the row.
- The run's MetricsOutput.tsv is a separate, single-shot import (tasks/import_dragen_tso500_metrics_output_task.py,
  tso500/dragen_metrics_output_parser.py + _records.py): it has no variants and no coordinates, and recognises itself
  by a banner line ending 'Metrics Output' (with the module version, as the CVO's has it). One file covers a whole run
  - the pipeline sends Results/MetricsOutput_orig.tsv, the copy the lab's run wrapper leaves untouched - and a 2.6.2
  column is a *pair* (the CVO's Pair ID) carrying both arms, not a sample. `sequencing_run` upload metadata is
  required and is part of the key, since a pair column alone does not identify a pair across runs. It writes one
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
- Failure is one-way: UploadStep.error_exception → upload/models/models.py:UploadPipeline.error sets ERROR, marks the
  VCF/samples ImportStatus.ERROR, logs an Event and reports to Rollbar. Later steps see status != PROCESSING and mark
  themselves SKIPPED; BulkGenotypeVCFProcessor.check_pipeline_for_failures bails mid-file. Running steps are not killed.
- Retry (upload/upload_processing.py:process_upload_pipeline, upload/tasks/vcf/genotype_vcf_tasks.py:reload_vcf_task)
  deletes only steps with origin IMPORT_TASK_FACTORY and resets the rest (USER_ADDITION steps re-run as is), so the
  create-data-from-header tasks must stay idempotent.
- A VCF whose build cannot be resolved gets ImportStatus.REQUIRES_USER_INPUT and the pipeline TERMINATED_EARLY
  (upload/tasks/vcf/genotype_vcf_tasks.py:ImportCreateVCFModelForGenotypeVCFTask); declare genome_build/source as
  upload metadata (upload/upload_metadata.py:validate_upload_metadata).
- Sample ImportStatus only moves with its VCF, through snpdb/import_status.py:set_vcf_and_samples_import_status
  (success from ImportGenotypeVCFSuccessTask, error from UploadPipeline.error). A finish task list that closes the
  pipeline before the success task leaves the VCF Importing forever - the success step is SKIPPED, not run.
- Queues (celery_settings.py:CELERY_TASK_ROUTES, keyed by dotted class path): web_workers reads uploaded files,
  variant_id_single_worker inserts variants and zygosity counts, scheduling_single_worker runs
  schedule_pipeline_stage_steps; everything else lands on db_workers.
- The retry button (upload/views/views.py:upload_retry_import) only shows for the file's owner or a superuser, with the
  input file still on disk, UPLOAD_ENABLED, and the URL visible (Shariant hides it). `manage.py reload_imports` is the
  way in without it - `--uploaded_file_type`, `--status`, `--upload_pipeline_id`, `--dry-run`. An import whose failure
  was before UploadPipeline.objects.create has only an orphan FileUpload to show for itself and nothing to reload;
  classification imports get back in through `manage.py classification_rematch_stuck`, which re-derives the coordinates
  and runs new pipelines.
- settings.UPLOAD_ENABLED=False makes InsertUnknownVariantsTask raise; IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS
  wipes the processing dir on success, so inspect a failed pipeline's files before retrying it.
- upload/tasks/vcf/import_vcf_step_task.py:pipeline_success_task is the only thing that closes a VCF pipeline, and it is
  guarded on status == PROCESSING (FINISH can be scheduled twice). A FINISH task that sets SUCCESS itself therefore
  silently swallows the timings, the success Event and the file cleanup - #928 was exactly that. Leave the status alone
  in get_finish_task_classes tasks.
- library/django_utils/django_file_utils.py:get_import_processing_dir creates the directory; use
  import_processing_dir_path when you only want to name one, and remove_import_processing_dir to remove it.
  manage.py import_processing_cleanup --dry-run reports what is reclaimable under settings.IMPORT_PROCESSING_DIR.
Tests:
- Whole pipeline in-process: create a FileUpload, then process_uploaded_file(file_upload, run_async=False) under
  CELERY_TASK_ALWAYS_EAGER (library/django_utils/unittest_utils.py:URLTestCase sets it) —
  upload/tests/test_import_patient_records.py:TestPatientUploadImport is the pattern. The VCF pipeline needs bcftools
  and a reference fasta, so tests stop short of preprocess_vcf (upload/tests/vcf/test_vcf_preprocess.py mocks run_pipe).
- Processors without celery: upload/tests/vcf/test_vcf_processors.py:TestVCFProcessors builds an UploadStep + VCF from
  a cyvcf2 reader via create_vcf_from_vcf and feeds process_entry directly (needs Sequence rows for GATC).
- Test VCFs: upload/test_data/vcf (grch38_brca1.vcf, no_genotype.GRCh37.vcf, symbolic_alt/, detect_build_by_header/)
  and upload/test_data/multiallele.vcf. Real import on a dev box: manage.py import_vcf <vcf> --name X --user Y
  (upload/management/commands/import_vcf.py:Command).
- URL test: upload/tests/test_urls.py:Test (owner vs non-owner on view_uploaded_file / view_upload_pipeline and the
  modified-variants datatable). API upload and dedupe-by-sha256: upload/tests/test_api.py.
Deep reference: claude/research/upload.md · claude/maps/models.md#upload
