# upload — agent notes
Owns: FileUpload, UploadPipeline/UploadStep, UploadedVCF + Uploaded* satellites, the celery VCF import pipeline (preprocess
→ unknown-variant insert → parallel genotype insert → annotation wait → finish), ModifiedImportedVariant, VCFImportInfo.
Start with:
- models/models.py — FileUpload (the model; there is no UploadedFile), UploadPipeline, UploadStep, UploadedVCF,
  UploadedVCFPipelineMaxVariant, UploadedVCFPendingAnnotation, ModifiedImportedVariant. models/models_enums.py has
  VCFPipelineStage and UploadedFileTypes.
- upload_processing.py — process_uploaded_file / process_upload_pipeline (retry) / process_vcf_file: entry points
  that pick a factory and fire the celery chain.
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
  get_post_data_insertion_classes / get_finish_task_classes
  (upload/import_task_factories/abstract_vcf_import_task_factory.py:AbstractVCFImportTaskFactory), or create the
  UploadStep inside a running step and call upload/models/models.py:UploadStep.launch_task.
- Preprocess is one shell pipe, not vt. upload/vcf/vcf_preprocess.py:_build_pipe_commands chains
  manage.py vcf_clean_and_filter | bcftools norm --multiallelics=- --old-rec-tag | vcf_clean_alts | split, retrying
  with bcftools sort when the unsorted marker file appears. Normalisation history rides in the INFO tag
  upload/models/models.py:ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG. Change record filtering in
  upload/management/commands/vcf_clean_and_filter.py:Command, not in the processors.
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
- Failure is one-way: UploadStep.error_exception → upload/models/models.py:UploadPipeline.error sets ERROR, marks the
  VCF/samples ImportStatus.ERROR, logs an Event and reports to Rollbar. Later steps see status != PROCESSING and mark
  themselves SKIPPED; BulkGenotypeVCFProcessor.check_pipeline_for_failures bails mid-file. Running steps are not killed.
- Retry (upload/upload_processing.py:process_upload_pipeline, upload/tasks/vcf/genotype_vcf_tasks.py:reload_vcf_task)
  deletes only steps with origin IMPORT_TASK_FACTORY and resets the rest (USER_ADDITION steps re-run as is), so the
  create-data-from-header tasks must stay idempotent.
- A VCF whose build cannot be resolved gets ImportStatus.REQUIRES_USER_INPUT and the pipeline TERMINATED_EARLY
  (upload/tasks/vcf/genotype_vcf_tasks.py:ImportCreateVCFModelForGenotypeVCFTask); declare genome_build/source as
  upload metadata (upload/upload_metadata.py:validate_upload_metadata).
- Queues (celery_settings.py:CELERY_TASK_ROUTES, keyed by dotted class path): web_workers reads uploaded files,
  variant_id_single_worker inserts variants and zygosity counts, scheduling_single_worker runs
  schedule_pipeline_stage_steps; everything else lands on db_workers.
- settings.UPLOAD_ENABLED=False makes InsertUnknownVariantsTask raise; IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS
  wipes the processing dir on success, so inspect a failed pipeline's files before retrying it.
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
