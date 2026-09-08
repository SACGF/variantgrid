# upload — research notes

Verified against 7c4408c62 on 2026-09-06

The upload app turns a file a user or API client sent into database rows, asynchronously, with a record of every step
taken. For most file types that is one celery task; for a VCF it is a multi-stage pipeline of `UploadStep` rows that
normalises the file with bcftools, creates whatever Variants do not exist yet through a single serialised worker,
bulk-loads genotypes with SQL COPY in parallel, waits for VEP annotation to catch up, and only then declares the VCF
imported. The rules for touching it are in `upload/CLAUDE.md`; this is the longer story of why those rules exist.
Models, URLs, tasks, commands and signals are enumerated in the generated maps ([models](../maps/models.md#upload),
[urls](../maps/urls.md#upload), [tasks](../maps/tasks.md#upload), [commands](../maps/commands.md),
[signals](../maps/signals.md#first-party)); this doc does not restate them.

## Flows

### From POST to a pipeline

The web page is a FilePond widget: `upload/views/views_json.py:upload_file` is its `process` endpoint and refuses
outright when `variantgrid/settings/components/default_settings.py:UPLOAD_ENABLED` is off. Both it and the token API
`upload/views/views_rest.py:APIFileUploadView.post` funnel into `upload/views/views_json.py:handle_file_upload`, which
creates the `upload/models/models.py:FileUpload` with `import_source=WEB_UPLOAD` (even for the API - the `path` a client
sends is a hint about *its* filesystem, not ours), works out the file type from the extension through
`upload/uploaded_file_type.py:get_uploaded_file_type`, validates any upload metadata while the client is still connected
(`upload/upload_metadata.py:validate_upload_metadata` - a bad key or an unresolvable build is a 400 now, not a failed
import three stages later), stores the sha256 (`FileUpload.store_sha256_hash`) and calls
`upload/upload_processing.py:process_uploaded_file`. The API de-duplicates by that hash before saving anything:
`APIFileUploadView._get_existing_file_upload` returns the newest successfully-processed upload with the same content the
user can view, so a client re-sending a file gets the old `file_upload_id` back unless it passes `force`. The two
companion endpoints, `upload/views/views_rest.py:APIUploadStatusView` and `APIAnnotatedDownloadView`, accept either the id
or the hash, because the hash is stable across servers.

`upload/upload_processing.py:process_uploaded_file` creates the one `upload/models/models.py:UploadPipeline` a FileUpload
may have and hands it to `process_upload_pipeline`, which is also the retry path: it deletes every step whose origin is
`IMPORT_TASK_FACTORY`, resets the rest to CREATED, then asks `get_upload_processing_task` for the factory whose
`get_uploaded_file_type` matches. Factories are the `upload/import_task_factories/import_task_factory.py:ImportTaskFactory`
subclasses collected by `get_import_task_factories`; each names its extensions, the `UploadData` classes it produces and
the metadata keys it accepts, and `create_import_task` returns a celery signature that `process_upload_pipeline` applies
(synchronously when `run_async=False`, which is how tests and `manage.py import_vcf` run the whole thing in-process).
Command-line and SeqAuto imports enter one level down at `upload/upload_processing.py:process_vcf_file`, which builds the
FileUpload from a server path with `import_source` set to where it really came from.

Single-shot file types (BED, gene list, PED, patient records, gene coverage, variant tags, wiki, analysis) return one task
built on `upload/tasks/import_task.py:ImportTask`: `ImportTask.run` starts the pipeline, calls `process_items(file_upload)`,
and records success with the item count or error with the traceback. Everything VCF-shaped - plain VCF, variants-only,
gene-level, manual variant entry, liftover, ClinVar, DRAGEN TSO500 fusions - subclasses
`upload/import_task_factories/abstract_vcf_import_task_factory.py:AbstractVCFImportTaskFactory` instead, and that is
where the multi-stage machinery starts.

### The VCF step graph

`AbstractVCFImportTaskFactory.create_import_task` builds two things. The first is a celery chain of the steps that can
run straight away: `upload/tasks/vcf/import_vcf_step_task.py:pipeline_start_task`, an optional `get_pre_vcf_task` (the
fusion loader uses it to write a VCF from a CSV first), the "Create Data from VCF Header" step and the "Preprocess VCF"
step. The second is a set of inert `upload/models/models.py:UploadStep` rows whose `pipeline_stage_dependency` names the
stage they wait for: `CheckStartAnnotationTask` and `ScheduleMultiFileOutputTasksTask` on PRE_DATA_INSERTION, the
factory's `get_post_data_insertion_classes` on DATA_INSERTION, and `get_finish_task_classes` on FINISH. For a genotype
VCF (`upload/import_task_factories/import_task_factories.py:GenotypeVCFImportFactory`) the DATA_INSERTION set is
`VCFCheckAnnotationTask`, `UpdateVariantZygosityCountsTask`, `SampleLocusCountsTask` and, when `settings.SOMALIER` is
enabled, `SomalierVCFTask`; FINISH is whatever `settings.FINISH_IMPORT_VCF_STEP_TASKS_CLASSES` names (a deployment hook,
empty by default) followed by `ImportGenotypeVCFSuccessTask`. A pipeline with no VCF file at all gets a
`DoNothingVCFTask` step in DATA_INSERTION purely so the dependency chain has something to complete.

Every step task is an `upload/tasks/vcf/import_vcf_step_task.py:ImportVCFStepTask`, a celery `Task` subclass registered as
an instance with `app.register_task` at the bottom of its module so the UploadStep can store its dotted class path in
`script` and `schedule_pipeline_stage_steps` can `import_class` it later. `ImportVCFStepTask.run` loads the step,
refuses to run twice (a step already holding another `celery_task` id raises), marks itself SKIPPED without doing
anything when the pipeline is no longer PROCESSING, otherwise calls `process_items(upload_step)` and compares the count
returned against `items_to_process` when both are set. Whatever happened, it stamps `end_date`, closes sub-steps and then
- deliberately after its own row is closed - calls `ImportVCFStepTask.check_pipeline_stage`.

`check_pipeline_stage` is the whole scheduling algorithm: if this step was the last unfinished step *of its stage* and
there are CREATED steps waiting on that stage, it fires `schedule_pipeline_stage_steps(pipeline, stage)`; separately,
if nothing but FINISH-dependent steps remain, it fires the FINISH stage. The race where two parallel steps both finish
"last" is resolved by running `upload/tasks/vcf/import_vcf_step_task.py:schedule_pipeline_stage_steps` on the
`scheduling_single_worker` queue and having it set `start_date` on each waiting step as it launches it, so a second
invocation finds nothing to launch. FINISH steps run as a chain terminating in `pipeline_success_task` rather than a
chord, because the chord version produced endless `celery.chord_unlock` retries (#175).

### Header first: the VCF model

`upload/tasks/vcf/genotype_vcf_tasks.py:ImportCreateVCFModelForGenotypeVCFTask` reads only the header with cyvcf2 and
calls `upload/vcf/vcf_import.py:create_vcf_from_vcf`, which reuses or creates the `UploadedVCF`, creates the
`snpdb` VCF with the raw header saved immediately (so a crash later still leaves something to inspect), resolves the
build, configures INFO/FORMAT/FILTER models and samples from the header (`configure_vcf_from_header`), links SeqAuto
data when the deployment has it (`create_backend_vcf_links`, `link_samples_and_vcfs_to_sequencing`) and assigns any
extraction metadata (`assign_sample_extractions`). `upload/vcf/vcf_import.py:resolve_genome_build` orders the evidence:
contig lengths in the header (`vcf_detect_genome_build_from_header`), then the build the submitter declared as upload
metadata, then the build configured for the VCF's `##source` in `VCFSourceSettings`, then the only annotated build on a
single-build server. Detected and declared disagreeing is a `GenomeBuildMismatchException`, not a guess. Nothing
resolving leaves `genome_build` null, and the task sets the VCF to `REQUIRES_USER_INPUT` and the pipeline to
`TERMINATED_EARLY` so every later step skips itself. The same task creates the cohort and its
`CohortGenotypeCollection` (`create_cohort_genotype_collection_from_vcf`) and adds the ANNOTATION_COMPLETE-dependent
"Calculate VCF Stats" step plus one mutational-signature step per somatic-only sample - steps a running step creates
are picked up by the scheduler like any other.

### Preprocess: one shell pipe

`upload/tasks/vcf/import_vcf_tasks.py:PreprocessAndAnnotateVCFTask` (the genotype factory's choice; plain
`PreprocessVCFTask` for the others, and `GeneLevelPreprocessVCFTask` splits only, since gene-level variants have no
reference base for bcftools to check) calls `upload/vcf/vcf_preprocess.py:preprocess_vcf`. It writes a cleaned header
first, then `_build_pipe_commands` assembles `zcat | manage.py vcf_clean_and_filter | bcftools norm | bcftools view
--no-header | manage.py vcf_clean_alts | split`, and `run_pipe` runs it under `set -o pipefail` (#3813) with the split
filter re-attaching the header and bgzipping each chunk of `settings.VCF_IMPORT_FILE_SPLIT_ROWS` lines. The stages that
carry a versioned tool become child `UploadStep` rows via `create_sub_step`, which is where `ToolVersion` and the
per-stage stdout/stderr live. `upload/management/commands/vcf_clean_and_filter.py:Command` is the gatekeeper bcftools
needs: it renames contigs to the reference's names, drops non-standard contigs and out-of-range positions (counted into
`VCFSkippedContigs`), drops records matching `settings.VCF_IMPORT_SKIP_RECORD_REGEX` and reference spans, moves FILTER
values the header never declared into INFO so bcftools does not abort on them (#1711), and watches sort order. An
unsorted file writes a marker and fails the pipe; `preprocess_vcf` catches that, `_reset_for_retry` clears the partial
outputs, and the pipe is rebuilt with `bcftools sort` spliced in (#127). `bcftools norm` splits multi-allelics
(`--multiallelics=-`), left-aligns, replaces `N` reference bases from the fasta (`--check-ref=s`, #888) and records what
it changed in the INFO tag `upload/models/models.py:ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG`; `--rm-dup` is
deliberately absent because it discards that tag (#985), so duplicates are handled later in Python.
`upload/management/commands/vcf_clean_alts.py` runs after the split so it can reject or convert one bad alt without
losing the rest of the record - symbolic alts, oversize explicit alts and `<DUP:TANDEM>` (#1247) are decided here.

`upload/vcf/vcf_preprocess.py:schedule_split_file_steps` then creates one PRE_DATA_INSERTION "Separate Unknown Variants
Task" step per chunk and launches it at once, optionally an "Annotate gnomAD AF" step
(`upload/tasks/vcf/unknown_variants_task.py:AnnotateImportedVCFTask`, `bcftools annotate` from the per-build file in
`settings.VCF_IMPORT_COMMON_FILTERS` into `INFO/VG_GNOMAD_AF`), and an `UploadStepMultiFileOutput` row per chunk. Those
rows are the hand-off: `ScheduleMultiFileOutputTasksTask` later fans the DATA_INSERTION step out over them. The last
chunk's size is unknown, so their `items_to_process` is left null and the count check is skipped.

### Unknown variants on the single worker

`upload/tasks/vcf/unknown_variants_task.py:SeparateUnknownVariantsTask` reads a chunk and feeds every record's
`VariantCoordinate` to `snpdb/variant_pk_lookup.py:VariantPKLookup`, which hashes coordinates and resolves them against
the Variant table in batches of `settings.SQL_BATCH_INSERT_SIZE`. Whatever is unknown after `VariantPKLookup.batch_check`
is written to a CSV by `handle_unknown_variants`, which also creates and launches a "Create Unknown Loci and Variants"
step. That step runs `InsertUnknownVariantsTask`, routed to `variant_id_single_worker` and additionally guarded by the
cache lock `insert-unknown-variants-lock`, and it re-checks every coordinate before inserting
(`InsertUnknownVariantsTask.process_items_in_lock`, `VariantPKLookup.batch_check` with `insert_unknown=True`) because a
parallel chunk or another pipeline may have created it since the CSV was written. The invariant this buys is one Variant
row per coordinate without a database unique constraint on the hot path - which is why nothing outside that worker may
create Locus or Variant rows. `InsertUnknownVariantsTask` also raises when `UPLOAD_ENABLED` is false: the setting exists
to stop variant creation when a deployment is out of disk, and it is enforced here rather than only at the form.

When the last PRE_DATA_INSERTION step ends, `upload/tasks/vcf/import_vcf_tasks.py:CheckStartAnnotationTask` pokes the
annotation scheduler if the pipeline inserted anything (`UploadStep.pipeline_inserted_unknown_variants` looks for that
step name), and `ScheduleMultiFileOutputTasksTask` creates the "Process VCF File" steps via
`ImportVCFStepTask._schedule_steps`. A VCF with no records left after filtering has no chunks; `_handle_no_vcf_records`
records a warning and marks every DATA_INSERTION- and ANNOTATION_COMPLETE-dependent step SKIPPED so the pipeline still
reaches FINISH (#1116).

### Genotype insert with SQL COPY

`upload/tasks/vcf/genotype_vcf_tasks.py:ProcessGenotypeVCFDataTask` runs once per chunk on `web_workers`, building a
processor with `upload/vcf/vcf_import.py:genotype_vcf_processor_factory` (`BulkGenotypeVCFProcessor` when the VCF has
samples, `upload/vcf/bulk_no_genotype_vcf_processor.py` otherwise) and looping `upload/vcf/vcf_import.py:import_vcf_file`.
`upload/vcf/bulk_genotype_vcf_processor.py:BulkGenotypeVCFProcessor.process_entry` converts each record to its internal
form (`AbstractBulkVCFProcessor.get_ref_alt_svlen`, which also forces negative SVLEN on `<DEL>` because Manta writes
positive ones, #1245), and holds records until the locus changes (`finished_locus`) because a decomposed multi-allelic
spreads one sample's allele depths over several rows and the reference depth has to be summed across them.
`BulkGenotypeVCFProcessor.batch_process_check` resolves the accumulated hashes to variant ids in one query, raises the
pipeline's max variant (`AbstractBulkVCFProcessor.set_max_variant`), and calls `process_cohort_genotypes`, which drops
duplicate variant ids (a `RMDUP` ModifiedImportedVariant each, since bcftools could not), splits rows into the common and
rare partitions by the gnomAD AF the preprocess stamped in - except for ids in `_get_uncommon_variant_ids`, variants with
a classification that must stay findable however common they are - writes a CSV with
`upload/vcf/sql_copy_files.py:write_sql_copy_csv` and creates a child `ImportCohortGenotypeSQLCopyTask` step
(`upload/tasks/vcf/import_sql_copy_task.py`) to COPY it in. Modified-variant rows go the same way through
`ImportModifiedImportedVariantSQLCopyTask`. `BulkGenotypeVCFProcessor.check_pipeline_for_failures` re-reads the pipeline
status once per row at most and throws `PipelineFailedJobTerminateEarlyException` so a dead pipeline stops burning a
worker. `VCFImporter` records importer version, cyvcf2 version and git hash on the UploadedVCF; the version list in
`BulkGenotypeVCFProcessor.get_vcf_importer_version` is the changelog of import semantics, and bumping it is how a
reload of an old VCF is later told apart from a fresh one.

`import_vcf_file` ends with `upload/vcf/vcf_import.py:update_uploaded_vcf_max_variant`, an upsert of one
`UploadedVCFPipelineMaxVariant` row per annotation pipeline type that only ever raises the id. Chunks run in parallel,
so the maximum has to be monotone across them; and it is per pipeline type (standard VEP, SV, AnnotSV) because those are
annotated by separate runs (#720, #1656).

### Waiting for annotation, then finishing

`upload/tasks/vcf/genotype_vcf_tasks.py:VCFCheckAnnotationTask` runs when DATA_INSERTION completes. It replaces any
`UploadedVCFPendingAnnotation` from an earlier attempt and calls
`upload/models/models.py:UploadedVCFPendingAnnotation.attempt_schedule_annotation_stage_steps`, which asks
`UploadedVCF.is_fully_annotated` whether, for every pipeline type with a max-variant row, the lowest unannotated variant
id in the active VAV is past that maximum. If so it marks itself finished and sends
`schedule_pipeline_stage_steps(pipeline, ANNOTATION_COMPLETE)` by task name (a string, to dodge a circular import).
If not, nothing polls: `upload/signals/signal_handlers.py:annotation_run_complete_signal_handler`, connected in
`upload/apps.py:UploadConfig.ready`, re-runs the check for every unfinished pending row of that build - narrowed to
VCFs that have variants of the pipeline type that just completed - each time an annotation run finishes. A VCF is
therefore only "imported" once its variants are annotated, and the annotation side never needs to know about pipelines.

The ANNOTATION_COMPLETE steps compute the VCF's stats (`CalculateVCFStatsTask`, warning when VEP skipped variants, #1409)
and mutational signatures. FINISH runs `ImportGenotypeVCFSuccessTask`, which sets `ImportStatus.SUCCESS` on the VCF and
samples, clears data-archive fields if this was a restore (#1536), sends `backend_vcf_import_success_signal` for SeqAuto
and `upload/signals/signals.py:vcf_import_success_signal` for everyone else, and writes the summary message; then
`pipeline_success_task` totals wall and CPU seconds from the step rows and `UploadPipeline.success` deletes the
processing directory when `IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS` is set.

### Failure and retry

Failure is one-way and top-down: `UploadStep.error_exception` stores the traceback and calls
`upload/models/models.py:UploadPipeline.error`, which sets ERROR, pushes `ImportStatus.ERROR` onto the VCF and samples
(`_set_related_data_import_status`), logs an Event and reports to Rollbar. Nothing kills running steps; they notice on
their next `check_pipeline_for_failures` or when they start and find the pipeline not PROCESSING. A step that hangs can be
closed by hand with `UploadStep.mark_timed_out`, which fails the pipeline the same way.

Retry (`upload/views/views.py:upload_retry_import` → `upload/uploaded_file_type.py:retry_upload_pipeline`) wipes the
processing directory and, for anything that loads a VCF, queues `upload/tasks/vcf/genotype_vcf_tasks.py:reload_vcf_task`
on the single worker: it subtracts the VCF's zygosity counts, deletes the genotype data (`VCF.delete_internal_data`)
and calls `process_upload_pipeline` on the same pipeline row. Because the UploadedVCF and VCF survive, the header task
reuses them (`create_vcf_from_vcf` reloads the existing UploadedVCF; `ImportCreateVCFModelForGenotypeVCFTask` re-detects
formats on an existing VCF) and everything downstream must be safe to run twice - which is the reason for the
`update_or_create` calls, the pending-annotation delete and the max-variant upsert. Steps a person added by hand
(`origin=USER_ADDITION`) are reset and re-run rather than deleted. `upload/tests/test_retry_import.py` is the contract.

## Why it is shaped this way

Steps are rows, not celery state, because an import outlives any one worker process: the page at
`upload/views/views.py:view_upload_pipeline` and `APIUploadStatusView` read progress straight from `UploadStep`, a
crashed worker leaves a row that says where it died, and retry can reason about what to delete. The stage/dependency
split (`upload/models/models_enums.py:VCFPipelineStage`) is what lets steps be created during the run - by preprocess
for each chunk, by the header task per sample - and still be scheduled correctly, since the scheduler only ever asks
"is anything in this stage still open" rather than following a fixed list.

Variants are found by hash, not by unique index, because `snpdb_variant` is the biggest hot table in every deployment
and a unique constraint over (locus, alt, svlen) with the insert rate of a whole-genome VCF was the bottleneck;
`VariantPKLookup` computes the same hash the SQL side can, so lookups are index scans and inserts are COPY. The cost is
the serialised `InsertUnknownVariantsTask`, which is why chunking keeps the unknown-variant CSVs small and the genotype
phase (which only needs existing ids) is the parallel one.

Genotypes go through CSV and COPY rather than the ORM because a cohort genotype row packs every sample's calls into
arrays (`snpdb/models/models_cohort.py:CohortGenotype`), and the common/rare split exists so that analyses filtering on
population frequency never scan the bulk of a genome's calls - the split is decided at import from the gnomAD AF
annotated in preprocess, which is why that annotation step sits in the pipe rather than waiting for VEP.

Normalisation history is kept as `ModifiedImportedVariant` because users search for the coordinate that was in their
file, not the left-aligned one: `ModifiedImportedVariant.get_variant_for_unnormalized_variant` and
`get_other_loci_variants_by_multiallelic` answer "where did my variant go", and `ModifiedImportedVariants.get_for_pipeline`
hangs the collection off the normalize sub-step so the tool version is recorded with it. Rows written before the bcftools
move still hold vt-format `OLD_VARIANT` strings, which `_vt_split_old_variant` parses, so the model reads both.

## History

The pipeline began as `vt decompose | vt normalize` and moved to bcftools over 2023-25: multi-allelic splitting was
moved ahead of normalisation to handle gVCFs (#1236), `N` reference bases are replaced from the fasta (#888), and
`--rm-dup` was dropped once it was found to strip the old-record tag (#985); the last vt references left with the
unsorted-VCF work (#127, 2026-08), which also added the `bcftools sort` retry. Standard contigs only and the
clean-and-filter command arrived with the T2T build (#814, #1082). SVLEN handling was hardened for Manta (#1245) and
DRAGEN (#1268), and symbolic alt conversion added (#1247), each bumping `VCF_IMPORTER_VERSION`.

Annotation waiting was rewritten per pipeline type in 2026-07 (#1656) after the AnnotSV lane was split from the SV VEP run
(#720): a single `max_variant` meant a stalled SV run blocked VCFs whose short variants were long done. The token API
for upload, status and annotated download landed the same month (#1640), with sha256 de-duplication (md5 was retired
in #1092) and the upload-metadata mechanism for declaring build and source (#1711) and, for SA Pathology, sample
extractions (#1707). The web widget moved from jQuery File Upload to FilePond (#1566, 2026-05); `UploadedFile` was
renamed `FileUpload` in 2026-08 because it collided with Django's own class. Gene fusions from DRAGEN TSO500
AllFusions.csv became a VCF-shaped pipeline (#1506), which is why gene-level variants have their own preprocess.
Guards for steps deleted mid-flight (#1593), the re-runnable retry, pipefail in the pipe (#3813) and the private-repo
security hardening of paths (#3825) are all 2026.

## Traps

The pipeline reports SUCCESS from `pipeline_success_task`, but the VCF's `import_status` is set by
`ImportGenotypeVCFSuccessTask` in the same FINISH chain; a deployment class listed in
`FINISH_IMPORT_VCF_STEP_TASKS_CLASSES` runs before it and can fail the pipeline after all data is in.

`is_fully_annotated` compares against the *active* VAV at check time (`AnnotationVersion.latest`). Creating a new VAV while
imports are pending makes every pending VCF wait for the full re-annotation; the `annotation_run_complete_signal` will
eventually release them but the upload page shows them stuck for hours.

`schedule_pipeline_stage_steps` catches `AttributeError` and reports "you probably didn't register your Celery class":
that is the symptom of a task module missing from `CELERY_IMPORTS` or a class not passed through `app.register_task`,
and it shows up only at the stage transition, not at import time.

`process_upload_pipeline` refuses a FileUpload that already has a pipeline; a second `process_uploaded_file` on the same
file is a bug, and the API's hash de-duplication is what keeps clients from triggering it.

A chunk's `SeparateUnknownVariantsTask` and its `InsertUnknownVariantsTask` are both PRE_DATA_INSERTION, so the stage
cannot complete while any insert is queued behind the lock - a long queue on `variant_id_single_worker` looks like a
stalled preprocess. `vg status` shows the queue depth.

`IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS` deletes the split chunks and per-stage logs; a failed pipeline keeps
them under `UploadPipeline.get_pipeline_processing_dir`, and retry deletes them before re-running, so read the
sub-step `output_text` on the pipeline page first.

`create_backend_vcf_links` raises when a SeqAuto deployment receives a VCF whose `path` matches no registered
`SingleSampleVCF` or `JointCalledVCF`; on non-SeqAuto deployments the same `path` is ignored entirely, so a client
sending `path` behaves differently per server.
