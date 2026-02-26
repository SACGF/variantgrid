# Upload App Research

## Purpose

The `upload` Django app implements the complete VCF/file import pipeline for VariantGrid. It handles VCF files, BED files, pedigree files, gene lists, and other data types through a multi-stage pipeline that includes variant normalization and database insertion.

## Core Models

### UploadedFile

Base model representing any file submitted to the system.

**Fields:**
- `path`: Filesystem path to the file
- `uploaded_file`: FileField using PrivateUploadStorage
- `sha256_hash`: Content hash for integrity verification
- `file_type`: UploadedFileTypes enum:
  - `V` = VCF
  - `B` = BED
  - `G` = GeneList
  - `O` = GeneCoverage
  - `P` = Pedigree
  - `M` = ManualVariantEntry
  - `R` = PatientRecords
  - `S` = Classifications
  - `A` = VariantTags
  - `Y` = VCF-VariantsOnly
  - `W` = WikiVariant
  - `w` = WikiGene
- `import_source`: Either `WEB_UPLOAD` or filesystem
- `user`: FK to user who uploaded the file

**Methods:**
- `can_view()`: Permission check
- `check_can_view()`: Permission check, raises exception if not allowed
- `get_filename()`: Returns the base filename
- `store_sha256_hash()`: Computes and stores the SHA-256 hash

### UploadPipeline

Tracks the overall import job for an uploaded file.

**Fields:**
- `status`: ProcessingStatus enum:
  - `CREATED`
  - `PROCESSING`
  - `SUCCESS`
  - `ERROR`
  - `TERMINATED_EARLY`
  - `TIMED_OUT`
- `uploaded_file`: OneToOne FK to UploadedFile
- `progress_status`: Human-readable status string
- `progress_percent`: Completion percentage
- `items_processed`: Count of processed items
- `processing_seconds_wall_time`: Wall clock time
- `processing_seconds_cpu_time`: CPU time consumed
- `celery_task`: Reference to the Celery task

**Methods:**
- `start()`: Transitions status to PROCESSING
- `success()`: Transitions status to SUCCESS
- `error()`: Transitions status to ERROR, sets import status on related data, creates an Event with LogLevel.ERROR, reports to Rollbar

### UploadStep

Represents an individual processing stage within a pipeline.

**Fields:**
- `upload_pipeline`: FK to UploadPipeline
- `parent_upload_step`: Self-referential FK (nullable) for nested steps
- `sort_order`: Execution ordering
- `status`: Processing status
- `origin`: Either `IMPORT_TASK_FACTORY` or `USER_ADDITION`
- `task_type`: Either `CELERY`, `SQL`, or `TOOL`
- `pipeline_stage`: VCFPipelineStage enum:
  - `PRE_DATA_INSERTION`
  - `DATA_INSERTION`
  - `ANNOTATION_COMPLETE`
  - `FINISH`
- `pipeline_stage_dependency`: Gates execution until the specified stage is reached
- `input_filename`: Input file for this step
- `output_filename`: Output file from this step
- `script`: Task class name
- `import_variant_table`: Temporary table used during import
- `celery_task`: Reference to the Celery task

**Methods:**
- `start()`: Begins the step
- `error_exception(e)`: Records exception, sets status to ERROR, propagates error up to parent and pipeline
- `mark_timed_out()`: Marks the step as timed out
- `launch_task(task_class)`: Dispatches the Celery task
- `close_sub_steps()`: Finalizes all child steps

### VCFImporter

Audit record of the import method and version used.

**Fields:**
- `name`: Importer name (e.g. "VariantGrid")
- `version`: Version string
- `vcf_parser`: Parser library name
- `vcf_parser_version`: Parser version
- `code_git_hash`: Git hash of the code at import time

Has a composite unique key across all fields.

### UploadedVCF

Links an uploaded file to the processed VCF record.

**Fields:**
- `uploaded_file`: OneToOne FK to UploadedFile
- `vcf`: OneToOne FK to snpdb.VCF (nullable)
- `upload_pipeline`: OneToOne FK to UploadPipeline (nullable)
- `max_variant`: FK to Variant (highest variant ID, used for annotation tracking)
- `vcf_importer`: FK to VCFImporter (nullable)

### UploadedVCFPendingAnnotation

Tracks annotation completion for an uploaded VCF.

**Fields:**
- `uploaded_vcf`: OneToOne FK to UploadedVCF
- `created`: Creation timestamp
- `finished`: Completion timestamp
- `schedule_pipeline_stage_steps_celery_task`: Reference to the scheduling task

**Methods:**
- `attempt_schedule_annotation_stage_steps()`: Checks annotation status and, when complete, schedules the next pipeline stage steps

### BackendVCF

Links a processed VCF to its filesystem source (used by SeqAuto).

**Fields:**
- `uploaded_vcf`: OneToOne FK to UploadedVCF
- `vcf_file`: OneToOne FK to seqauto.VCFFile (nullable)
- `combo_vcf`: OneToOne FK to seqauto.SampleSheetCombinedVCFFile (nullable)

**Properties:**
- `filesystem_vcf`: Returns the appropriate filesystem VCF record
- `sample_sheet`: Returns the associated sample sheet
- `variant_caller`: Returns the variant caller used

**Methods:**
- `get_samples_by_sequencing_sample()`: Maps sequencing samples to VCF samples

### VCFImportInfo and Subclasses

Import messages and warnings attached to a VCF import.

**SimpleVCFImportInfo:** Has `type`, `message_string`, `count` fields. `add_message_count()` performs batch count updates efficiently.

**ModifiedImportedVariants:** Tracks variant normalization events with a `count` field and links to individual `ModifiedImportedVariant` instances.

**ModifiedImportedVariant:** Records a single variant transformation.
- `variant`: FK to the final normalized variant
- `operation`: Either `NORMALIZATION` or `RMDUP`
- `old_variant`: Original variant string
- `old_variant_formatted`: Formatted original
- `old_multiallelic`: Original multiallelic representation

**VCFSkippedContigs:** Records chromosomes excluded during import.

### UploadedFile Subclasses

All are related via OneToOne FK to UploadedFile:

| Subclass | Links To |
|---|---|
| `UploadedGeneList` | GeneList |
| `UploadedBed` | GenomicIntervalsCollection |
| `UploadedPedFile` | PedFile |
| `UploadedPatientRecords` | PatientRecords |
| `UploadedGeneCoverage` | GeneCoverageCollection + Sample |
| `UploadedManualVariantEntryCollection` | ManualVariantEntryCollection |
| `UploadedClassificationImport` | ClassificationImport |
| `UploadedLiftover` | LiftoverRun |
| `UploadedWikiCollection` | (wiki data) |
| `UploadedClinVarVersion` | ClinVarVersion |
| `UploadedVariantTags` | VariantTagsImport |

### UploadSettings

Per-user preferences for the upload interface.

**Fields:**
- `time_filter_method`: Either `DAYS` or `RECORDS`
- `time_filter_value`: Numeric value for the filter
- `file_types`: M2M to file type choices via `UploadSettingsFileType`

## VCF Import Pipeline End-to-End

The pipeline proceeds through the following stages:

**Stage 1: File Receipt**
User uploads VCF via web form or API. `handle_file_upload()` is called, which creates an `UploadedFile` record and calls `process_uploaded_file()`.

**Stage 2: Pipeline Initialization**
`UploadPipeline` is created. `ImportTaskFactory.create_import_task()` instantiates the chain of `UploadStep` records that will drive the process.

**Stage 3: Preprocessing (PRE_DATA_INSERTION)**
`PreprocessVCFTask` runs two normalization tools against the raw VCF:
- `vt decompose`: Splits multi-allelic records into biallelic records
- `vt normalize`: Left-aligns and normalizes indels

The preprocessed VCF is then split by chromosome into individual per-chromosome files.

**Stage 4: Parallel Chromosome Processing**
`ScheduleMultiFileOutputTasksTask` launches parallel child `UploadStep` tasks, one per chromosome file.

**Stage 5: Bulk Insertion (DATA_INSERTION)**
`BulkGenotypeVCFProcessor` (or a no-genotype variant for variant-only VCFs) reads each chromosome VCF using `cyvcf2`. It builds SQL COPY-formatted files and bulk-inserts variants and genotypes into the database.

**Stage 6: Annotation Monitoring**
An `UploadedVCFPendingAnnotation` record is created. The annotation engine is notified of the new variants (via `max_variant`). The `attempt_schedule_annotation_stage_steps()` method is polled until annotation finishes.

**Stage 7: Post-Annotation Steps (ANNOTATION_COMPLETE)**
Steps with `pipeline_stage_dependency=ANNOTATION_COMPLETE` are triggered. These include linking alleles, ClinGen allele registration, liftover creation, and similar post-annotation operations.

**Stage 8: Finalization (FINISH)**
`UploadPipelineFinishedTask` sets the pipeline status to SUCCESS and performs optional cleanup of temporary files and intermediate data.

## Variant Normalization Tracking

When `vt decompose` and `vt normalize` transform variants, the tools embed the original coordinates in the VCF INFO field as `OLD_VARIANT`. The pipeline captures these in `ModifiedImportedVariant` records.

For multiallelic splits, `vt` writes `OLD_MULTIALLELIC` in the format `1:100:A:C,G|1`, where `|1` indicates which allele index this record came from. The pipeline parses these to create two `ModifiedImportedVariant` records, one per allele.

This enables a reverse lookup: given an original (unnormalized) coordinate from the source VCF, find the final normalized `Variant` record in the database.

## Celery Tasks

| Task | Description |
|---|---|
| `PreprocessVCFTask` | Runs vt decompose/normalize, splits by chromosome |
| `CheckStartAnnotationTask` | Triggers annotation engine for new variants |
| `ScheduleMultiFileOutputTasksTask` | Launches parallel per-chromosome tasks |
| `ProcessVCFSetMaxVariantTask` | Records the highest variant ID after insertion |
| `ProcessVCFLinkAllelesSetMaxVariantTask` | Links alleles and sets max variant |
| `ProcessVCFClinGenAlleleTask` | Registers variants with ClinGen allele registry |
| `LiftoverCreateVCFTask` | Creates liftover VCF for alternate genome builds |
| `ImportCreateUploadedVCFTask` | Creates the UploadedVCF record |
| `UploadPipelineFinishedTask` | Finalizes pipeline status |
| `import_gene_list_task` | Imports gene list from uploaded file |
| `import_bedfile_task` | Imports BED file regions |
| `import_ped_task` | Imports pedigree file |
| `import_patient_records_task` | Imports patient record data |
| `import_gene_coverage_task` | Imports gene coverage data |

## Error Propagation

Errors propagate upward through the pipeline hierarchy:

1. An exception in a task calls `UploadStep.error_exception(e)`
2. `error_exception()` sets the step `status=ERROR` and calls `upload_pipeline.error(message)`
3. `upload_pipeline.error()` sets the pipeline `status=ERROR`
4. The pipeline error handler sets `ImportStatus` on any related data objects (e.g., the snpdb.VCF record)
5. An `Event` is created with `LogLevel.ERROR` for the audit log
6. The error is reported to Rollbar for external monitoring

## Views

**Web views:**
- Upload form (file selection and submission)
- `jfu_upload`: jQuery File Upload POST endpoint
- `jfu_upload_delete`: Delete an uploaded file
- `view_uploaded_file`: Details page for a single uploaded file
- `view_upload_pipeline`: Pipeline status page with AJAX progress polling
- `view_upload_pipeline_warnings_and_errors`: Warnings and error details page
- `upload_retry_import`: Retry a failed import
- `accept_vcf_import_info_tag`: Accept/dismiss an import warning

**REST API:**
- `APIFileUploadView`: Programmatic file upload endpoint
