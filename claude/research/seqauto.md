# SeqAuto App Research

## Purpose

The `seqauto` Django app implements sequencing automation for VariantGrid. It scans the filesystem for raw genomics data (FASTQs, BAMs, VCFs), creates corresponding database records, orchestrates QC analysis, and bridges raw sequencing files with the VariantGrid analysis pipeline.

## Core Models

### SeqAutoRun

Represents one complete filesystem scan execution.

**Fields:**
- `status`: SeqAutoRunStatus enum:
  - `CREATED`
  - `SCANNING_FILES`
  - `CREATE_MODELS`
  - `SCRIPTS_AND_JOBS`
  - `FINISHED`
  - `ERROR`
- `task_id`: Celery task identifier
- `scan_start`: Timestamp when scanning began
- `create_models_start`: Timestamp when model creation began
- `scripts_and_jobs_start`: Timestamp when job script generation began
- `finish_date`: Completion timestamp
- `error_exception`: Stored exception if status is ERROR
- `job_launch_script_filename`: Path to the generated job launch script
- `fake_data`: FK to fake data configuration (for testing)

**Methods:**
- `get_status()`: Returns current status
- `get_scan_resources_dir()`: Returns the directory for intermediate scan outputs
- `get_job_scripts_dir()`: Returns the directory for generated job scripts
- `remove_scan_resources_dir()`: Cleans up intermediate scan data
- `get_last_success_datetime()` (static): Returns the datetime of the most recent successful run

### SeqAutoRecord

Polymorphic base model for all filesystem-discovered records. Uses `InheritanceManager` for polymorphic queries.

**Fields:**
- `path`: Absolute filesystem path
- `sequencing_run`: FK to SequencingRun (nullable)
- `file_last_modified`: Float timestamp of file modification time
- `hash`: Content or path hash
- `is_valid`: Validation result boolean
- `data_state`: DataState enum:
  - `COMPLETE`
  - `RUNNING`
  - `DELETED`
  - `SKIPPED`
  - `NON_EXISTENT`
  - `ERROR`

**Methods:**
- `validate()`: Validates the record, sets `is_valid`
- `add_messages()`: Attaches SeqAutoMessage records
- `_close_messages_with_code()`: Closes open messages matching a code

### SeqAutoMessage

Stores validation results and informational messages for scan records.

**Fields:**
- `seqauto_run`: FK to SeqAutoRun (at least one of seqauto_run/record must be set)
- `record`: FK to SeqAutoRecord (at least one of seqauto_run/record must be set)
- `severity`: LogLevel enum: `ERROR` / `WARNING` / `INFO` / `DEBUG`
- `code`: Machine-readable message code
- `message`: Human-readable message text
- `open`: Bool indicating whether the condition is still present

### SequencingRun

The primary entity representing one sequencing run. The primary key is the run name.

**Fields:**
- `name` (PK): Unique run identifier string
- `date`: Run date
- `sequencer`: FK to Sequencer
- `enrichment_kit`: FK to EnrichmentKit (nullable)
- `experiment`: FK to Experiment (nullable)
- `gold_standard`: Bool for reference-quality runs
- `bad`: Bool flagging problematic runs
- `hidden`: Bool hiding the run from normal views
- `legacy`: Bool indicating the run should not be auto-updated
- `has_basecalls`: Bool indicating basecall data is present
- `has_interop`: Bool indicating InterOp data is present
- `fake_data`: Bool for synthetic test data

**Methods:**
- `get_params()`: Returns a dict for string substitution with keys: `sequencing_run`, `flowcell_id`, `sequencing_run_dir`, `enrichment_kit`, `experiment`. Used when building job script command lines.
- `get_current_sample_sheet()`: Returns the active SampleSheet
- `get_original_illumina_sequencing_run()`: Returns the original Illumina run record if linked

### SampleSheet

CSV configuration file defining the samples in a run.

**Fields:**
- `sequencing_run`: FK to SequencingRun
- `path`: Filesystem path to the SampleSheet.csv file
- `hash`: Content hash for change detection

**Methods:**
- `get_sequencing_samples_by_name()`: Returns a dict mapping sample name to SequencingSample
- `set_as_current_sample_sheet()`: Makes this sample sheet active and signals downstream change detection
- `get_params()`: Returns substitution parameters

### SequencingRunCurrentSampleSheet

A join table enforcing one active SampleSheet per SequencingRun. Each side is a OneToOne relationship.

### SequencingSample

Represents a single sample row within a SampleSheet.csv.

**Fields:**
- `sample_sheet`: FK to SampleSheet
- `sample_id`: Sample identifier from the sheet
- `sample_name`: Sample name
- `sample_project`: Project name
- `sample_number`: Numeric position in the sheet
- `lane`: Lane number (nullable for non-lane-split runs)
- `barcode`: Index barcode sequence
- `enrichment_kit`: FK to EnrichmentKit (nullable, overrides run-level kit)
- `is_control`: Bool for control samples
- `failed`: Bool marking failed samples
- `automatically_process`: Bool controlling whether SeqAuto will auto-process this sample

**Methods:**
- `get_params()`: Builds a substitution dict combining run-level and sample-level parameters

### SampleFromSequencingSample

Links a VariantGrid `snpdb.Sample` to its originating `SequencingSample`.

**Fields:**
- `sample`: OneToOne FK to snpdb.Sample
- `sequencing_sample`: FK to SequencingSample

**Properties:**
- `sequencing_run`: Convenience access to the parent SequencingRun

### IlluminaFlowcellQC

Flowcell-level quality metrics loaded from Illumina InterOp data.

**Fields:**
- `sample_sheet`: OneToOne FK to SampleSheet
- `mean_cluster_density`: Mean cluster density (K/mm^2)
- `mean_pf_cluster_density`: Mean passing-filter cluster density
- `total_clusters`: Total cluster count
- `total_pf_clusters`: Passing-filter cluster count
- `percentage_of_clusters_pf`: Percent of clusters passing filter
- `aligned_to_phix`: Percent alignment to PhiX control

Loaded via `illuminate_report.load_from_file()`.

### ReadQ30

Per-read Q30 percentage metrics.

**Fields:**
- `illumina_flowcell_qc`: FK to IlluminaFlowcellQC
- `read`: Read identifier (R1 / R2 / I1 / I2)
- `percent`: Q30 percentage value
- `is_index`: Bool distinguishing index reads from sequence reads

### Fastq

Represents a single FASTQ file.

**Fields:**
- `sequencing_sample`: FK to SequencingSample
- `name`: FASTQ filename
- `read`: Read direction (R1 or R2)

**Static methods:**
- `get_pair_paths_from_sequencing_sample()`: Returns the R1/R2 path pair for a sample

### FastQC

FastQC report metrics for a FASTQ file.

**Fields:**
- `fastq`: OneToOne FK to Fastq
- `total_sequences`: Total sequence count
- `filtered_sequences`: Filtered sequence count
- `gc`: GC content percentage

Parsed from `fastqc_data.txt` output files.

### UnalignedReads

Groups an R1/R2 FASTQ pair for a sample.

**Fields:**
- `sequencing_sample`: FK to SequencingSample
- `fastq_r1`: FK to Fastq (R1)
- `fastq_r2`: FK to Fastq (R2)

**Properties:**
- `data_state`: Combined state derived from both R1 and R2 Fastq states

### BamFile

Represents a BAM alignment file.

**Fields:**
- `unaligned_reads`: FK to UnalignedReads
- `name`: BAM filename
- `aligner`: FK to Aligner

### Flagstats

Samtools flagstat mapping statistics for a BAM file.

**Fields:**
- `bam_file`: OneToOne FK to BamFile
- `total`: Total read count
- `read1`: Read 1 count
- `read2`: Read 2 count
- `mapped`: Mapped read count
- `properly_paired`: Properly paired count

**Properties:**
- `mapped_percent`: Percentage of reads mapped
- `properly_paired_percent`: Percentage of reads properly paired

### VCFFile

Per-sample VCF produced from a BAM file.

**Fields:**
- `bam_file`: FK to BamFile
- `variant_caller`: FK to VariantCaller

`load_from_file()` triggers VCF import into VariantGrid if `SEQAUTO_IMPORT_VCF=True` in settings.

### SampleSheetCombinedVCFFile

Joint-called multi-sample VCF covering an entire sample sheet.

**Fields:**
- `sample_sheet`: FK to SampleSheet
- `variant_caller`: FK to VariantCaller

**Methods:**
- `needs_to_be_linked()`: Returns whether this VCF still needs to be linked to VariantGrid samples
- `get_paths_from_sample_sheet()`: Derives expected VCF paths from the sample sheet

### QC

Groups all QC records for a single sample.

**Fields:**
- `bam_file`: FK to BamFile
- `vcf_file`: FK to VCFFile

**Properties:**
- `genome_build`: Genome build derived from the enrichment kit

### QCExecSummary

Comprehensive QC metrics parsed from an `exec_summary.txt` file. Contains 30+ numeric fields including:

**Key fields:**
- `deduplicated_reads`: Read count after deduplication
- `indels_dbsnp_percent`: Percent of indels in dbSNP
- `mean_coverage_across_genes`: Mean depth across gene targets
- `percent_10x_goi`: Percent of gene-of-interest bases at 10x coverage
- `percent_20x_kit`: Percent of kit target bases at 20x coverage
- `percent_100x_goi`: Percent of gene-of-interest bases at 100x coverage
- `ts_to_tv_ratio`: Transition/transversion ratio
- `percent_softclip`: Percent of soft-clipped bases
- `uniformity_of_coverage`: Uniformity metric
- `percent_read_enrichment`: On-target read percentage
- `percent_duplication`: PCR duplication rate

### ExecSummaryReferenceRange

Expected value ranges for QC metrics using DecimalRangeField types.

**Fields:**
- `percent_20x` range: Expected range for 20x coverage percentage
- `percent_10x` range: Expected range for 10x coverage percentage
- `mean_coverage` range: Expected mean coverage range
- `min_mean_coverage_across_kit`: Minimum acceptable mean coverage
- `min_percent_20x_kit`: Minimum acceptable 20x coverage percentage

### QCGeneCoverage

Links QC to a gene coverage collection.

**Fields:**
- `qc`: OneToOne FK to QC
- `gene_coverage_collection`: OneToOne FK to GeneCoverageCollection

### QCGeneList

Gene-of-interest list used for QC calculations.

**Fields:**
- `qc`: FK to QC
- `custom_text_gene_list`: OneToOne FK to CustomTextGeneList (nullable)
- `sample_gene_list`: FK to SampleGeneList (nullable)

**Methods:**
- `create_gene_list()`: Creates the underlying gene list
- `link_samples_if_exist()`: Connects to existing samples if available
- `load_from_file()`: Parses a gene list file and populates the record

## Scan Flow

The scan process is initiated by `scan_run_jobs()` via Celery and proceeds through three phases:

**Phase 1: Scan (SCANNING_FILES)**
External scripts are invoked to discover filesystem data. These scripts write their output as JSON or TSV files into `scan_resources_dir`. Discovery covers:
- Sample sheets
- FASTQ files
- BAM files
- Flagstat files
- VCF files
- QC metric files (exec_summary, FastQC, InterOp)

**Phase 2: Create Models (CREATE_MODELS)**
The scan outputs are parsed to create or update `SeqAutoRecord` subclass instances. Each record's `validate()` method is called to check for issues. `SeqAutoMessage` records are created for any warnings or errors found.

**Phase 3: Scripts and Jobs (SCRIPTS_AND_JOBS)**
PBS or Slurm job scripts are generated for any samples that need processing (alignment, variant calling, QC). The launch script is executed to submit jobs to the cluster scheduler.

## QC Metric Hierarchy

```
SequencingRun
├── IlluminaFlowcellQC (via SampleSheet)
│   └── ReadQ30 (per read: R1/R2/I1/I2)
└── SampleSheet
    └── SequencingSample
        ├── UnalignedReads
        │   ├── Fastq (R1) + FastQC
        │   └── Fastq (R2) + FastQC
        └── BamFile (via UnalignedReads)
            ├── Flagstats
            └── QC
                ├── QCExecSummary
                ├── ExecSummaryReferenceRange
                ├── QCGeneCoverage → GeneCoverageCollection
                └── QCGeneList
```

## Scheduled Tasks and Management Commands

**Celery beat task:**
- `scan_run_jobs()`: Runs the full scan cycle at a configurable frequency

**Management commands:**
- `import_sequencing_info`: Imports sequencing run information from external sources
- `scan_run_jobs`: Manual trigger for the scan cycle
- `job_script_complete`: Marks a job script as completed
- `job_script_submitted`: Records job script submission to the scheduler
- `reload_illumina_flowcell_qc`: Re-parses and reloads InterOp QC data
- `reload_qc_gene_coverage`: Re-parses and reloads gene coverage data
- `set_gold_standard_runs`: Flags specified runs as gold standard references
