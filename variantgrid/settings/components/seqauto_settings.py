import os
from variantgrid.settings.components.default_settings import BASE_DIR, MANAGE_COMMAND, PRIVATE_DATA_ROOT

SEQAUTO_ENABLED = False
SEQAUTO_ALLOW_NON_STAFF_MANUAL_RUN = True
SEQAUTO_DIR = os.path.join(BASE_DIR, "seqauto")
SEQAUTO_SCRIPTS_DIR = os.path.join(SEQAUTO_DIR, "scripts", "test")
SEQAUTO_SCAN_RESOURCES_DIR = os.path.join(PRIVATE_DATA_ROOT, 'scan_resources')
SEQAUTO_USER = 'seqauto'
SEQAUTO_GROUP = None

TAU_DIR = os.path.expanduser('~/localwork/tau')
SEQAUTO_VIRTUALENV_RUNNER = None
SEQAUTO_SCRIPT_PARAMS = {"virtualenv_runner": SEQAUTO_VIRTUALENV_RUNNER or '',
                         "tau_pipeline_dir": os.path.join(TAU_DIR, "scripts"),
                         "tau_scripts_dir": os.path.join(TAU_DIR, "archive"),
                         "pythonpath": TAU_DIR}

SEQAUTO_FLOWCELL_SCRIPT = 'find_flowcells.sh'
SEQAUTO_FASTQ_SCRIPT = 'find_fastqs.sh'
SEQAUTO_FASTQC_SCRIPT = 'find_fastqc.sh'
SEQAUTO_ILLUMINATE_QC = 'find_illuminate_qc.sh'
SEQAUTO_BAM_SCRIPT = 'find_bams.sh'
SEQAUTO_FLAGSTATS_SCRIPT = None  # Don't do this anymore 'find_flagstats.sh'
SEQAUTO_VCF_SCRIPT = 'find_vcfs.sh'
SEQAUTO_QC_SCRIPT = 'find_qc.sh'

SEQAUTO_SKIP_FLOWCELLS_FILE = None
SEQAUTO_SKIP_FLOWCELLS_PATTERNS = []
SEQAUTO_SKIP_INDIVIDUAL_FLOWCELL_FILE = ".variantgrid_skip_flowcell"

SEQAUTO_CONTROL_SAMPLE_REGEX = None
SEQAUTO_ALIGNED_BASE_DIR = os.path.join(SEQAUTO_DIR, "test_data", "clinical_hg38")
SEQAUTO_SCRATCH_BASE_DIR = "/tmp"
SEQAUTO_GOLD_BASE_DIR = None

SEQAUTO_ALIGNED_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_BASE_DIR, "%(enrichment_kit)s", "%(lab_run_id)s_%(original_sequencing_run)s")
SEQAUTO_BAM_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_DIR_PATTERN, "1_BAM")
SEQAUTO_VCF_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_DIR_PATTERN, "2_variants")
SEQAUTO_QC_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_DIR_PATTERN, "4_QC")

SEQAUTO_MISEQ_ALIGNED_PATTERN = "%(sample_name_underscores)s_S%(sample_number)s"
SEQAUTO_HISEQ_ALIGNED_PATTERN = "%(sample_id)s"

SEQAUTO_BAM_PATTERN = "%(sample_name)s.hg38.bam"
SEQAUTO_VCF_PATTERN = "gatk_per_sample/%(sample_name)s.hg38.vcf.gz"
SEQAUTO_COMBINED_VCF_PATTERN = "%(lab_run_id)s_%(original_sequencing_run)s.gatk.hg38.vcf.gz"

SEQAUTO_QC_EXEC_SUMMARY_PATTERN = "exec_stats/%(sample_name)s_qc_summary.txt"
SEQAUTO_QC_EXEC_SUMMARY_TSV_PATTERN = "exec_stats/%(sample_name)s_stats.txt"
SEQAUTO_QC_GENE_LIST_PATTERN = "exec_stats/%(sample_name)s_QC.txt"
SEQAUTO_QC_GENE_COVERAGE_PATTERN = "bam_stats/samples/%(sample_name)s.per_gene_coverage.tsv.gz"

SEQAUTO_JOB_SCRIPTS_BASE_DIR = os.path.join(PRIVATE_DATA_ROOT, "job_scripts")
SEQAUTO_JOB_SCRIPTS_OUT_DIR = os.path.join(SEQAUTO_JOB_SCRIPTS_BASE_DIR, "out")
SEQAUTO_JOB_SCRIPTS_DIR = os.path.join(SEQAUTO_JOB_SCRIPTS_BASE_DIR, "scripts")

# PBS related
SEQAUTO_USE_PBS = False
SEQAUTO_PBS_QUEUE = None
SEQAUTO_PBS_EMAIL = None

TAU_PYTHON_COMMAND = "python"

SEQAUTO_LAUNCH_SCRIPT_PATHS_LIST = None
manage_command_str = " ".join(MANAGE_COMMAND)
SEQAUTO_JOB_SUBMITTED = f"{manage_command_str} job_script_submitted"
SEQAUTO_JOB_COMPLETE = f"{manage_command_str} job_script_complete"
SEQAUTO_SAMPLE_SHEET_COMMAND = "%(tau_scripts_dir)s/basecall_sample_sheet.sh %(sequencing_run_dir)s %(sample_sheet)s"
SEQAUTO_ILLUMINA_FLOWCELL_QC_COMMAND = "%(virtualenv_runner)s %(tau_scripts_dir)s/illuminate_qc.sh %(sequencing_run_dir)s"
SEQAUTO_FASTQC_COMMAND = "%(tau_scripts_dir)s/fastqc.sh %(fastq)s"
SEQAUTO_BAM_COMMAND = TAU_PYTHON_COMMAND + ' %(tau_pipeline_dir)s/BWA_2_hap.py -r "%(sequencing_run)s" -p "%(enrichment_kit)s" -s "%(full_sample_name)s"'
SEQAUTO_FLAGSTATS_COMMAND = "%(tau_scripts_dir)s/bam_index_flagstats.sh %(bam)s"
SEQAUTO_VCF_COMMAND = None  # Done as part of BAM script
SEQAUTO_COMBINED_VCF_COMMAND = TAU_PYTHON_COMMAND + ' %(tau_pipeline_dir)s/hap_2_VCF.py -r "%(sequencing_run)s" -p "%(enrichment_kit)s"'
SEQAUTO_QC_COMMAND = TAU_PYTHON_COMMAND + ' %(tau_pipeline_dir)s/single_sample_QC.py -r "%(sequencing_run)s" -p "%(enrichment_kit)s" -s "%(full_sample_name)s" -c -u'
SEQAUTO_MIGRATE_COMMAND = TAU_PYTHON_COMMAND + ' %(tau_scripts_dir)s/delete_temp_and_archive.py -r "%(sequencing_run)s" -p "%(enrichment_kit)s"'

SEQAUTO_IMPORT_VCF = False
SEQAUTO_IMPORT_COMBO_VCF = True
SEQAUTO_MIN_COVERAGE = 20
SEQAUTO_LOAD_GENE_COVERAGE = True
SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS = []

SEQAUTO_COVERAGE_ENRICHMENT_KITS = []
