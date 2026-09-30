import os

from variantgrid.settings.components.default_settings import BASE_DIR

SEQAUTO_ENABLED = False
SEQAUTO_DIR = os.path.join(BASE_DIR, "seqauto")
SEQAUTO_USER = 'seqauto'
SEQAUTO_GROUP = None
# Writes to the seqauto REST API (/seqauto/api/) need a superuser or a member of this group - add the
# user whose DRF token the sequencing pipeline posts with. Reads are open to any logged-in user
SEQAUTO_API_WRITE_GROUP = "seqauto_api_write"

# Path patterns derive a default path for QC / gene coverage / GOI records when the API omits one
SEQAUTO_ALIGNED_BASE_DIR = os.path.join(SEQAUTO_DIR, "test_data", "clinical_hg38")
SEQAUTO_ALIGNED_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_BASE_DIR, "%(enrichment_kit)s", "%(sequencing_run)s")
SEQAUTO_GOI_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_DIR_PATTERN, "0_goi")
SEQAUTO_BAM_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_DIR_PATTERN, "1_BAM")
SEQAUTO_QC_DIR_PATTERN = os.path.join(SEQAUTO_ALIGNED_DIR_PATTERN, "4_QC")
SEQAUTO_RUN_PARAMETERS_SUB_DIR = "4_QC/sequencing_stats"  # Subdir of SequencingRun

SEQAUTO_MISEQ_ALIGNED_PATTERN = "%(sample_name_underscores)s_S%(sample_number)s"
SEQAUTO_HISEQ_ALIGNED_PATTERN = "%(sample_id)s"

SEQAUTO_BAM_PATTERN = "%(sample_name)s.hg38.bam"
SEQAUTO_GOI_LIST_PATTERN = "%(sequencing_run)s_%(sample_name)s.txt"
SEQAUTO_QC_EXEC_SUMMARY_PATTERN = "exec_stats/%(sample_name)s_qc_summary.txt"
SEQAUTO_QC_GENE_COVERAGE_PATTERN = "bam_stats/samples/%(sample_name)s.per_gene_coverage.tsv.gz"

SEQAUTO_QC_GENE_COVERAGE_STORE_ALL = False
SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL = True

SEQAUTO_MIN_COVERAGE = 20
SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS = []

SEQAUTO_COVERAGE_ENRICHMENT_KITS = []
SEQAUTO_FAKE_SAMPLE_ENRICHMENT_KITS_DF = None  # Used to be able to load gake test data - leave None to load from system

# Links from a SequencingRun page out to external systems (e.g. QC dashboards). Each entry is a dict with:
#   label: short link text to display (e.g. "RAVEN")
#   url_pattern: %()s template formatted with SequencingRun.get_params() — supports %(sequencing_run)s,
#                %(original_sequencing_run)s, %(flowcell_id)s, %(enrichment_kit)s, %(experiment)s
#   min_date: optional ISO date string "YYYY-MM-DD" — runs with no date or an earlier date are skipped
#   exclude_enrichment_kits: optional list of enrichment kit names (case-insensitive exact match);
#                            if the run's enrichment_kit name matches any, the link is skipped
SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS = []
