# so we don't pull in 'variantgrid.py' in this dir with Python 2.7.2

from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import


# import all the base settings #
ROLLBAR['enabled'] = False

USE_DJANGO_DEBUG_TOOLBAR = False
if USE_DJANGO_DEBUG_TOOLBAR:
    INSTALLED_APPS += ['debug_toolbar']
    MIDDLEWARE = MIDDLEWARE + ('debug_toolbar.middleware.DebugToolbarMiddleware',)
    INTERNAL_IPS = [
        '127.0.0.1',
    ]

SEQAUTO_ENABLED = True
#SEQAUTO_LOAD_GENE_COVERAGE=False

SEQAUTO_CONTROL_SAMPLE_REGEX = NO_DNA_CONTROL_REGEX
VCF_IMPORT_NO_DNA_CONTROL_SAMPLE_REGEX = NO_DNA_CONTROL_REGEX

ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")
ANNOTATION_REFERENCE_BASE_DIR = "/media/dlawrence/SpinningIron/reference"
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "VEP")
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "ensembl-vep")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")
_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "fasta")

ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'

ANNOTATION[BUILD_GRCH37].update({
    "annotation_consortium": "RefSeq",
    "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
})
ANNOTATION[BUILD_GRCH38].update({
    "annotation_consortium": "RefSeq",
    "enabled": True,
    "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.39_GRCh38.p13_genomic.fna.gz")
})

COMPANY = "SA_Pathology"  # Used for gene lists

SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS = ["SAPOrderNumber", "Sex", "Panel", "R1kDVersion", "RunID", "CaptureID", "RunReference"]

_SA_PATH_ENRICHMENT_KITS = [{"name": "roche_1k_disease", "version": 6}, {"name": "medical_exomes"}]
SEQAUTO_COVERAGE_ENRICHMENT_KITS = _SA_PATH_ENRICHMENT_KITS
GENE_GRID_DEFAULT_ENRICHMENT_KITS = _SA_PATH_ENRICHMENT_KITS
# Fields must be from GoldCoverageSummary and COLUMNS + LABELS must line up!
GENE_GRID_ENRICHMENT_KIT_COLUMNS = ['depth_20x_5th_percentile']
GENE_GRID_ENRICHMENT_KIT_COLUMN_TOOL_TIPS = ["original_transcript_id"]
GENE_GRID_ENRICHMENT_KIT_COLUMN_LABELS = ["% at 20x*"]
GENE_GRID_ENRICHMENT_KIT_COLUMN_LABEL_TOOL_TIPS = ["% at 20x for 5th percentile of gold runs (ie expected worst case)"]

PATHOLOGY_TEST_SORTED_ENRICHMENT_KITS = _SA_PATH_ENRICHMENT_KITS
GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 4

LIFTOVER_NCBI_REMAP_ENABLED = True
LIFTOVER_NCBI_REMAP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

SOMALIER["enabled"] = True
SOMALIER["annotation_base_dir"] = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "somalier")

#CLINGEN_ALLELE_REGISTRY_DOMAIN = "http://reg.test.genome.network"

#DEBUG = False
COMPRESS_ENABLED = False

_SHARIANT_MODE = False
_SAPATHOLOGY_MODE = True
_RUNX1_MODE = False

if _SHARIANT_MODE:
    SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
    STATICFILES_DIRS = (SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

    DISCORDANCE_ENABLED = True
    SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
    TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
    SITE_NAME = "Shariant"
    VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN = False
    VARIANT_MANUAL_CREATE_BY_NON_ADMIN = False

    SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN = True

    # Don't show annotation or samples on variant page - don't want to be responsible for it
    VARIANT_DETAILS_SHOW_ANNOTATION = False
    VARIANT_DETAILS_SHOW_SAMPLES = False
    VARIANT_SHOW_CANONICAL_HGVS = False
    VIEW_GENE_SHOW_WIKI = False
    VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH = True

    URLS_APP_REGISTER.update({"analysis": False,
                              "expression": False,
                              "pathtests": False,
                              "pedigree": False,
                              "seqauto": False})

    URLS_NAME_REGISTER.update({"classification_dashboard": True,
                               "classification_import_tool": True})

elif _SAPATHOLOGY_MODE:

    SAPATH_APP = 'sapath.apps.SapathConfig'
    INSTALLED_APPS += [SAPATH_APP]

    PATHOLOGY_TEST_EXTERNAL_CODE = "SAPOrderNumber"
    SAPATHOLOGY_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "sapathology_static")
    if os.path.exists(SAPATHOLOGY_STATIC_FILES_DIR):
        STATICFILES_DIRS = (SAPATHOLOGY_STATIC_FILES_DIR,) + STATICFILES_DIRS

    SAPATHOLOGY_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/sapathology_templates")
    if os.path.exists(SAPATHOLOGY_TEMPLATES_DIR):
        TEMPLATES[0]["DIRS"].insert(0, SAPATHOLOGY_TEMPLATES_DIR)

    INITIAL_USER_DATA_PREFIX_KWARGS = {"prefix": '/tau',
                                       "replacement": '\\\\frgeneseq01.imvs.sa.gov.au\\tau'}

    SAPATH_HELIX_USER = get_secret("SAPATH.HELIX.user")
    SAPATH_HELIX_PASSWORD = get_secret("SAPATH.HELIX.password")
    SAPATH_HELIX_GENERATE_CSV_FROM_SQL = False

    PATHOLOGY_TESTS_ENABLED = True
    PATHOLOGY_TEST_REQUESTS_REDIRECT_URL = "sapathology_test_requests"

    PATHOLOGY_TEST_EXTERNAL_CODE = "SAPOrderNumber"
    PATHOLOGY_TEST_CASE_EXTERNAL_CODE = "HelixID"

elif _RUNX1_MODE:

    ANALYSIS_TEMPLATES_RUNX1 = True

    RUNX1_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "runx1_static")
    STATICFILES_DIRS = (RUNX1_STATIC_FILES_DIR,) + STATICFILES_DIRS

    RUNX1_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/runx1_templates")
    TEMPLATES[0]["DIRS"].insert(0, RUNX1_TEMPLATES_DIR)

    PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB = True
    PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD = "GRCh37"
    PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID = 19
    PUBLIC_SAMPLE_GENE_MATRIX_SHOW_PRIVATE_SAMPLES = True
    PUBLIC_SAMPLE_GENE_MATRIX_TYPE = 'RUNX1_classified_damage'
    PUBLIC_SAMPLE_GENE_MATRIX_HIGHLIGHT_GENE_SYMBOLS = ["RUNX1"]

    LOGIN_REDIRECT_URL = '/snpdb/index'
    COHORT_SAMPLE_GENE_DAMAGE_COUNTS = 'annotation.tasks.cohort_sample_gene_damage_counts.CalculateCohortSampleGeneDamageCountsTask'
    print(f"before: FINISH_IMPORT_VCF_STEP_TASKS_CLASSES: {FINISH_IMPORT_VCF_STEP_TASKS_CLASSES}")
    FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = [COHORT_SAMPLE_GENE_DAMAGE_COUNTS]
    SITE_NAME = "RUNX1db"

    PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB = True
    VCF_DOWNLOAD_ADMIN_ONLY = True
