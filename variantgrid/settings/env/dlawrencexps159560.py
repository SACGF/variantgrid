from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

# import all the base settings #
ROLLBAR['enabled'] = False

#SEQAUTO_ENABLED = True
AVATAR_CACHE_ENABLED = False  # So we can change etc
SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS = ["SAPOrderNumber", "Sex", "Panel", "R1kDVersion", "RunID", "CaptureID", "RunReference"]


ANALYSIS_TEMPLATES_SAPATH = True
ANALYSIS_TEMPLATES_RUNX1 = True

COMPANY = "SA_Pathology"  # Used for gene lists

# Automatically run after import
FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = [
    'annotation.tasks.cohort_sample_gene_damage_counts.CalculateCohortSampleGeneDamageCountsTask',
]

ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'

ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.expanduser("~/bin/perlbrew_runner.sh")
ANNOTATION_VEP_ARGS = ["--buffer_size", "1000"]  # default = 5000
ANNOTATION[BUILD_GRCH38]["enabled"] = True
ANNOTATION[BUILD_GRCH38]["annotation_consortium"] = "RefSeq"

# On laptop I get all kinds of errors using BigWig files, so just turn off
ANNOTATION[BUILD_GRCH37]["vep_config"].update({
    'phastcons100way': None,
    'phastcons46way': None,
    'phylop100way': None,
    'phylop46way': None,
    "spliceai_snv": "annotation_data/GRCh37/spliceai_scores.raw.snv.head_100.hg19.vcf.gz",
    "spliceai_indel": "annotation_data/GRCh37/spliceai_scores.raw.indel.head_100.hg19.vcf.gz",
})
ANNOTATION[BUILD_GRCH38]["vep_config"].update({
    'phastcons100way': None,
    'phastcons30way': None,
    'phylop100way': None,
    'phylop30way': None,
    "spliceai_snv": "annotation_data/GRCh38/spliceai_scores.raw.snv.head_100.hg38.vcf.gz",
    "spliceai_indel": "annotation_data/GRCh38/spliceai_scores.raw.indel.head_100.hg38.vcf.gz",
})

LIFTOVER_NCBI_REMAP_ENABLED = True
LIFTOVER_NCBI_REMAP_PERLBREW_RUNNER_SCRIPT = None  # Use system Perl

_SA_PATH_ENRICHMENT_KITS = [{"name": "idt_gmp_focus"}, {"name": "idt_exome"}]
SEQAUTO_COVERAGE_ENRICHMENT_KITS = _SA_PATH_ENRICHMENT_KITS
GENE_GRID_DEFAULT_ENRICHMENT_KITS = _SA_PATH_ENRICHMENT_KITS
PATHOLOGY_TEST_SORTED_ENRICHMENT_KITS = _SA_PATH_ENRICHMENT_KITS


PEDIGREE_MADELINE2_COMMAND = "madeline2"

URLS_NAME_REGISTER["view_patient_contact_tab"] = True

# CLINGEN_ALLELE_REGISTRY_DOMAIN = "http://reg.test.genome.network"
GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1

_SAPATHOLOGY_MODE = True
_SHARIANT_MODE = False
_RUNX1_MODE = False

SOMALIER["enabled"] = True
# SOMALIER["annotation_base_dir"] = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "somalier")

USE_OIDC = False

if _SAPATHOLOGY_MODE:
    SEQAUTO_ENABLED = True
    SAPATH_APP = 'sapath.apps.SapathConfig'
    INSTALLED_APPS += [SAPATH_APP]
    CELERY_IMPORTS += ('sapath.tasks.import_helix_task',)
    PATHOLOGY_TEST_REQUESTS_REDIRECT_URL = "sapathology_test_requests"
    SAPATHOLOGY_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "sapathology_static")
    if os.path.exists(SAPATHOLOGY_STATIC_FILES_DIR):
        STATICFILES_DIRS = (SAPATHOLOGY_STATIC_FILES_DIR,) + STATICFILES_DIRS

    SAPATHOLOGY_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/sapathology_templates")
    if os.path.exists(SAPATHOLOGY_TEMPLATES_DIR):
        TEMPLATES[0]["DIRS"].insert(0, SAPATHOLOGY_TEMPLATES_DIR)
    PATHOLOGY_TESTS_ENABLED = True
elif _SHARIANT_MODE:

    VARIANT_DETAILS_SHOW_ANNOTATION = False
    VARIANT_CLASSIFICATION_STATS_USE_SHARED = True  # False=Use visible to user. True = Shared

    SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
    STATICFILES_DIRS = (SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

    SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
    TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
    SITE_NAME = "Shariant"
    LOGIN_REDIRECT_URL = '/snpdb/index'

    SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN = True  # set True to only find classified variants
    VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN = False
    VARIANT_MANUAL_CREATE_BY_NON_ADMIN = False

    VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH = True
    VIEW_GENE_SHOW_HOTSPOT_GRAPH = False
    USER_SETTINGS_SHOW_GROUPS = False

    URLS_NAME_REGISTER.update({"classification_dashboard": True,
                               "classification_import_tool": True})

    URLS_NAME_REGISTER.update({  # Disable selected urls
        # Variants
        "variants": False,
        "variant_tags": False,
        "manual_variant_entry": False,
        "variantopedia_wiki": False,
    })

elif _RUNX1_MODE:
    RUNX1_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "runx1_static")
    STATICFILES_DIRS = (RUNX1_STATIC_FILES_DIR,) + STATICFILES_DIRS

    RUNX1_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/runx1_templates")
    TEMPLATES[0]["DIRS"].insert(0, RUNX1_TEMPLATES_DIR)

    PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD = "GRCh37"
    PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID = None  # Was 1
    PUBLIC_SAMPLE_GENE_MATRIX_SHOW_PRIVATE_SAMPLES = True
    PUBLIC_SAMPLE_GENE_MATRIX_TYPE = 'RUNX1_classified_damage'
    PUBLIC_SAMPLE_GENE_MATRIX_HIGHLIGHT_GENE_SYMBOLS = ["RUNX1"]

    LOGIN_REDIRECT_URL = '/snpdb/index'
    COHORT_SAMPLE_GENE_DAMAGE_COUNTS = 'annotation.tasks.cohort_sample_gene_damage_counts.CalculateCohortSampleGeneDamageCountsTask'
    FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = [COHORT_SAMPLE_GENE_DAMAGE_COUNTS]
    SITE_NAME = "RUNX1db"

    PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB = True
