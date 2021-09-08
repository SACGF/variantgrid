"""
See https://bitbucket.org/sacgf/variantgrid/wiki/Annotation%20Setup

"""

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

SYNC_DETAILS = get_secrets("SYNC", ["enabled", "username", "password", "host"])

# import all the base settings #
# Overwrite settings for your system below
ANALYSIS_TEMPLATES_RUNX1 = True

ANNOTATION_ENTREZ_EMAIL = 'your@email.com'
ALLOWED_HOSTS = ["runx1db.runx1-fpd.org", "runx1db.runx1.com", "144.6.226.231", "43.240.97.68", "localhost"]
DEBUG = False
SEND_EMAILS = True

ANNOTATION_BASE_DIR = "/mnt/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "ensembl-vep")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")
_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "fasta")

ANNOTATION[BUILD_GRCH37].update({
    "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.25_GRCh37.p13_genomic.fna.gz"),
})

GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1  # MedEx

PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB = True
PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD = "GRCh37"
PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID = None  # was 5
PUBLIC_SAMPLE_GENE_MATRIX_SHOW_PRIVATE_SAMPLES = True
PUBLIC_SAMPLE_GENE_MATRIX_TYPE = 'RUNX1_classified_damage'
PUBLIC_SAMPLE_GENE_MATRIX_HIGHLIGHT_GENE_SYMBOLS = ["RUNX1"]

LOGIN_REDIRECT_URL = '/snpdb/index'
SITE_NAME = "RUNX1db"
SITE_ID = 7

PEDIGREE_MADELINE2_COMMAND = "madeline2"

SOMALIER["enabled"] = True
SOMALIER["annotation_base_dir"] = os.path.join(ANNOTATION_BASE_DIR, "somalier")

RUNX1_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "runx1_static")
STATICFILES_DIRS = (RUNX1_STATIC_FILES_DIR,) + STATICFILES_DIRS

#AUTHENTICATION (RUNX1 ONLY)
AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
        'OPTIONS': {
            'min_length': 10,
        }
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = [
    'annotation.tasks.cohort_sample_gene_damage_counts.CalculateCohortSampleGeneDamageCountsTask',
]

RUNX1_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/runx1_templates")
TEMPLATES[0]["DIRS"].insert(0, RUNX1_TEMPLATES_DIR)

# Lock down menu - hide some VariantGrid urls / menu
URLS_NAME_REGISTER.update({
    "sequencing_data": False,

    # Variants
    "variants": False,
    "manual_variant_entry": False,
    "variantopedia_wiki": False,

    # patients
    "cases": False,
})
