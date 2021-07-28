"""

Global settings - these are loaded then overwritten in hostname specific files

See __init__.py in this dir for details

"""

from collections import defaultdict
import os
import socket

from library.django_utils.django_secret_key import get_or_create_django_secret_key
from library.git import Git

# if certain user settings are not relevant for the environment, list the columns in this
from variantgrid.settings.components.secret_settings import get_secret, get_secrets

CSRF_FAILURE_VIEW = 'variantgrid.views.csrf_error'

# used by
# python3 manage.py collectstatic_js_reverse
# after you need to refer to a JavaScript file, please checkin the resulting reverse.js
JS_REVERSE_OUTPUT_PATH = './variantgrid/static_files/default_static/django_js_reverse'

# Up 2 more dirs than normal (as we're in variantgrid.settings.components dir)
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SETTINGS_DIR = os.path.dirname(_THIS_DIR)
BASE_DIR = os.path.dirname(os.path.dirname(SETTINGS_DIR))

PRIVATE_DATA_ROOT = os.path.join(BASE_DIR, "data")
UPLOAD_RELATIVE_PATH = "data/uploads"  # Needed for FileSystemStorage
UPLOAD_DIR = os.path.join(BASE_DIR, UPLOAD_RELATIVE_PATH)
UPLOAD_ENABLED = True  # This disables uploading files or creating variants (used eg if Redis is out of sync)

# Absolute filesystem path to the directory that will hold GLOBALLY VISIBLE user-uploaded files.
# Example: "/var/www/example.com/media/"
MEDIA_ROOT = os.path.join(BASE_DIR, 'media_root')
MEDIA_URL = '/media/'

PYTHON_COMMAND = "python3.8"
MANAGE_COMMAND = [PYTHON_COMMAND, os.path.join(BASE_DIR, "manage.py")]

# Need 5x as many as largest cohort for CohortNode zygosity query
DATA_UPLOAD_MAX_NUMBER_FIELDS = 5000

# Nightly task to fix missing GRCh37/38 representations
ALLELE_VALIDATION = False

# if None, discordance emails wont be sent
DISCORDANCE_EMAIL = None
ACCOUNTS_EMAIL = None
# If you change this value you should run 'recalc' for all ClinicalContexts in admin
DISCORDANCE_ENABLED = False
# How long you have to update flags after a discordance is closed for them to still
# be considered in the report
DISCORDANCE_REPORT_LEEWAY = 14

DEBUG = True
# If SEND_EMAILS is False, all emails that would go through EmailLog will still be record but wont be sent
# Good for test environments where you don't want to spam test accounts
SEND_EMAILS = False
# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ["localhost"]

ADMINS = (
    # ('Your Name', 'your_email@example.com'),
)

AVATAR_PROVIDERS = (
    'avatar.providers.PrimaryAvatarProvider',
    'library.django_utils.avatar.SpaceThemedAvatarProvider',
)

MANAGERS = ADMINS

BACKEND_ENGINE = "postgres"

CONN_MAX_AGE = 60  # Reuse DB connections

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': get_secret('DB.name'),
        'USER': get_secret('DB.user'),
        'PASSWORD': get_secret('DB.password'),
        'HOST': get_secret('DB.host'),
        # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
        'PORT': get_secret('DB.port'),  # Set to empty string for default.
    }
}

CACHE_HOURS = 48
TIMEOUT = 60 * 60 * CACHE_HOURS
REDIS_PORT = 6379
CACHE_VERSION = 31  # increment to flush caches (eg if invalid due to upgrade)
CACHES = {
    'default': {
        "BACKEND": "redis_cache.RedisCache",
        "LOCATION": "redis://127.0.0.1:%d/1" % REDIS_PORT,
        'TIMEOUT': TIMEOUT,
        'VERSION': CACHE_VERSION,
    },
    'debug-panel': {
        "BACKEND": "redis_cache.RedisCache",
        "LOCATION": "redis://127.0.0.1:%d/1" % REDIS_PORT,
        'TIMEOUT': TIMEOUT,
        'OPTIONS': {
            'MAX_ENTRIES': 200
        }
    },
}

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
TIME_ZONE = 'Australia/Adelaide'
DATE_FORMAT = "%-d/%-m/%y"

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = False

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# If you set this to True (and setup everything else required) the app
# will do auth through OIDC
USE_OIDC = False
OIDC_REQUIRED_GROUP = None
OIDC_USER_SERVICES = None

VARIANTGRID_APP_DIR = os.path.join(BASE_DIR, "variantgrid")

### annotation

VG_REFERENCE_DIR = os.path.join(VARIANTGRID_APP_DIR, "data", "reference")
ANNOTATION_BASE_DIR = "/data/annotation"
ANNOTATION_VEP_FAKE_VERSION = False  # Overridden in unit tests to not call VEP to get version
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = None  # os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

# I've had VEP hang on me when running --fork so by default we run in small batches
# This causes a small amount of overhead obtaining an AnnotationRangeLock
# If you get ERROR: Forked process(es) died: read-through of cross-process communication detected
# You may want to set buffer_size in ANNOTATION_VEP_ARGS below
# @see https://github.com/Ensembl/ensembl-vep/issues/150
ANNOTATION_VEP_FORK = 1
# get_unannotated_count_min_max does quick queries to try and get VEP batch sizes within a range
# If it gets below min, it does a slower query to get range lock.
# The variant table is usually ~55% alt variants but may be different due to data or if you've deleted records
ANNOTATION_VEP_BATCH_MIN = 5000  # Dont' set too low due to overhead of running pipeline etc
ANNOTATION_VEP_BATCH_MAX = 100000  # Set to None to do all in 1 job (probably want to set FORK higher)
ANNOTATION_VEP_ARGS = []  # ["--buffer_size", "1000"] # default = 5000
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "ensembl-vep")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")
# @see https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick_order
ANNOTATION_VEP_PICK_ORDER = None
ANNOTATION_VEP_DISTANCE = None  # VEP --distance arg (default=5000) - how far up/down to assign to a transcript

_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "fasta")

BUILD_GRCH37 = "GRCh37"
BUILD_GRCH38 = "GRCh38"

ANNOTATION = {
    BUILD_GRCH37: {
        "enabled": True,
        "annotation_consortium": "Ensembl",
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.25_GRCh37.p13_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "cytoband": os.path.join(VG_REFERENCE_DIR, "hg19", "cytoband.hg19.txt.gz"),

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "cosmic": "annotation_data/GRCh37/CosmicCodingMuts.normal.grch37.vcf.gz",
            "dbnsfp": "annotation_data/GRCh37/dbNSFP4.0a.grch37.stripped.gz",
            "dbscsnv": "annotation_data/GRCh37/dbscSNV1.1_GRCh37.txt.gz",
            "gnomad2": "annotation_data/GRCh37/gnomad2.1.1_GRCh37_combined_af.vcf.bgz",
            "mastermind": "annotation_data/GRCh37/mastermind_cited_variants_reference-2021.04.02-grch37.vcf.gz",
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh37/hg19.100way.phastCons.bw",
            'phastcons30way': None,  # n/a for GRCh37
            'phastcons46way': "annotation_data/GRCh37/hg19.phastCons46way.placental.bw",
            'phylop100way': "annotation_data/GRCh37/hg19.100way.phyloP100way.bw",
            'phylop30way': None,  # n/a for GRCh37
            'phylop46way': "annotation_data/GRCh37/hg19.phyloP46way.placental.bw",
            "repeatmasker": "annotation_data/GRCh37/repeatmasker_hg19.bed.gz",
            "spliceai_snv": "annotation_data/GRCh37/spliceai_scores.raw.snv.hg19.vcf.gz",
            "spliceai_indel": "annotation_data/GRCh37/spliceai_scores.raw.indel.hg19.vcf.gz",
            "topmed": "annotation_data/GRCh37/TOPMED_GRCh37.vcf.gz",
            "uk10k": "annotation_data/GRCh37/UK10K_COHORT.20160215.sites.vcf.gz",
        }
    },
    # GRCh38 is NOT enabled by default - overwrite "enabled" in your server settings to use
    BUILD_GRCH38: {
        "enabled": False,
        "annotation_consortium": "Ensembl",
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "cytoband": os.path.join(VG_REFERENCE_DIR, "hg38", "cytoband.hg38.txt.gz"),

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "cosmic": "annotation_data/GRCh38/CosmicCodingMuts.normal.grch38.vcf.gz",
            "dbnsfp": "annotation_data/GRCh38/dbNSFP4.0a.grch38.stripped.gz",
            "dbscsnv": "annotation_data/GRCh38/dbscSNV1.1_GRCh38.txt.gz",
            "gnomad2": "annotation_data/GRCh38/gnomad2.1.1_GRCh38_combined_af.vcf.bgz",
            "gnomad3": "annotation_data/GRCh38/gnomad3.1_GRCh38_merged.vcf.bgz",
            "mastermind": "annotation_data/GRCh38/mastermind_cited_variants_reference-2021.04.02-grch38.vcf.gz",
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh38/hg38.phastCons100way.bw",
            'phastcons30way': "annotation_data/GRCh38/hg38.phastCons30way.bw",
            'phylop100way': "annotation_data/GRCh38/hg38.phyloP100way.bw",
            'phylop30way': "annotation_data/GRCh38/hg38.phyloP30way.bw",
            "repeatmasker": "annotation_data/GRCh38/repeatmasker_hg38.bed.gz",
            "spliceai_snv": "annotation_data/GRCh38/spliceai_scores.raw.snv.hg38.vcf.gz",
            "spliceai_indel": "annotation_data/GRCh38/spliceai_scores.raw.indel.hg38.vcf.gz",
            "topmed": "annotation_data/GRCh38/TOPMED_GRCh38_20180418.vcf.gz",
            "uk10k": "annotation_data/GRCh38/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz",
        }
    },
}

# Build independent config
ANNOTATION_VCF_DUMP_DIR = os.path.join(PRIVATE_DATA_ROOT, 'annotation_dump')
# Admin email used in PubMed queries to contact before throttling/banning
ANNOTATION_ENTREZ_EMAIL = get_secret("ENTREZ.email")  # Automatically set in in annotation.apps.AnnotationConfig
ANNOTATION_ENTREZ_API_KEY = get_secret("ENTREZ.api_key")
ANNOTATION_PUBMED_GENE_SYMBOL_COUNT_CACHE_DAYS = 30
ANNOTATION_PUBMED_SEARCH_TERMS_ENABLED = False

MUTATIONAL_SIGNATURE_CALCULATOR = "Mutational Signature Calculator"
MUTATIONAL_SIGNATURE_ITERATIONS = 100
MUTATIONAL_SIGNATURE_SAMPLING_FRACTION = 0.8
MUTATIONAL_SIGNATURE_DATA_DIR = os.path.join(VG_REFERENCE_DIR, "mutational_signatures")
MUTATIONAL_SIGNATURE_DATA_FILE = os.path.join(MUTATIONAL_SIGNATURE_DATA_DIR, "signatures_probabilities.sorted.txt")
MUTATIONAL_SIGNATURE_INFO_FILE = os.path.join(MUTATIONAL_SIGNATURE_DATA_DIR, "signature_analysis_data.formatted.txt")

CACHED_WEB_RESOURCE_CLINGEN_DISEASE_VALIDITY = "ClinGenDiseaseValidity"
CACHED_WEB_RESOURCE_GNOMAD_GENE_CONSTRAINT = "GnomADGeneConstraint"
CACHED_WEB_RESOURCE_HGNC = "HGNC"
CACHED_WEB_RESOURCE_PANEL_APP_AUSTRALIA_PANELS = "PanelApp Australia Panels"
CACHED_WEB_RESOURCE_PANEL_APP_ENGLAND_PANELS = "Genomics England PanelApp Panels"
CACHED_WEB_RESOURCE_PFAM = "Pfam"
CACHED_WEB_RESOURCE_REFSEQ_GENE_SUMMARY = "RefSeq Gene Summary"
CACHED_WEB_RESOURCE_UNIPROT = "UniProt"

ANNOTATION_CACHED_WEB_RESOURCES = [
    CACHED_WEB_RESOURCE_GNOMAD_GENE_CONSTRAINT,
    CACHED_WEB_RESOURCE_HGNC,
    CACHED_WEB_RESOURCE_PANEL_APP_AUSTRALIA_PANELS,
    CACHED_WEB_RESOURCE_PANEL_APP_ENGLAND_PANELS,
    CACHED_WEB_RESOURCE_PFAM,
    CACHED_WEB_RESOURCE_REFSEQ_GENE_SUMMARY,
    CACHED_WEB_RESOURCE_UNIPROT,
    CACHED_WEB_RESOURCE_CLINGEN_DISEASE_VALIDITY,
]

VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES = ['refseq_transcript_accession', 'ensembl_transcript_accession']
VARIANT_ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS = not DEBUG
# If true, then if we don't have a specific transcript version, we'll match it to the closest one we can
VARIANT_TRANSCRIPT_VERSION_BEST_ATTEMPT = True

VARIANT_ZYGOSITY_GLOBAL_COLLECTION = "global"

PREFER_ALLELE_LINKS = False

CLINGEN_ALLELE_REGISTRY_DOMAIN = "http://reg.genome.network"
CLINGEN_ALLELE_REGISTRY_MAX_RECORDS = 2000
CLINGEN_ALLELE_REGISTRY_LOGIN = get_secret("CLINGEN_ALLELE_REGISTRY.login")
CLINGEN_ALLELE_REGISTRY_PASSWORD = get_secret("CLINGEN_ALLELE_REGISTRY.password")
CLINGEN_ALLELE_REGISTRY_MAX_MANUAL_REQUESTS = 10_000  # On nodes and VCFs

NO_DNA_CONTROL_REGEX = "(^|[^a-zA-Z])NDC([^a-zA-Z]|$)"  # No DNA Control - eg _NDC_ or -NDC_

VCF_DOWNLOAD_ADMIN_ONLY = False
VCF_IMPORT_DELETE_TEMP_FILES_ON_SUCCESS = not DEBUG
VCF_IMPORT_CREATE_COHORT_FROM_MULTISAMPLE_VCFS = True
VCF_IMPORT_NO_DNA_CONTROL_SAMPLE_REGEX = None
VCF_IMPORT_FILE_SPLIT_ROWS = 50000
VCF_IMPORT_STORE_GVCF_NON_VAR_BLOCKS = False
VCF_IMPORT_VT_COMMAND = "vt"  # Needs to be installed and in path

COMPANY = None  # Used for gene list categories

GENERATED_DIR = os.path.join(MEDIA_ROOT, 'generated')

PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB = False
IMPORT_PROCESSING_DIR = os.path.join(PRIVATE_DATA_ROOT, 'import_processing')

# @see https://github.com/SACGF/variantgrid/wiki/Liftover
LIFTOVER_CLASSIFICATIONS = True
LIFTOVER_TO_CHROMOSOMES_ONLY = True  # False = Liftover to alt/patches
LIFTOVER_DBSNP_ENABLED = False  # Default=False - doesn't work so well due to dbSNP IDs being for loci
LIFTOVER_NCBI_REMAP_ENABLED = False
LIFTOVER_NCBI_REMAP_PERLBREW_RUNNER_SCRIPT = None  # os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

PANEL_APP_CACHE_DAYS = 7  # Automatically re-check after this time
PANEL_APP_CHECK_ENABLED = False

# Non-authenticated (no login) sample gene matrix
PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD = None  # Must be set if system has multiple genome builds
PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID = None  # Sample gene matrix of these genes
PUBLIC_SAMPLE_GENE_MATRIX_SHOW_PRIVATE_SAMPLES = False  # True = show all samples
PUBLIC_SAMPLE_GENE_MATRIX_TYPE = 'Damage'
PUBLIC_SAMPLE_GENE_MATRIX_HIGHLIGHT_GENE_SYMBOLS = []
PUBLIC_LAB_LOCATION_PAGE_VISIBLE = False

PROCESSED_BED_FILES_DIR = os.path.join(PRIVATE_DATA_ROOT, 'processed_bed_files')
INTERSECT_BED_SCRIPT = os.path.join(BASE_DIR, 'scripts', 'intersect_bed_and_upload_variant_collection.sh')

PUBLIC_GROUP_NAME = "public"
LOGGED_IN_USERS_GROUP_NAME = "logged_in_users"

# key/value = Organization.group_name : lab group name pattern
# Org must already exist, lab pattern is filled with User values (if you want to create a group for each user)
USER_CREATE_ORG_LABS = {
    # "test_organization": "test_lab",
    # "test_organization": "user_group_%(username)s",
}


# To use SeqAuto, your settings need to have:
# "from variantgrid.settings.defaults.seqauto_default_settings import *"
# after including this file
SEQAUTO_ENABLED = False

# Occasionally turn off if it helps with testing
# (Server Status page will freeze up if this is true and celery is not on)
CELERY_ENABLED = True

GENE_GRID_DEFAULT_ENRICHMENT_KITS = []
# Fields must be from GoldCoverageSummary and COLUMNS + LABELS must line up!
GENE_GRID_ENRICHMENT_KIT_COLUMNS = ['mean', 'min_mean', 'depth_20x_5th_percentile']
GENE_GRID_ENRICHMENT_KIT_COLUMN_TOOL_TIPS = ["original_transcript"] * len(GENE_GRID_ENRICHMENT_KIT_COLUMNS)
GENE_GRID_ENRICHMENT_KIT_COLUMN_LABELS = ["Mean", "Min Mean", "Depth 20x (5th percentile)"]
GENE_GRID_ENRICHMENT_KIT_COLUMN_LABEL_TOOL_TIPS = [None, None, None]

GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = None

DEFAULT_COLUMNS_NAME = 'Default columns'

DEFAULT_ENRICHMENT_KIT_LEFT_PADDING = 0
DEFAULT_ENRICHMENT_KIT_RIGHT_PADDING = 0

ANALYSIS_DUAL_SCREEN_MODE_FEATURE_ENABLED = False  # Currently broken
ANALYSIS_TEMPLATES_AUTO_SAMPLE = "Sample tab auto analysis"
ANALYSIS_WARN_IF_NO_QC_GENE_LIST_MESSAGE = None  # disabled by default
ANALYSIS_NODE_CACHE_Q = True

VARIANT_STANDARD_BASES_ONLY = True  # True to reject anything other than A, C, G, T
VARIANT_SHOW_CANONICAL_HGVS = True

VARIANT_CLASSIFICATION_SUPPORTED_TRANSCRIPTS = {"NR", "NM", "NC", "ENST"}
VARIANT_CLASSIFICATION_MATCH_VARIANTS = True  # exists only so we can turn it off during testing
VARIANT_CLASSIFICATION_REQUIRE_OVERWRITE_NOTE = True
VARIANT_CLASSIFICATION_AUTOFUZZ_AGE = False

VARIANT_CLASSIFICATION_DASHBOARD_SIZE = 50
VARIANT_CLASSIFICATION_RECLASSIFICATION_EMAIL = True
VARIANT_CLASSIFICATION_ID_FILTER = True
VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME = True
VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN = True  # Should Allele origin (e.g. germline/somatic) be shown on the grid
VARIANT_CLASSIFICATION_STATS_USE_SHARED = False  # False=Use visible to user. True = Shared
VARIANT_CLASSIFICATION_GRID_SHOW_PHGVS = True
VARIANT_CLASSIFICAITON_SHOW_SPECIMEN_ID = True
VARIANT_CLASSIFICATION_NEW_GROUPING = False

# Require people to click "my sample's not here" (ie encourage them to find it)
VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE = True
VARIANT_CLASSIFICATION_WEB_FORM_CREATE = True
VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN = True
VARIANT_CLASSIFICATION_WEB_FORM_CREATE_ALLOW_NO_VARIANT = True  # Can create purely from HGVS

VARIANT_CLASSIFICATION_FILE_ATTACHMENTS = True  # allow users to attach files to classifications

VARIANT_CLASSIFICATION_MAX_FULL_ALLELE_LENGTH = 100  # Used for MVL export, for general display limit is 10

PATHOLOGY_TESTS_ENABLED = False
PATHOLOGY_TEST_REQUESTS_REDIRECT_URL = None
PATHOLOGY_TEST_EXTERNAL_CODE = None
PATHOLOGY_TEST_SORTED_ENRICHMENT_KITS = []  # Tests get automatically compared to
PATHOLOGY_TEST_CASE_EXTERNAL_CODE = None

PEDIGREE_MIN_COHORT_SAMPLE_MATCHES_FOR_AUTO_MATCH = 3
PEDIGREE_MADELINE2_COMMAND = None  # Install https://madeline.med.umich.edu/madeline/ set this to "madeline2"

INITIAL_USER_DATA_PREFIX_KWARGS = {}  # Create UserDataPrefix object to setup IGV for new users

USER_SETTINGS_SHOW_GROUPS = True

REDIS_PIPELINE_SIZE = 100000
SQL_BATCH_INSERT_SIZE = 50000
SQL_SCRIPTS_DIR = os.path.join(BASE_DIR, "dbscripts")
SITE_NAME = "VariantGrid"

SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN = False  # set True to only find classified variants
SEARCH_SUMMARY = True
SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR = True
SILENCED_SYSTEM_CHECKS = ['models.E006']  # 'captcha.recaptcha_test_key_error']
SITE_ID = 2
SITE_MESSAGE = None  # displayed at the top of all logged-in pages

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = os.path.join(VARIANTGRID_APP_DIR, "sitestatic", "static")

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/static/'

# modify this configuration per environment
rollbar_access_token = get_secret("ROLLBAR.access_token")
rollbar_client_access_token = get_secret("ROLLBAR.client_access_token")

ROLLBAR = {
    'access_token': rollbar_access_token,
    'client_access_token': rollbar_client_access_token,
    'environment': socket.gethostname().lower().split('.')[0].replace('-', ''),
    'enabled': rollbar_access_token and rollbar_client_access_token,  # set to false in environments to disable rollbar
    'branch': 'master',
    'root': BASE_DIR,
    'capture_username': True,
    'code_version': Git(BASE_DIR).hash,
}

SLACK = {
    'enabled': get_secret('SLACK.enabled'),
    'admin_callback_url': get_secret('SLACK.admin_callback_url')
}
# if true, automated health checks will post to Slack if enabled
HEALTH_CHECK_ENABLED = True

SERVER_MIN_DISK_WARNING_GIGS = 1
USER_FEEDBACK_ENABLED = True  # note that Rollbar enabled must also be true to enable user feedback

HEARTBEAT_URL = None  # URL to ping to prove scheduler is still alive

STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "default_static")

# Additional locations of static files
STATICFILES_DIRS = (
    STATIC_FILES_DIR,
    # Use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    #    'django.contrib.staticfiles.finders.DefaultStorageFinder',
    'compressor.finders.CompressorFinder'
)

STATICFILES_STORAGE = "django.contrib.staticfiles.storage.ManifestStaticFilesStorage"

# Needs to be unique and not checked into source control, so make
# django_secret_key.txt in this dir (which is hidden via .gitignore)
SECRET_KEY = get_or_create_django_secret_key(SETTINGS_DIR)

TAG_REQUIRES_CLASSIFICATION = "RequiresClassification"  # tags can't have spaces

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            os.path.join(VARIANTGRID_APP_DIR, "templates/default_templates"),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.request',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',
                "snpdb.processors.settings_context_processor",
            ],
        },
    },
]

# only used if you install and enable querycount.middleware.QueryCountMiddleware
QUERYCOUNT = {
    'DISPLAY_DUPLICATES': 5,
}

MIDDLEWARE = (
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'global_login_required.GlobalLoginRequiredMiddleware',  # Must be after other auth middleware
    'library.django_utils.rollbar_middleware.CustomRollbarNotifierMiddleware',
    'htmlmin.middleware.HtmlMinifyMiddleware',
    'htmlmin.middleware.MarkRequestMiddleware',
    'threadlocals.middleware.ThreadLocalMiddleware'
    # 'querycount.middleware.QueryCountMiddleware',
    # 'mozilla_django_oidc.middleware.SessionRefresh',
    # 'debug_panel.middleware.DebugPanelMiddleware',
)
HTML_MINIFY = True

ROOT_URLCONF = 'variantgrid.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'variantgrid.wsgi.application'

# We want a class to be able to register itself as a factory
# The following will be imported, and looked for ImportTaskFactory.__subclasses__()
IMPORT_TASK_FACTORY_IMPORTS = (
    'annotation.import_task_factories',
    'upload.import_task_factories.import_task_factories',
)

# To be able to automatically do things after VCF import - eg damage counts
FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = []

# Turn ON in production!
CACHE_GENERATED_FILES = True

REST_FRAMEWORK = {
    'DEFAULT_PERMISSION_CLASSES': (
        'rest_framework.permissions.IsAuthenticated',
    ),
    'EXCEPTION_HANDLER': 'rollbar.contrib.django_rest_framework.post_exception_handler'
    #    'PAGE_SIZE': 10,
    #    'DEFAULT_PAGINATION_CLASS': 'rest_framework.pagination.LimitOffsetPagination',
}

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',  # default
    'guardian.backends.ObjectPermissionBackend',
)
ANONYMOUS_USER_ID = -1

INBOX_ENABLED = True
HELP_URL = "https://variantgrid.readthedocs.io/en/latest/"
LOGIN_URL = '/accounts/login/'
LOGIN_REDIRECT_URL = "dashboard"  # Remember to put @terms_required around this view
LOGIN_USERNAME_PLACEHOLDER = None
LOGOUT_REDIRECT_URL = '/'

ACCOUNT_ACTIVATION_DAYS = 7  # One-week activation window

INSTALLED_APPS = [
    'compressor',
    'avatar',
    'vcauth',
    'datetimeutc',
    'django_admin_json_editor',
    'django.contrib.humanize',
    'django.contrib.sites',
    'captcha',
    'registration',
    'django.contrib.auth',
    'mozilla_django_oidc',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    # 3rd party libraries
    'dal',  # Django Autocomplete Light v3
    'dal_select2',  # DAL Plugin
    'django_messages',
    'django_dag',
    'django_js_reverse',
    'django_starfield',
    'django_extensions',
    'djgeojson',
    'easy_thumbnails',
    'guardian',
    'jfu',
    'leaflet',
    'rest_framework',
    'termsandconditions',
    'crispy_forms',  # used to make bootstrap compatible forms
    #    'debug_toolbar',
    #    'debug_panel',
    # Internal apps
    'analysis.apps.AnalysisConfig',
    'annotation.apps.AnnotationConfig',
    'ontology',
    'eventlog',
    'expression',
    'flags',
    'manual',
    'email_manager',
    'sync',
    'genes.apps.GenesConfig',
    'pathtests',
    'patients',
    'pedigree',
    'seqauto.apps.SeqautoConfig',
    'uicore.apps.UiCoreConfig',
    'snpdb.apps.SnpdbConfig',
    'upload.apps.UploadConfig',
    'classification.apps.ClassificationConfig',
    'variantopedia',
    # Uncomment the next line to enable admin documentation:
    # 'django.contrib.admindocs',
]

# https://django-crispy-forms.readthedocs.io/en/latest/install.html
# CRISPY_ALLOWED_TEMPLATE_PACKS = ('bootstrap4_neat', 'bootstrap4') # need to do this if you make an alternative template pack
CRISPY_TEMPLATE_PACK = 'bootstrap4'

# https://github.com/iatkin/leaflet-svgicon
LEAFLET_CONFIG = {
    'PLUGINS': {
        'forms': {
            'js': ['/static/js/lib/leaflet-svgicon/svg-icon.js'],
            'auto-include': True,
        },
    }
}

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'root': {
        'level': 'INFO',
        'handlers': ['console'],
    },
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
    },
    'handlers': {
        'null': {
            'level': 'DEBUG',
            'class': 'logging.NullHandler',
        },
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },
        'db': {
            'level': 'DEBUG',
            'class': 'eventlog.loggers.EventLogHandler',
            'formatter': 'verbose'
        },

        #         'mail_admins': {
        #             'level': 'ERROR',
        #             'class': 'django.utils.log.AdminEmailHandler',
        #         }
    },
    'loggers': {
        'django': {
            'handlers': ['console', 'db'],
            'propagate': True,
            'level': 'INFO',
        },
        #        'django.request': {
        #            'handlers': ['mail_admins'],
        #            'level': 'ERROR',
        #            'propagate': False,
        #        },
        'django.db.backends': {
            'handlers': ['console'],
            'propagate': False,
            #            'level':'DEBUG',
        },
    }
}

# Instead of @login_required we use GlobalLoginRequiredMiddleware to block non-auth by
# default then whitelist with @login_not_required or PUBLIC_PATHS setting below.

# Django REST Framework:
# DRF overrides Django's internals, so request.user is AnonymousUser during
# the middleware check (and will fail with GlobalLoginRequiredMiddleware).
# @see https://github.com/encode/django-rest-framework/issues/760#issuecomment-391127616
# Thus for DRF we need to exclude it and add its own check (at view time)
PUBLIC_PATHS = [
    r'^/accounts/.*',  # allow public access to all django registration views,
    r'^/oidc/.*',  # all oidc URLs
    r'^/classification/api/.*'  # REST framework used by command line tools
]

# Both need to be set to enable - and use get_secret in server settings files to keep out of source control
RECAPTCHA_PUBLIC_KEY = ""  # get_secret('RECAPTCHA.public_key')
RECAPTCHA_PRIVATE_KEY = ""  # get_secret('RECAPTCHA.private_key')
REGISTRATION_FORM = "variantgrid.forms.ReCaptchaSignupForm"

SESSION_ENGINE = "django.contrib.sessions.backends.cache"


# Somalier config - see https://github.com/brentp/somalier

SOMALIER = {
    "enabled": False,
    "admin_only": False,
    "vcf_base_dir": os.path.join(PRIVATE_DATA_ROOT, "somalier"),  # Private data
    "report_base_dir": os.path.join(MEDIA_ROOT, "somalier"),  # Static files served by Nginx
    "annotation_base_dir": os.path.join(ANNOTATION_BASE_DIR, "somalier"),
    "annotation": {  # All annotation paths relative to "annotation_base_dir"
        "command": "somalier",
        "ancestry_labels": "ancestry-labels-1kg.tsv",
        "ancestry_somalier_dir": "1kg-somalier",
        "sites": {
            "GRCh37": "sites.GRCh37.vcf.gz",  # No chr prefix on chromosome names
            "GRCh38": "sites.hg38.nochr.vcf.gz",  # No chr prefix on chromosome names
        },
    },
    # Minimums for related samples to appear at bottom of view_sample page
    "relatedness": {
        "min_relatedness": 0.1,
        "min_shared_hets": 1000,
        "min_shared_hom_alts": 1000,
    }
}


# @see https://github.com/SACGF/variantgrid/wiki/URL---Menu-configuration
# Before URLs are registered, the URLS_APP_REGISTER and URLS_NAME_REGISTER are looked up
# To make a whitelist - change the default to False, then add overrides, eg 'url_name' : True for allowed
# To make a blacklist - leave default as True, add 'url_name' : False

# Use this to handle url registration at the app level (to eg block and entire app)
_URLS_APP_REGISTER_DEFAULT = True
_URLS_APP_REGISTER_OVERRIDE = {}  # Keys are app names (eg "seqauto" or "snpdb")
URLS_APP_REGISTER = defaultdict(lambda: _URLS_APP_REGISTER_DEFAULT, _URLS_APP_REGISTER_OVERRIDE)

# This only works with NAMED urls (but we should give them all an URL)
# This is added to context as "url_name_visible", use as per:
# {% if url_name_visible.patients %}
#     <li><a href="{% url 'patients' %}">patients</a></li>
# {% endif %}

_URLS_NAME_REGISTER_DEFAULT = True
# Keys are url names (eg 'view_sample' or 'node_edit')
_URLS_NAME_REGISTER_OVERRIDE = {
    "view_patient_contact_tab": False,
    "classification_import_tool": False,
    "condition_aliases": False,
    "clinvar_exports": False,
    "condition_matchings": False,
    "condition_match_test": False
}
URLS_NAME_REGISTER = defaultdict(lambda: _URLS_NAME_REGISTER_DEFAULT, _URLS_NAME_REGISTER_OVERRIDE)

DEFAULT_TERMS_SLUG = 'site-terms'
ACCEPT_TERMS_PATH = '/terms/accept/'
# TERMS_EXCLUDE_URL_PREFIX_LIST = {'/admin', '/terms'}
# TERMS_EXCLUDE_URL_LIST = {'/', '/termsrequired/', '/logout/', '/securetoo/','/external_help','/external_help?page=about.html',
#    '/external_help?page=terms.html'}
# TERMS_EXCLUDE_URL_CONTAINS_LIST = {}  # Useful if you are using internationalization and your URLs could change per language
TERMS_CACHE_SECONDS = 0
TERMS_EXCLUDE_USERS_WITH_PERM = 'auth.can_skip_t&c'
TERMS_BASE_TEMPLATE = 'base_tc.html'

# Lock down menu
URLS_NAME_REGISTER.update({"classification_dashboard": False,
                           "keycloak_admin": False,
                           "version_diffs": False,
                           "evidence_keys_logged_in": False,
                           "classification_import_upload": False})

VARIANT_DETAILS_SHOW_ANNOTATION = True  # also doubles as GENE_SHOW_ANNOTATION
VARIANT_DETAILS_SHOW_SAMPLES = True
VARIANT_DETAILS_NEARBY_RANGE = 50
VARIANT_VCF_DB_PREFIX = "vg"
VARIANT_MANUAL_CREATE = True
VARIANT_MANUAL_CREATE_BY_NON_ADMIN = True
VARIANT_PK_HASH_USE_REDIS = False

VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH = False
VIEW_GENE_SHOW_HOTSPOT_GRAPH = True
VIEW_GENE_SHOW_WIKI = True

KEY_CLOAK_REALM = None

UNSHARED_FLAG_ENABLED = False

## COMPRESSION SETTINGS
COMPRESS_ENABLED = True

# Standard SECRETS
# Config used to talk directly to keycloak
KEYCLOAK_SYNC_DETAILS = None
# Config used to upload to Shariant
SYNC_DETAILS = None
# Config used if using AWS email
AWS_SES_ACCESS_KEY_ID = None
AWS_SES_SECRET_ACCESS_KEY = None
AWS_SES_REGION = None

# Command line tool to unzip and cat a file to stdout
# for macOS need to set this to gzcat as the default zcat has short-comings
BASH_ZCAT = 'zcat'
# If True, will run a series of bash commands as one long string with Shell=True
# otherwise will pipe each command into the next more safely with Shell=False
POPEN_SHELL = True  # For vcf split - todo put back...

# Bootstrapped themed messages
from django.contrib.messages import constants as messages

# @see https://docs.djangoproject.com/en/3.1/ref/settings/#message-tags
MESSAGE_TAGS = {
    messages.DEBUG: 'alert-info',
    messages.INFO: 'alert-info',
    messages.SUCCESS: 'alert-success',
    messages.WARNING: 'alert-warning',
    messages.ERROR: 'alert-danger',
}


def get_keycloak_sync_secrets() -> dict:
    """
    KEYCLOAK_SYNC_DETAILS = get_keycloak_sync_secrets()
    """
    return get_secrets("KEYCLOAK", ["username", "password", "host", "oauth_url", "client_id"])


def get_shariant_sync_secrets() -> dict:
    """
    SYNC_DETAILS = get_shariant_sync_secrets()
    """
    return get_secrets("SYNC", ["enabled", "username", "password", "host", "oauth_url", "client_id"])


def get_aws_secrets() -> dict:
    """
    aws_ses_dict = get_aws_secrets()
    AWS_SES_ACCESS_KEY_ID, AWS_SES_SECRET_ACCESS_KEY, AWS_SES_REGION = \
        [aws_ses_dict[k] for k in ('AWS_SES_ACCESS_KEY_ID', 'AWS_SES_SECRET_ACCESS_KEY', 'AWS_SES_REGION')]
    """
    return {
        "AWS_SES_ACCESS_KEY_ID": get_secret("AWS.SES.access_key_id"),
        "AWS_SES_SECRET_ACCESS_KEY": get_secret("AWS.SES.secret_access_key"),
        "AWS_SES_REGION": get_secret("AWS.SES.region")
    }


def get_s3_secrets() -> dict:
    """
    aws_s3_dict = get_s3_secrets()
    AWS_S3_ACCESS_KEY_ID, AWS_S3_SECRET_ACCESS_KEY = \
        [aws_s3_dict[k] for k in ("AWS_S3_ACCESS_KEY_ID", "AWS_S3_SECRET_ACCESS_KEY")]
    :return:
    """
    return {
        "AWS_S3_ACCESS_KEY_ID": get_secret("AWS.S3.access_key_id"),
        "AWS_S3_SECRET_ACCESS_KEY": get_secret("AWS.S3.secret_access_key")
    }
