"""

Global settings - these are loaded then overwritten in hostname specific files

See __init__.py in this dir for details

"""

import os
import re
import socket
import sys
from collections import defaultdict

from library.django_utils.django_secret_key import get_or_create_django_secret_key
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.git import Git
# if certain user settings are not relevant for the environment, list the columns in this
from variantgrid.settings.components.secret_settings import get_secret, get_secrets
from variantgrid.settings.components.settings_paths import BASE_DIR, PRIVATE_DATA_ROOT, SETTINGS_DIR, ANNOTATION_BASE_DIR

CSRF_FAILURE_VIEW = 'variantgrid.views.csrf_error'

# used by
# python3 manage.py collectstatic_js_reverse
# after you need to refer to a JavaScript file, please checkin the resulting reverse.js
JS_REVERSE_OUTPUT_PATH = './variantgrid/static_files/default_static/django_js_reverse'

UPLOAD_RELATIVE_PATH = "data/uploads"  # Needed for FileSystemStorage
UPLOAD_DIR = os.path.join(BASE_DIR, UPLOAD_RELATIVE_PATH)
UPLOAD_ENABLED = True  # This disables uploading files or creating variants (eg if out of disk)

# Absolute filesystem path to the directory that will hold GLOBALLY VISIBLE user-uploaded files.
# Example: "/var/www/example.com/media/"
MEDIA_ROOT = os.path.join(BASE_DIR, 'media_root')
MEDIA_URL = '/media/'

PYTHON_COMMAND = sys.executable
MANAGE_COMMAND = [PYTHON_COMMAND, os.path.join(BASE_DIR, "manage.py")]

# Need 5x as many as largest cohort for CohortNode zygosity query
DATA_UPLOAD_MAX_NUMBER_FIELDS = 5000

# if None, discordance emails wont be sent
DISCORDANCE_EMAIL = None

# if None, Emails to admins wont be sent, is the FROM email for emails sent form this server to admins
ADMIN_EMAIL_NOTIFICATION = None
# Enable to have contact us functionality
CONTACT_US_ENABLED = False

ACCOUNTS_EMAIL = None
# If you change this value you should run 'recalc' for all ClinicalContexts in admin
DISCORDANCE_ENABLED = False

ALLELE_ORIGIN_NOT_PROVIDED_BUCKET = "U"
# If allele origin isn't provided, what bucket do we group it in
# "U" for Unknown, "G" for Germline, and "S" for Somatic

UNIT_TEST = sys.argv[1:2] == ['test']
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
AVATAR_THUMB_FORMAT = "PNG"

MANAGERS = ADMINS

CONN_MAX_AGE = 60  # Reuse DB connections

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"

DATABASES = {
    'default': {
        'ENGINE': 'psqlextra.backend',
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
CACHE_VERSION = 39  # increment to flush caches (eg if invalid due to upgrade)
CACHES = {
    'default': {
        "BACKEND": "django.core.cache.backends.redis.RedisCache",
        "LOCATION": "redis://127.0.0.1:%d/1" % REDIS_PORT,
        'TIMEOUT': TIMEOUT,
        'VERSION': CACHE_VERSION,
    },
    'debug-panel': {
        "BACKEND": "django.core.cache.backends.redis.RedisCache",
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

AVAILABLE_TZS = [
    "Australia/ACT",
    "Australia/Adelaide",
    "Australia/Brisbane",
    "Australia/Broken_Hill",
    "Australia/Canberra",
    "Australia/Currie",
    "Australia/Darwin",
    "Australia/Eucla",
    "Australia/Hobart",
    "Australia/LHI",
    "Australia/Lindeman",
    "Australia/Lord_Howe",
    "Australia/Melbourne",
    "Australia/NSW",
    "Australia/North",
    "Australia/Perth",
    "Australia/Queensland",
    "Australia/South",
    "Australia/Sydney",
    "Australia/Tasmania",
    "Australia/Victoria",
    "Australia/West",
    "Australia/Yancowinna"
]

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
MAINTENANCE_MODE = get_secret("SECURITY.maintenance_mode")  # If true, only non-bot admin users will be able to login, currently only works for ODIC
OIDC_REQUIRED_GROUP = None
OIDC_USER_SERVICES = None

VARIANTGRID_APP_DIR = os.path.join(BASE_DIR, "variantgrid")
_VARIANTGRID_REPO_REFERENCE_DIR = os.path.join(VARIANTGRID_APP_DIR, "data", "reference")

MUTATIONAL_SIGNATURE_CALCULATOR = "Mutational Signature Calculator"
MUTATIONAL_SIGNATURE_ITERATIONS = 100
MUTATIONAL_SIGNATURE_SAMPLING_FRACTION = 0.8
MUTATIONAL_SIGNATURE_DATA_DIR = os.path.join(_VARIANTGRID_REPO_REFERENCE_DIR, "mutational_signatures")
MUTATIONAL_SIGNATURE_DATA_FILE = os.path.join(MUTATIONAL_SIGNATURE_DATA_DIR, "signatures_probabilities.sorted.txt")
MUTATIONAL_SIGNATURE_INFO_FILE = os.path.join(MUTATIONAL_SIGNATURE_DATA_DIR, "signature_analysis_data.formatted.txt")


VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES = ['lrg_identifier', 'refseq_transcript_accession', 'ensembl_transcript_accession']
# Use highest TranscriptVersion canonical, set False to use representative transcript (ie VEP pick = variant annotation)
VARIANT_TRANSCRIPT_USE_TRANSCRIPT_CANONICAL = True

VARIANT_ZYGOSITY_GLOBAL_COLLECTION = "global"

PREFER_ALLELE_LINKS = False

# ClinGen Allele Registry paper - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6519371/
CLINGEN_ALLELE_REGISTRY_DOMAIN = "http://reg.genome.network"
CLINGEN_ALLELE_REGISTRY_BATCH_SIZE = 1000
CLINGEN_ALLELE_REGISTRY_MAX_RECORDS = 2000
# From paper: The maximal nucleotide (transcript or genomic) allele size is 10,000 bp
CLINGEN_ALLELE_REGISTRY_LOGIN = get_secret("CLINGEN_ALLELE_REGISTRY.login")
CLINGEN_ALLELE_REGISTRY_PASSWORD = get_secret("CLINGEN_ALLELE_REGISTRY.password")
CLINGEN_ALLELE_REGISTRY_MAX_MANUAL_REQUESTS = 10_000  # On nodes and VCFs
CLINGEN_ALLELE_REGISTRY_REQUIRE_REF_ALLELE = True
CLINGEN_ALLELE_REGISTRY_REATTEMPT_WITH_ACTUAL_REF = True
CLINGEN_ALLELE_REGISTRY_MAX_CACHE_DAYS = 180  # Set to None to last forever. 180 = ~6 months

NO_DNA_CONTROL_REGEX = "(^|[^a-zA-Z])NDC([^a-zA-Z]|$)"  # No DNA Control - e.g. _NDC_ or -NDC_

VCF_DOWNLOAD_ADMIN_ONLY = False
VCF_IMPORT_CREATE_COHORT_FROM_MULTISAMPLE_VCFS = True
VCF_IMPORT_NO_DNA_CONTROL_SAMPLE_REGEX = None
VCF_IMPORT_FILE_SPLIT_ROWS = 50000
VCF_IMPORT_SKIP_RECORD_REGEX = {
    "Fusion": "VARTYPE=fusion",
}

VCF_IMPORT_COMMON_FILTER_INFO = "VG_GNOMAD_AF"
VCF_IMPORT_COMMON_FILTERS = {
    # If the filenames don't start with "/" they're relative to ANNOTATION_VEP_BASE_DIR
    "GRCh37": {
        "gnomad_af_filename": "annotation_data/GRCh37/gnomad_GRCh37_af_greater_than_5.contigs.vcf.bgz",
        "gnomad_version": "2.1.1",
        "gnomad_af_min": 0.05,
        "clinical_significance_max": "3",
    },
    "GRCh38": {
        "gnomad_af_filename": "annotation_data/GRCh38/gnomad4.0_GRCh38_af_greater_than_5.stripped.vcf.gz",
        "gnomad_version": "4.0",
        "gnomad_af_min": 0.05,
        "clinical_significance_max": "3",
    },
    "T2T-CHM13v2.0": {
        "gnomad_af_filename": "annotation_data/T2T-CHM13v2.0/gnomad4.1.t2t_liftover_T2T-CHM13v2.0_af_greater_than_5.stripped.vcf.gz",
        "gnomad_version": "4.1.t2t_liftover",
        "gnomad_af_min": 0.05,
        "clinical_significance_max": "3",
    }
}


COMPANY = None  # Used for gene list categories

GENERATED_DIR = os.path.join(MEDIA_ROOT, 'generated')

HGVS_DEFAULT_METHOD = "biocommons_hgvs"  # HGVSConverterType (any case) ie "pyhgvs", "biocommons_hgvs", "combo"
HGVS_MAX_REF_ALLELE_LENGTH = 10  # Set to 0 for "del" instead of "delC" etc
HGVS_VALIDATE_REFSEQ_TRANSCRIPT_LENGTH = True
HGVS_MAX_SEQUENCE_LENGTH = 50_000
HGVS_MAX_SEQUENCE_LENGTH_REPRESENTATIVE_TRANSCRIPT = 200_000
HGVS_RETRIEVE_TRANSCRIPT_SEQUENCE = False  # Biocommons only - attempt to retrieve transcript sequence (slower)

PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB = False
IMPORT_PROCESSING_DIR = os.path.join(PRIVATE_DATA_ROOT, 'import_processing')
IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS = not DEBUG

# @see https://github.com/SACGF/variantgrid/wiki/Liftover
LIFTOVER_CLASSIFICATIONS = True
LIFTOVER_TO_CHROMOSOMES_ONLY = True  # False = Liftover to alt/patches
LIFTOVER_DBSNP_ENABLED = False  # Default=False - doesn't work so well due to dbSNP IDs being for loci

LIFTOVER_BCFTOOLS_ENABLED = True
# 2025-02-13 - bcftools liftover currently gives a warning about symbolic variants
LIFTOVER_BCFTOOLS_SYMBOLIC = False
LIFTOVER_BCFTOOLS_MAX_LENGTH = 1000
LIFTOVER_BCFTOOLS_PLUGIN_DIR = "/usr/share/bcftools/plugins"
BCFTOOLS_COMMAND = "bcftools"  # if not absolute, needs to be in path

PANEL_APP_CACHE_DAYS = 7  # Automatically re-check after this time
GENE_RELATION_PANEL_APP_LIVE_UPDATE = False  # Use GenCC cached result if False, poll panel app if True

# "PanelApp Australia"

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
GENE_GRID_FAKE_GOLD_COVERAGE = False  # For user testing

GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = None
VIEW_GENE_SYMBOL_SHOW_GENE_COVERAGE = False
VIEW_TRANSCRIPT_SHOW_CLASSIFICATIONS = True
VIEW_TRANSCRIPT_VERSION_SHOW_CLASSIFICATIONS = False

DEFAULT_COLUMNS_NAME = 'Default columns'

DEFAULT_ENRICHMENT_KIT_LEFT_PADDING = 0
DEFAULT_ENRICHMENT_KIT_RIGHT_PADDING = 0

ANALYSIS_DUAL_SCREEN_MODE_FEATURE_ENABLED = False  # Currently broken
ANALYSIS_TEMPLATES_AUTO_SAMPLE = "Sample tab auto analysis"
ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT = "Cohort VCF Export auto analysis"
ANALYSIS_WARN_IF_NO_QC_GENE_LIST_MESSAGE = None  # disabled by default
ANALYSIS_NODE_CACHE_Q = True
ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX = 1000

VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT = True  # For analysis Grid/CSV export. VCF export is always unit
VARIANT_SHOW_CANONICAL_HGVS = True

# if True, CR_lab_id will be used in all instances
CLASSIFICATION_ID_OVERRIDE_PREFIX = False
CLASSIFICATION_OMNI_IMPORTER_APP_DIR = None  # location of OmniImporter (if present) to map imported classifications
CLASSIFICATION_OMNI_IMPORTER_DATA_DIR = os.path.join(PRIVATE_DATA_ROOT, "lab_classification_files")  # location of directory for saving and mapping files to send to OmniImporter (can by any directory with write access)
CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL = "lab"  # change this to logged_in_users for prod environments
CLASSIFICATION_OMNI_IMPORTER_INCLUDE_SOURCE = False  # change this to True for dev environments (too dangerous to set to True by default)
CLASSIFICATION_OMNI_IMPORTER_PYTHON_COMMAND = PYTHON_COMMAND
CLASSIFICATION_OMNI_IMPORTER_PARSERS = ["vg_tags", "mvl_patch", "curio"]

CLASSIFICATION_SUPPORTED_TRANSCRIPTS = {"NR", "NM", "NC", "ENST", "LRG_", "XR"}
CLASSIFICATION_MATCH_VARIANTS = True  # exists only so we can turn it off during testing
CLASSIFICATION_REQUIRE_OVERWRITE_NOTE = True
CLASSIFICATION_AUTOFUZZ_AGE = False

CLASSIFICATION_DASHBOARD_SIZE = 50
CLASSIFICATION_RECLASSIFICATION_EMAIL = True
CLASSIFICATION_ID_FILTER = True
CLASSIFICATION_GRID_SHOW_USERNAME = True

CLASSIFICATION_STATS_USE_SHARED = False  # False=Use visible to user. True = Shared
CLASSIFICATION_GRID_SHOW_PHGVS = True
CLASSIFICATION_GRID_SHOW_SAMPLE = True
CLASSIFICATION_GRID_MULTI_LAB_FILTER = False
CLASSIFICATION_SHOW_SPECIMEN_ID = True
CLASSIFICATION_NEW_GROUPING = False

# Require people to click "my sample's not here" (ie encourage them to find it)
CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE = True
CLASSIFICATION_WEB_FORM_CREATE = True
CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN = True
CLASSIFICATION_WEB_FORM_CREATE_ALLOW_NO_VARIANT = True  # Can create purely from HGVS

# If change is synonymous, and molecular consequence is splicing, then change from p.= to p.?
CLASSIFICATION_AUTO_POPULATE_P_HGVS_SYNONYMOUS_SPLICE_CHANGE_TO_UNKNOWN = False

CLASSIFICATION_FILE_ATTACHMENTS = True  # allow users to attach files to classifications

CLASSIFICATION_MAX_REFERENCE_LENGTH = 100  # Used for MVL export, general display use HGVS_MAX_REF_ALLELE_LENGTH

CLASSIFICATION_REDCAP_EXPORT = True
CLASSIFICATION_NON_ACMG_ASSERTION_METHOD = None  # when calculating ACMG points, even if we have ACMG criteria, are they a little too trnaslated to be useful

CLASSIFICATION_ALLOW_DELETE = True
"""
Is a hard-delete offered to classification owners (if false admins will have to delete from admin screen)
"""
CLASSIFICATION_ALLOW_UNKNOWN_KEYS = True  # default to true for the sake of environments syncing from other environments

ONTOLOGY_EXTERNAL_LINKS = False  # Generate external or internal links for ontology terms

PATHOLOGY_TESTS_ENABLED = False
PATHOLOGY_TEST_REQUESTS_REDIRECT_URL = None
PATHOLOGY_TEST_EXTERNAL_CODE = None
PATHOLOGY_TEST_SORTED_ENRICHMENT_KITS = []  # Tests get automatically compared to
PATHOLOGY_TEST_CASE_EXTERNAL_CODE = None

PEDIGREE_MIN_COHORT_SAMPLE_MATCHES_FOR_AUTO_MATCH = 3
PEDIGREE_MADELINE2_COMMAND = None  # Install https://madeline.med.umich.edu/madeline/ set this to "madeline2"

REQUESTS_DISABLE_IPV6 = True  # we've found a performance penalty w/IPv6

INITIAL_USER_DATA_PREFIX_KWARGS = {}  # Create UserDataPrefix object to setup IGV for new users
ISSUE_TRACKER_URL = "https://github.com/SACGF/variantgrid/issues"
ISSUE_TRACKER_TEXT = "GitHub issue tracker"

USER_SETTINGS_SHOW_GROUPS = True

SQL_BATCH_INSERT_SIZE = 50000
SQL_SCRIPTS_DIR = os.path.join(BASE_DIR, "dbscripts")
SITE_NAME = "VariantGrid"

# TODO instead of make settings for admin and non-admin enabled searches, have a search value that can be
# DISABLED, ADMIN_ONLY, ENABLED
SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN = False  # set True to only find classified variants
SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR = True
SEARCH_HGVS_GENE_SYMBOL = True
SEARCH_HGVS_GENE_SYMBOL_USE_MANE = True
SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS = False
SEARCH_HGVS_GENE_SYMBOL_REPORT_FAILURES = True
SEARCH_USER_ADMIN_ONLY = False
SEARCH_COSMIC_ENABLED = True
SEARCH_COSMIC_TRANSCRIPT_MESSAGES = False
SEARCH_CONTIG_GENOME_BUILD_ADMIN_ONLY = False

SILENCED_SYSTEM_CHECKS = [
    'models.E006',  # 'captcha.recaptcha_test_key_error'
]
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
    'enabled': bool(rollbar_access_token and rollbar_client_access_token),  # set to false in environments to disable rollbar
    'branch': 'master',
    'root': BASE_DIR,
    'capture_username': True,
    'code_version': Git(BASE_DIR).hash,
    'ignorable_404_urls': (
        re.compile(r'.*\.map'),
        re.compile(r'.*\.ico')
    ),
}

SLACK = {
    'enabled': get_secret('SLACK.enabled'),
    'admin_callback_url': get_secret('SLACK.admin_callback_url'),
    'emoji': ':dna:'  # overwrite the emoji for different environments
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
    'auditlog.middleware.AuditlogMiddleware',
    #'rollbar.contrib.django.middleware.RollbarNotifierMiddleware',
    # if you want to always avoid 404, use
    # 'rollbar.contrib.django.middleware.RollbarNotifierMiddlewareExcluding404'

    'htmlmin.middleware.HtmlMinifyMiddleware',
    'htmlmin.middleware.MarkRequestMiddleware',
    'threadlocals.middleware.ThreadLocalMiddleware',
    # 'querycount.middleware.QueryCountMiddleware',
    # 'mozilla_django_oidc.middleware.SessionRefresh',
)
HTML_MINIFY = True
EXCLUDE_FROM_MINIFYING = (
    '^media/',  # Died in Django 4: FileResponse instance has no 'content' attribute. Use 'streaming_content' instead.
)
EXCLUDE_TAGS_FROM_MINIFYING = ("pre", "script", "textarea", "nomin", "code")

ROOT_URLCONF = 'variantgrid.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'variantgrid.wsgi.application'

# We want a class to be able to register itself as a factory
# The following will be imported, and looked for ImportTaskFactory.__subclasses__()
IMPORT_TASK_FACTORY_IMPORTS = (
    'annotation.import_task_factories',
    'upload.import_task_factories.import_task_factories',
)

# To be able to automatically do things after VCF import - e.g. damage counts
FINISH_IMPORT_VCF_STEP_TASKS_CLASSES = []

# Turn ON in production!
CACHE_GENERATED_FILES = True

REST_FRAMEWORK = {
    # NOTE: Middleware is run first - so GlobalLoginRequiredMiddleware will reject tokens w/o logins
    # before DRF even sees it. You need to add your APIs to PUBLIC_PATHS
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'rest_framework.authentication.BasicAuthentication',  # Needed to classification export clients
        'rest_framework.authentication.SessionAuthentication',
        'rest_framework.authentication.TokenAuthentication',
    ],
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
LOGO_VIEW_NAME = "index"
LOGIN_USERNAME_PLACEHOLDER = None
LOGOUT_REDIRECT_URL = '/'

ACCOUNT_ACTIVATION_DAYS = 7  # One-week activation window

INSTALLED_APPS = [
    'compressor',
    'avatar',
    'vcauth',
    'datetimeutc',
    'django_json_widget',
    'django.contrib.humanize',
    'django.contrib.sites',
    'django_recaptcha',
    'registration',
    'django.contrib.auth',
    'mozilla_django_oidc',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    "django.contrib.postgres",
    # 3rd party libraries
    'auditlog',
    'dal',  # Django Autocomplete Light v3
    'dal_select2',  # DAL Plugin
    'django_messages',
    'django_dag',
    'django_js_reverse',
    'django_starfield',
    'django_extensions',
    'djgeojson',
    'easy_thumbnails',
    'fontawesomefree',
    'guardian',
    'jfu',
    'leaflet',
    "psqlextra",
    'rest_framework',
    'rest_framework.authtoken',
    'termsandconditions',
    'crispy_forms',  # used to make bootstrap compatible forms
    'crispy_bootstrap4',
    # Internal apps
    'analysis.apps.AnalysisConfig',
    'annotation.apps.AnnotationConfig',
    'ontology',
    'eventlog',
    'expression',  # Empty but referenced by other migrations - can delete after squash
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
    'review'
    # Uncomment the next line to enable admin documentation:
    # 'django.contrib.admindocs',
]

# https://django-crispy-forms.readthedocs.io/en/latest/install.html
# CRISPY_ALLOWED_TEMPLATE_PACKS = ('bootstrap4_neat', 'bootstrap4') # need to do this if you make an alternative template pack
CRISPY_TEMPLATE_PACK = 'bootstrap4'
CRISPY_CLASS_CONVERTERS = {
    "blanknullbooleanselect": "form-control",
    "select": "form-control",
    "numberinput": "form-control",
    "textinput": "form-control"
}

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
        # 'file': {
        #     'level': 'DEBUG',
        #     'class': 'logging.FileHandler',
        #     'filename': '/tmp/django_debug.log',
        # },
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
            'propagate': False,
            'level': 'INFO',
        },
        'hgvs': {
            'level': 'WARNING',
        },
        #        'django.request': {
        #            'handlers': ['mail_admins'],
        #            'level': 'ERROR',
        #            'propagate': False,
        #        },
        'django.db.backends': {
            'handlers': ['console'],
            'propagate': False,
            # 'level': 'DEBUG',
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
    r'^/classification/api/.*',  # REST framework used by command line tools
    r'^/seqauto/api/.*',
    r'^/upload/api/.*',
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
            "T2T-CHM13v2.0": "sites.chm13v2.T2T.vcf.gz",
        },
    },
    # Minimums for related samples to appear at bottom of view_sample page
    "relatedness": {
        "min_relatedness": 0.1,
        "min_shared_hets": 1000,
        "min_shared_hom_alts": 200,
    }
}


# @see https://github.com/SACGF/variantgrid/wiki/URL---Menu-configuration
# Before URLs are registered, the URLS_APP_REGISTER and URLS_NAME_REGISTER are looked up
# To make a whitelist - change the default to False, then add overrides, e.g. 'url_name' : True for allowed
# To make a blacklist - leave default as True, add 'url_name' : False

# Use this to handle url registration at the app level (to e.g. block and entire app)
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
    "hgvs_resolution_tool": False,
    "classification_view_metrics": False,
    "clinvar_key_summary": False,
    "clinvar_match": False,
    "condition_matchings": False,
    "condition_match_test": False,
    "discordance_reports": False,
    "vus": False
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
                           "classification_upload_unmapped": False})

VARIANT_DETAILS_SHOW_ANNOTATION = True  # also doubles as GENE_SHOW_ANNOTATION
VARIANT_DETAILS_SHOW_GENE_COVERAGE = False
VARIANT_DETAILS_SHOW_SAMPLES = True
VARIANT_DETAILS_NEARBY_RANGE = 50
VARIANT_DETAILS_NEARBY_SHOW_GENE = False

# Can only use block or allow list. Use "text" value of link (see vc_links.js)
VARIANT_DETAILS_QUICK_LINKS_ALLOW_LIST = []
VARIANT_DETAILS_QUICK_LINKS_BLOCK_LIST = []

VARIANT_VCF_DB_PREFIX = "vg"
VARIANT_MANUAL_CREATE = True
VARIANT_MANUAL_CREATE_BY_NON_ADMIN = True

# Below this size, variants are stored with ref/alt sequences. Above this threshold, they become
# structural variants and use symbolic
VARIANT_SYMBOLIC_ALT_ENABLED = True
VARIANT_SYMBOLIC_ALT_SIZE = 1000
VARIANT_SYMBOLIC_ALT_VALID_TYPES = {VCFSymbolicAllele.CNV, VCFSymbolicAllele.DEL,
                                    VCFSymbolicAllele.DUP, VCFSymbolicAllele.INV}

# VCF 4.3 defines SVLEN as negative for DELs, but 4.4+ it's always positive
# This is for our internal storage in Variant records
# See https://github.com/SACGF/variantgrid/issues/1344
VARIANT_SYMBOLIC_ALT_SVLEN_ALWAYS_POSITIVE = False

VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS = False
VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS_PREFER_CANONICAL_WITH_DIFF_VERSION = True
VIEW_GENE_HOTSPOT_GRAPH = True
VIEW_GENE_HOTSPOT_GRAPH_PREFER_CANONICAL_WITH_DIFF_VERSION = True
VIEW_GENE_WIKI = True

KEY_CLOAK_REALM = None

UNSHARED_FLAG_ENABLED = False

## COMPRESSION SETTINGS
COMPRESS_ENABLED = True

# Standard SECRETS
# Config used to talk directly to keycloak
KEYCLOAK_SYNC_DETAILS = None
# Config used to upload to Shariant, keys are SyncDestination.config['sync_details']
SYNC_DETAILS = {}
# Config used if using AWS email
AWS_SES_ACCESS_KEY_ID = None
AWS_SES_SECRET_ACCESS_KEY = None
AWS_SES_REGION = None

CLINVAR_EXPORT = None


# Command line tool to unzip and cat a file to stdout
# for macOS need to set this to gzcat as the default zcat has short-comings
BASH_ZCAT = 'zcat'
# If True, will run a series of bash commands as one long string with Shell=True
# otherwise will pipe each command into the next more safely with Shell=False
VCF_IMPORT_PREPROCESS_POPEN_SHELL = True  # For vcf split
VCF_EXPORT_VERSION = "4.3"

CLASSIFICATION_DOWNLOADABLE_JSON_LITERATURE_CITATIONS = False
CLASSIFICATION_DOWNLOADABLE_NOTES_AND_EXPLAINS = True
CLASSIFICATION_DOWNLOADABLE_FIELDS = "*"

# Bootstrapped themed messages
from django.contrib.messages import constants as messages

# @see https://docs.djangoproject.com/en/3.1/ref/settings/#message-tags
# Provide both Bootstrap classes (for within the product) and default (for within the admin)
MESSAGE_TAGS = {
    messages.DEBUG: 'debug alert-info',
    messages.INFO: 'info alert-info',  # annoyingly this is branded same as success in admin
    messages.SUCCESS: 'success alert-success',
    messages.WARNING: 'warning alert-warning',
    messages.ERROR: 'error alert-danger',
}

def get_clinvar_export_secrets() -> dict:
    return get_secrets("CLINVAR_EXPORT", ["mode", "api_key", "org_id"])


def get_keycloak_sync_secrets() -> dict:
    """
    KEYCLOAK_SYNC_DETAILS = get_keycloak_sync_secrets()
    """
    return get_secrets("KEYCLOAK", ["username", "password", "host", "oauth_url", "client_id"])


def get_shariant_sync_secrets() -> dict:
    """
    SYNC_DETAILS = get_shariant_sync_secrets()
    """
    sync_fields = ["enabled", "username", "password", "host", "oauth_url", "client_id"]
    sync = get_secret("SYNC")
    # New SYNC keys are SyncDestination.config['sync_details'] - old one was only 1 level deep
    # Check whether it's old sync field
    if set(sync.keys()) == set(sync_fields):
        raise ValueError("Old secret 'SYNC' detected - need to use SyncDestination.config['sync_details'] as keys")

    sync_all_fields = ["enabled", "username", "password", "host", "oauth_url", "client_id", "app_username", "app_password"]
    return {sd: get_secrets(f"SYNC.{sd}", sync_all_fields, False) for sd in sync}


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
