"""
Shariant - https://shariant.org.au

See https://github.com/sacgf/variantgrid/wiki/Annotation%20Setup

"""

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
import re

# for nginx providing HTTPS
USE_X_FORWARDED_HOST = True
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

SYNC_DETAILS = get_shariant_sync_secrets()

ALLELE_VALIDATION = True
AVATAR_PROVIDERS = (
    'library.django_utils.avatar.SpaceThemedAvatarProvider'
)
CLASSIFICATION_REDCAP_EXPORT = False

CLINVAR_EXPORT = get_clinvar_export_secrets()

aws_dict = get_aws_secrets()
AWS_SES_ACCESS_KEY_ID, AWS_SES_SECRET_ACCESS_KEY, AWS_SES_REGION = \
        [aws_dict[k] for k in ('AWS_SES_ACCESS_KEY_ID', 'AWS_SES_SECRET_ACCESS_KEY', 'AWS_SES_REGION')]

aws_s3_dict = get_s3_secrets()
AWS_S3_ACCESS_KEY_ID, AWS_S3_SECRET_ACCESS_KEY = \
    [aws_s3_dict[k] for k in ("AWS_S3_ACCESS_KEY_ID", "AWS_S3_SECRET_ACCESS_KEY")]
AWS_DEFAULT_ACL = None

KEYCLOAK_SYNC_DETAILS = get_keycloak_sync_secrets()

# HEARTBEAT_URL = 'https://heartbeat.uptimerobot.com/m788641874-4c58c98a716180f36670e551a0bd03fff47abfea'

INBOX_ENABLED = False
DISCORDANCE_ENABLED = True
DISCORDANCE_EMAIL = '"Shariant Updates" <updates@shariant.org.au>'
ACCOUNTS_EMAIL = 'accounts@shariant.org.au'
SEND_EMAILS = False
HEALTH_CHECK_ENABLED = True

CLASSIFICATION_GRID_MULTI_LAB_FILTER = False
CLASSIFICATION_NEW_GROUPING = True
CLASSIFICATION_ALLOW_DELETE = False
CLASSIFICATION_NON_ACMG_ASSERTION_METHOD = [
    re.compile(r'.*VCGS.*', flags=re.IGNORECASE),
    re.compile(r'.*Sherloc.*', flags=re.IGNORECASE)
]

# Keycloak

AUTHENTICATION_BACKENDS = (
    'auth.backend.VariantGridOIDCAuthenticationBackend',
    'django.contrib.auth.backends.ModelBackend',  # default
    'guardian.backends.ObjectPermissionBackend',
)

MIDDLEWARE += (
    'auth.session_refresh.VariantGridSessionRefresh',
    'auth.oidc_error_handler.HandleOIDC400Middleware',
)

REST_FRAMEWORK = {
    'DEFAULT_PERMISSION_CLASSES': (
        'rest_framework.permissions.IsAuthenticated',
    ),
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'mozilla_django_oidc.contrib.drf.OIDCAuthentication',
        'rest_framework.authentication.SessionAuthentication'
    ],
    'EXCEPTION_HANDLER': 'rollbar.contrib.django_rest_framework.post_exception_handler'
}

# user logging
MIDDLEWARE += ('eventlog.middleware.PageViewsMiddleware', )
LOG_ACTIVITY_APPS = {"classification", "variantopedia", "snpdb", "genes", "ontology"}

#OIDC_DRF_AUTH_BACKEND = 'auth.backend.VariantGridOIDCAuthenticationBackend'

# OIDC SETTINGS
USE_OIDC = True
LOGIN_URL = '/oidc_login/'
OIDC_RP_SIGN_ALGO = 'RS256'
# TODO, change this then move it to
OIDC_RP_CLIENT_SECRET = get_secret('OIDC.client_secret')
KEY_CLOAK_BASE = 'https://shariant.org.au/auth'
KEY_CLOAK_REALM = 'agha'
KEY_CLOAK_PROTOCOL_BASE = KEY_CLOAK_BASE + '/realms/' + KEY_CLOAK_REALM + '/protocol/openid-connect'
OIDC_OP_JWKS_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/certs'
OIDC_OP_AUTHORIZATION_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/auth'
OIDC_OP_TOKEN_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/token'
OIDC_OP_USER_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/userinfo'
OIDC_USER_SERVICES = KEY_CLOAK_BASE + '/realms/' + KEY_CLOAK_REALM + '/account'
OIDC_OP_LOGOUT_URL_METHOD = 'auth.backend.provider_logout'

# login failure is generally user is inactive, which is how prod distinguishes between prod and test logins

HELP_URL = "https://shariant.readthedocs.io/en/latest/"
# LOGIN_REDIRECT_URL = '/variantopedia/dashboard'
LOGIN_REDIRECT_URL_FAILURE = '/accounts/logout/'

EMAIL_BACKEND = 'django_amazon_ses.EmailBackend'
# Overwrite settings for your system below

DEBUG = False
ALLOWED_HOSTS = ['*']
CSRF_TRUSTED_ORIGINS=['https://test.shariant.org.au', 'https://shariant.org.au', 'https://www.shariant.org.au', 'https://demo.shariant.org.au']

ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED = False  # Only used for analysis optimisation
_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'annotation_scratch')
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

ANNOTATION[BUILD_GRCH37].update({
    "annotation_consortium": "RefSeq",
})
ANNOTATION[BUILD_GRCH38].update({
    "enabled": True,
    "annotation_consortium": "RefSeq",
})

LOGIN_REDIRECT_URL = '/classification/dashboard'
LOGO_VIEW_NAME = "classification_dashboard"

SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
STATICFILES_DIRS = (SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)

HGVS_DEFAULT_METHOD = "pyhgvs"

# SITE_MESSAGE = "Shariant is currently in pre-BETA. Please excuse bugs and missing features, and the site may be shut down for upgrades"

# "LRG_" has been disabled, see https://github.com/SACGF/shariant-admin/issues/126
CLASSIFICATION_OMNI_IMPORTER_APP_DIR = "/opt/shariant_omni_importer"
CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL = "logged_in_users"
CLASSIFICATION_OMNI_IMPORTER_INCLUDE_SOURCE = False

CLASSIFICATION_SUPPORTED_TRANSCRIPTS = {"NM", "NR", "ENST", "XR"}
CLASSIFICATION_REQUIRE_OVERWRITE_NOTE = False
CLASSIFICATION_AUTOFUZZ_AGE = True
CLASSIFICATION_GRID_SHOW_USERNAME = False  # In Shariant - this may be a lab's API user so hide it
CLASSIFICATION_STATS_USE_SHARED = True  # False=Use visible to user. True = Shared
CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE = False
CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN = False
CLASSIFICATION_GRID_SHOW_PHGVS = False  # would be nice to show but isn't populated consistently
CLASSIFICATION_GRID_SHOW_SAMPLE = False
CLASSIFICATION_FILE_ATTACHMENTS = False
CLASSIFICATION_ID_FILTER = False
CLASSIFICATION_SHOW_SPECIMEN_ID = False

VARIANT_SHOW_CANONICAL_HGVS = False
CLASSIFICATION_MAX_REFERENCE_LENGTH = 1000000  # Try to generate large values for sake of MVL

# Lock down Shariant menu - hide a lot of VariantGrid urls
# Completely hide URLS from these apps
URLS_APP_REGISTER.update({"analysis": False,
                          "expression": False,
                          "pathtests": False,
                          "pedigree": False,
                          "seqauto": False})

URLS_NAME_REGISTER.update({  # Disable selected snpdb urls
    "data": False,

    # Selected patient urls
    "cohorts": False,
    "patient_record_imports": False,
    "patient_term_approvals": False,
    "patients": False,

    "gene_lists": False,
    "genes": False,
    "gene_grid": False,

    # Variants
    "variants": False,
    "variant_tags": False,
    "manual_variant_entry": False,
    "variantopedia_wiki": False,

    # Annotation
    "version_diffs": False,

    # Settings
    "change_password": False,
    "custom_columns": False,
    "tag_settings": False,
    "igv_integration": False,
    "sequencing_software_versions": False,
    "classification_dashboard": True,
    "classification_import_tool": True,
    "hgvs_resolution_tool": True,
    "keycloak_admin": True,
    "hgvs_issues": True,

    # Upload - list all URLS (only want them visible by admin)
    "upload": False,
    "upload_poll": False,
    "view_uploaded_file": False,
    "view_upload_pipeline": False,
    "view_upload_pipeline_warnings_and_errors": False,
    "upload_retry_import": False,
    "upload_pipeline_modified_variants_grid": False,
    "view_upload_stats_detail": False,
    "accept_vcf_import_info_tag": False,
    "jfu_upload": False,
    "jfu_delete": False,
    "download_uploaded_file": False,

    # discordance
    "discordance_reports": True,

    "classification_upload_unmapped": True,
    "condition_matchings": True,
    "condition_match_test": True,
    "classification_view_metrics": True,

    # ClinVarExport
    "clinvar_key_summary": True,
    "clinvar_match": True,
    "vus": True
})

PREFER_ALLELE_LINKS = True

VARIANT_MANUAL_CREATE_BY_NON_ADMIN = False
# Don't show annotation or samples on variant page - don't want to be responsible for it
VARIANT_DETAILS_SHOW_ANNOTATION = False
VARIANT_DETAILS_SHOW_SAMPLES = False
VARIANT_VCF_DB_PREFIX = "stv"
VARIANT_SYMBOLIC_ALT_ENABLED = False

VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS = False
VIEW_GENE_WIKI = False

USER_SETTINGS_SHOW_GROUPS = False
USER_SETTINGS_SHOW_BUILDS = False

SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN = True
SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR = False
SEARCH_CLASSIFICATION_HGVS_SUFFIX = False
SEARCH_HGVS_GENE_SYMBOL_USE_MANE = False
SEARCH_HGVS_GENE_SYMBOL = False
SEARCH_USER_ADMIN_ONLY = True
SEARCH_CONTIG_GENOME_BUILD_ADMIN_ONLY = True

UNSHARED_FLAG_ENABLED = True
VIEW_GENE_HOTSPOT_GRAPH = False

GENE_RELATION_PANEL_APP_LIVE_UPDATE = True

# if False lab record ID will be used for records instead of CR_lab_id
CLASSIFICATION_ID_OVERRIDE_PREFIX = True
ALLELE_ORIGIN_NOT_PROVIDED_BUCKET = "G"
