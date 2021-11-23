"""
Shariant - https://shariant.org.au

See https://github.com/sacgf/variantgrid/wiki/Annotation%20Setup

"""

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

AVATAR_PROVIDERS = (
    'library.django_utils.avatar.SpaceThemedAvatarProvider'
)
VARIANT_CLASSIFICATION_REDCAP_EXPORT = False

CLINVAR_EXPORT = get_clinvar_export_secrets()

aws_dict = get_aws_secrets()
AWS_SES_ACCESS_KEY_ID, AWS_SES_SECRET_ACCESS_KEY, AWS_SES_REGION = \
        [aws_dict[k] for k in ('AWS_SES_ACCESS_KEY_ID', 'AWS_SES_SECRET_ACCESS_KEY', 'AWS_SES_REGION')]

aws_s3_dict = get_s3_secrets()
AWS_S3_ACCESS_KEY_ID, AWS_S3_SECRET_ACCESS_KEY = \
    [aws_s3_dict[k] for k in ("AWS_S3_ACCESS_KEY_ID", "AWS_S3_SECRET_ACCESS_KEY")]
AWS_DEFAULT_ACL = None

KEYCLOAK_SYNC_DETAILS = get_keycloak_sync_secrets()

HEARTBEAT_URL = 'https://heartbeat.uptimerobot.com/m788641874-4c58c98a716180f36670e551a0bd03fff47abfea'

# import all the base settings #
SITE_ID = 5  # shariant.org.au

ALLELE_VALIDATION = True
INBOX_ENABLED = False
DISCORDANCE_ENABLED = True
DISCORDANCE_EMAIL = '"Shariant Updates" <updates@shariant.org.au>'
ACCOUNTS_EMAIL = 'accounts@shariant.org.au'
SEND_EMAILS = True
VARIANT_CLASSIFICATION_NEW_GROUPING = True

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

#OIDC_DRF_AUTH_BACKEND = 'auth.backend.VariantGridOIDCAuthenticationBackend'

# OIDC SETTINGS
USE_OIDC = True
OIDC_REQUIRED_GROUP = '/variantgrid/shariant_production'
LOGIN_URL = '/oidc_login/'

OIDC_RP_SIGN_ALGO = 'RS256'
OIDC_RP_CLIENT_ID = 'shariant'
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
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=https%3A%2F%2Fshariant.org.au'
LOGIN_REDIRECT_URL_FAILURE = '/accounts/logout'

EMAIL_BACKEND = 'django_amazon_ses.EmailBackend'
# Overwrite settings for your system below

DEBUG = False
ALLOWED_HOSTS = ['*']

_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'annotation_scratch')
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

ANNOTATION[BUILD_GRCH37]["annotation_consortium"] = "RefSeq"
ANNOTATION[BUILD_GRCH38]["enabled"] = True
ANNOTATION[BUILD_GRCH38]["annotation_consortium"] = "RefSeq"

LIFTOVER_NCBI_REMAP_ENABLED = True
LIFTOVER_NCBI_REMAP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

LOGIN_REDIRECT_URL = '/classification/dashboard'

SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
STATICFILES_DIRS = (SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
SITE_NAME = "Shariant"

# SITE_MESSAGE = "Shariant is currently in pre-BETA. Please excuse bugs and missing features, and the site may be shut down for upgrades"

# "LRG_" has been disabled, see https://github.com/SACGF/shariant-admin/issues/126
VARIANT_CLASSIFICATION_SUPPORTED_TRANSCRIPTS = {"NM", "NR", "ENST", "XR"}
VARIANT_CLASSIFICATION_REQUIRE_OVERWRITE_NOTE = False
VARIANT_CLASSIFICATION_AUTOFUZZ_AGE = True
VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME = False  # In Shariant - this may be a lab's API user so hide it
VARIANT_CLASSIFICATION_STATS_USE_SHARED = True  # False=Use visible to user. True = Shared
VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE = False
VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN = False
VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN = False
VARIANT_CLASSIFICATION_GRID_SHOW_PHGVS = False  # would be nice to show but isn't populated consistently
VARIANT_CLASSIFICATION_FILE_ATTACHMENTS = False
VARIANT_CLASSIFICATION_ID_FILTER = False
VARIANT_CLASSIFICAITON_SHOW_SPECIMEN_ID = False

VARIANT_SHOW_CANONICAL_HGVS = False
VARIANT_CLASSIFICATION_MAX_FULL_ALLELE_LENGTH = 1000000  # Try to generate large values for sake of MVL

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
    "keycloak_admin": True,
    "hgvs_issues": True,

    # Upload - list all URLS (only want them visible by admin)
    "upload": False,
    "upload_poll": False,
    "view_uploaded_file": False,
    "view_upload_pipeline": False,
    "view_upload_pipeline_warnings_and_errors": False,
    "upload_retry_import": False,
    "upload_step_grid": False,
    "upload_pipeline_modified_variants_grid": False,
    "view_upload_stats_detail": False,
    "accept_vcf_import_info_tag": False,
    "jfu_upload": False,
    "jfu_delete": False,
    "download_uploaded_file": False,

    "classification_import_upload": True,
    "condition_matchings": True,
    "condition_match_test": True,
    # "condition_aliases": True

    # ClinVarExport
    "clinvar_key_summary": True
})

PREFER_ALLELE_LINKS = True

VARIANT_MANUAL_CREATE_BY_NON_ADMIN = False
# Don't show annotation or samples on variant page - don't want to be responsible for it
VARIANT_DETAILS_SHOW_ANNOTATION = False
VARIANT_DETAILS_SHOW_SAMPLES = False
VARIANT_VCF_DB_PREFIX = "stv"

VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH = False
VIEW_GENE_SHOW_WIKI = False

USER_SETTINGS_SHOW_GROUPS = False
USER_SETTINGS_SHOW_BUILDS = False

SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN = True
SEARCH_SUMMARY = True  # Little tips that show up on search
SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR = False
SEARCH_SUMMARY_VARIANT_SHOW_CLASSIFICATIONS = False

UNSHARED_FLAG_ENABLED = True
VIEW_GENE_SHOW_HOTSPOT_GRAPH = False

PANEL_APP_CHECK_ENABLED = True