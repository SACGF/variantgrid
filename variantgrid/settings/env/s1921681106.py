"""
See https://bitbucket.org/sacgf/variantgrid/wiki/Annotation%20Setup
(test)
"""

import json

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
"""
EMAIL_BACKEND = 'django_amazon_ses.EmailBackend'
aws_dict = get_aws_secrets()
AWS_SES_ACCESS_KEY_ID, AWS_SES_SECRET_ACCESS_KEY, AWS_SES_REGION = \
        [aws_dict[k] for k in ('AWS_SES_ACCESS_KEY_ID', 'AWS_SES_SECRET_ACCESS_KEY', 'AWS_SES_REGION')]

KEYCLOAK_SYNC_DETAILS = get_keycloak_sync_secrets()
"""
aws_s3_dict = get_s3_secrets()
AWS_S3_ACCESS_KEY_ID, AWS_S3_SECRET_ACCESS_KEY = \
    [aws_s3_dict[k] for k in ("AWS_S3_ACCESS_KEY_ID", "AWS_S3_SECRET_ACCESS_KEY")]
AWS_DEFAULT_ACL = None

CLINVAR_EXPORT = get_clinvar_export_secrets()

SYNC_DETAILS = get_shariant_sync_secrets()
VARIANT_CLASSIFICATION_NEW_GROUPING = True
COMPRESS_ENABLED = False
BASH_ZCAT = 'gzcat'
POPEN_SHELL = True

USER_SETTINGS_SHOW_GROUPS = False

# import all the base settings #
CELERY_ENABLED = False
DISCORDANCE_ENABLED = True
VARIANT_CLASSIFICATION_STATS_USE_SHARED = True
DISCORDANCE_EMAIL = 'discordance@test.shariant.org.au'
ACCOUNTS_EMAIL = 'accounts@shariant.org.au'

ROLLBAR['enabled'] = False

INTERNAL_IPS = [
    '127.0.0.1',
    '0.0.0.0',
    'localhost'
]

INSTALLED_APPS += [
    # 'debug_toolbar'
]

MIDDLEWARE += (
    # 'debug_toolbar.middleware.DebugToolbarMiddleware',
)

"""
AUTHENTICATION_BACKENDS = (
    'auth.backend.VariantGridOIDCAuthenticationBackend',
    'django.contrib.auth.backends.ModelBackend', # default
    'guardian.backends.ObjectPermissionBackend',
)

MIDDLEWARE += (
    'mozilla_django_oidc.middleware.SessionRefresh',
    'debug_toolbar.middleware.DebugToolbarMiddleware',
    'auth.oidc_error_handler.HandleOIDC400Middleware'
)

REST_FRAMEWORK = {
    'DEFAULT_PERMISSION_CLASSES': (
        'rest_framework.permissions.IsAuthenticated',
    ),
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'mozilla_django_oidc.contrib.drf.OIDCAuthentication',
        'rest_framework.authentication.SessionAuthentication'
    ],
}

OIDC_DRF_AUTH_BACKEND = 'auth.backend.VariantGridOIDCAuthenticationBackend'
USE_OIDC = True
#OIDC_REQUIRED_GROUP = '/variantgrid/shariant_production'
LOGIN_URL = '/oidc_login/'

OIDC_RP_SIGN_ALGO = 'RS256'
OIDC_RP_CLIENT_ID = 'variantgrid'
OIDC_RP_CLIENT_SECRET = get_secret('OIDC.client_secret')
KEY_CLOAK_BASE = 'http://localhost:8080/auth'
#KEY_CLOAK_BASE = 'http://localhost:8080/auth'

KEY_CLOAK_REALM = 'agha'
KEY_CLOAK_PROTOCOL_BASE = KEY_CLOAK_BASE + '/realms/' + KEY_CLOAK_REALM + '/protocol/openid-connect'
OIDC_OP_JWKS_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/certs'
OIDC_OP_AUTHORIZATION_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/auth'
OIDC_OP_TOKEN_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/token'
OIDC_OP_USER_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/userinfo'
OIDC_USER_SERVICES = KEY_CLOAK_BASE + '/realms/' + KEY_CLOAK_REALM + '/account'
OIDC_OP_LOGOUT_URL_METHOD = 'auth.backend.provider_logout'
LOGIN_REDIRECT_URL = '/variantopedia/dashboard'
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=http%3A%2F%2Fubuntu-18.04%3A8000'
"""

VARIANT_CLASSIFICATION_AUTOFUZZ_AGE = True
VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME = True

# Use legacy until we generate "global"
# VARIANT_ZYGOSITY_GLOBAL_COLLECTION = "legacy"

# Overwrite settings for your system below
ALLOWED_HOSTS = ["*"]
"""
if DB_FILE == 'local_login':
    INTERNAL_IPS = [
        '127.0.0.1',
        '10.211.55.2'
    ]
"""
ANNOTATION_ENTREZ_EMAIL = 'James.Andrews@sa.gov.au'

ANNOTATION_BASE_DIR = "/Users/jamesandrews/Projects/VariantGrid/data/annotation"
ANNOTATION_VCF_DUMP_DIR = os.path.join(ANNOTATION_BASE_DIR, 'annotation_scratch')
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "ensembl-vep")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")
_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "fasta")

ANNOTATION = {
    BUILD_GRCH37: {
        "enabled": True,
        "annotation_consortium": "Ensembl",
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"),
        "reference_fasta_has_chr": False,
        "cytoband": os.path.join(VG_REFERENCE_DIR, "hg19", "cytoband.hg19.txt.gz"),

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "dbnsfp": "annotation_data/GRCh37/dbNSFP4.0b2a.hg19.stripped.gz",
            "dbscsnv": "annotation_data/GRCh37/dbscSNV1.1_GRCh37.txt.gz",
            "gnomad": "annotation_data/GRCh37/gnomad_GRCh37_combined_af.vcf.bgz",
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh37/hg19.100way.phastCons.bw",
            'phastcons30way': None,  # n/a for GRCh37
            'phastcons46way': "annotation_data/GRCh37/hg19.phastCons46way.placental.bw",
            'phylop100way': "annotation_data/GRCh37/hg19.100way.phyloP100way.bw",
            'phylop30way': None,  # n/a for GRCh37
            'phylop46way': "annotation_data/GRCh37/hg19.phyloP46way.placental.bw",
            "repeatmasker": "annotation_data/GRCh37/repeatmasker_hg19.bed.gz",
            "topmed": "annotation_data/GRCh37/TOPMED_GRCh37.vcf.gz",
            "uk10k": "annotation_data/GRCh37/UK10K_COHORT.20160215.sites.vcf.gz",
        }
    },
    # GRCh38 is NOT enabled by default - overwrite "enabled" in your server settings to use
    BUILD_GRCH38: {
        "enabled": False,
        "annotation_consortium": "Ensembl",
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
        "reference_fasta_has_chr": False,
        "cytoband": os.path.join(VG_REFERENCE_DIR, "hg38", "cytoband.hg38.txt.gz"),

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "dbnsfp": "annotation_data/GRCh38/dbNSFP4.0b2a.hg38.stripped.gz",
            "dbscsnv": "annotation_data/GRCh38/dbscSNV1.1_GRCh38.txt.gz",
            "gnomad": "annotation_data/GRCh38/gnomad_GRCh38_combined_af.vcf.bgz",
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh38/hg38.phastCons100way.bw",
            'phastcons30way': "annotation_data/GRCh38/hg38.phastCons30way.bw",
            'phylop100way': "annotation_data/GRCh38/hg38.phyloP100way.bw",
            'phylop30way': "annotation_data/GRCh38/hg38.phyloP30way.bw",
            "repeatmasker": "annotation_data/GRCh38/repeatmasker_hg38.bed.gz",
            "topmed": "annotation_data/GRCh38/TOPMED_GRCh38_20180418.vcf.gz",
            "uk10k": "annotation_data/GRCh38/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz",
        }
    },
}

VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE = False

USER_SETTINGS_SHOW_BUILDS = False

_SHARIANT_MODE = True
if _SHARIANT_MODE:
    LOGIN_REDIRECT_URL = '/classification/dashboard'

    MIMIC = "shariant"
    SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", f"{MIMIC}_static")
    SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_test_static")
    STATICFILES_DIRS = (SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

    SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, f"templates/{MIMIC}_templates")
    TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
    SITE_NAME = MIMIC

SITE_NAME = "James Dev"

#INSTALLED_APPS.append('debug_toolbar')

#SAPATH_APP = 'sapath.apps.SapathConfig'
#INSTALLED_APPS += [SAPATH_APP]

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'root': {
        'level': 'DEBUG',
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
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },
        'db':{
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
            'level': 'WARNING',
        },
#        'django.request': {
#            'handlers': ['mail_admins'],
#            'level': 'ERROR',
#            'propagate': False,
#        },
        'snpdb': {
            'handlers': ['console'],
            'level': 'DEBUG',
        },
        'django.db.backends': {
            'handlers': ['console'],
            'propagate': False,
#            'level':'DEBUG',
        },
    }
}
"""
URLS_APP_REGISTER.update({"analysis" : False,
                          "expression" : False,
                          "pathtests" : False,
                          "pedigree" : False,
                          "seqauto" : False,
                          "upload" : False})

URLS_NAME_REGISTER.update({ # Disable selected snpdb urls
    "data" : False,
    "manual_variant_entry" : True,

    # Selected patient urls
    "cohorts" : False,
    "patient_record_imports" : False,
    "patient_term_approvals" : False,
    "patients": False,

    # Disabled selected gene urls (still need some)
    "gene_lists" : False,
    "genes" : True,
    "gene_grid" : False,

    # Annotation
    "variants": False,
    "version_diffs" : False,

    # Settings
    "change_password" : False,
    "columns" : False,
    "tag_settings" : False,
    "igv_integration" : False,
    "sequencing_software_versions" : False,
    "classification_dashboard" : True,

    "classification_import_tool": True,
    "keycloak_admin" : True
})
"""
URLS_NAME_REGISTER.update({  # Disable selected snpdb urls
    "classification_dashboard": True,
    "hgvs_issues": True,
    "classification_import_tool": True,
    "keycloak_admin": True,
    "classification_import_upload": True,
    "condition_aliases": True,
    "clinvar_key_summary": True,
    "condition_matchings": True,
    "condition_match_test": True
})

# mimic shariant
VARIANT_DETAILS_SHOW_ANNOTATION = True
VARIANT_DETAILS_SHOW_SAMPLES = True
VARIANT_CLASSIFICATION_GRID_SHOW_PHGVS = True

UNSHARED_FLAG_ENABLED = True
os.environ['OAUTHLIB_INSECURE_TRANSPORT'] = 'true'

VARIANT_CLASSIFICATION_ID_FILTER = True
VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH = True

ANNOTATION_CACHED_WEB_RESOURCES = []
SITE_MESSAGE = "This is the setting site message"
VARIANT_CLASSIFICATION_MATCH_VARIANTS = False

PREFER_ALLELE_LINKS = True
PANEL_APP_CHECK_ENABLED = True

CLINVAR_EXPORT_ENABLED = True