"""
See https://bitbucket.org/sacgf/variantgrid/wiki/Annotation%20Setup
(test)
"""

import json

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

aws_dict = get_aws_secrets()
AWS_SES_ACCESS_KEY_ID, AWS_SES_SECRET_ACCESS_KEY, AWS_SES_REGION = \
        [aws_dict[k] for k in ('AWS_SES_ACCESS_KEY_ID', 'AWS_SES_SECRET_ACCESS_KEY', 'AWS_SES_REGION')]

SYNC_DETAILS = get_shariant_sync_secrets()

KEYCLOAK_SYNC_DETAILS = get_keycloak_sync_secrets()


# import all the base settings #
DISCORDANCE_ENABLED = True
DISCORDANCE_EMAIL = None  # 'discordance@shariant.org.au'
ACCOUNTS_EMAIL = 'accounts@shariant.org.au'

ROLLBAR['enabled'] = False

"""
Settings for AWS SES (Simple Email Service)
"""
EMAIL_BACKEND = 'django_amazon_ses.EmailBackend'

AUTHENTICATION_BACKENDS = (
    'auth.backend.VariantGridOIDCAuthenticationBackend',
    'django.contrib.auth.backends.ModelBackend',  # default
    'guardian.backends.ObjectPermissionBackend',
)

MIDDLEWARE += (
    'mozilla_django_oidc.middleware.SessionRefresh',
    'debug_toolbar.middleware.DebugToolbarMiddleware'
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

VARIANT_CLASSIFICATION_AUTOFUZZ_AGE = True
VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME = False

OIDC_RP_SIGN_ALGO = 'RS256'
OIDC_RP_CLIENT_ID = 'variantgrid'
OIDC_RP_CLIENT_SECRET = get_secret('OIDC.client_secret')
KEY_CLOAK_BASE = 'http://ubuntu-18.04:8080/auth'
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

ANNOTATION_BASE_DIR = "/data"
ANNOTATION_VCF_DUMP_DIR = os.path.join(ANNOTATION_BASE_DIR, 'annotation_scratch')

VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE = False

_SHARIANT_MODE = True
if _SHARIANT_MODE:
    LOGIN_REDIRECT_URL = '/variantclassification/dashboard'

    SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
    SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_test_static")
    STATICFILES_DIRS = (SHARIANT_TEST_STATIC_FILES_DIR, SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

    SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
    TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
    SITE_NAME = "Shariant"

INSTALLED_APPS.append('debug_toolbar')

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
    "variant_classification_dashboard" : True,
    "clinical_contexts": True,

    "variant_classification_import_tool": True,
    "keycloak_admin" : True
})
"""
URLS_NAME_REGISTER.update({  # Disable selected snpdb urls
     "variant_classification_dashboard": True,
     "clinical_contexts": True,
     "hgvs_issues": True,
     "variant_classification_import_tool": True,
     "keycloak_admin": True
})

# mimic shariant
VARIANT_DETAILS_SHOW_ANNOTATION = False
VARIANT_DETAILS_SHOW_SAMPLES = False


UNSHARED_FLAG_ENABLED = True
os.environ['OAUTHLIB_INSECURE_TRANSPORT'] = 'true'
