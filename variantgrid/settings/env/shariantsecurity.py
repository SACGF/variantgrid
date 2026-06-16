"""
Shariant Demo - https://demo.shariant.org.au
"""

from variantgrid.settings.env.shariantcommon import *  # pylint: disable=wildcard-import, unused-wildcard-import

# we import Shariant settings, all that's left is the overrides

# import all the base settings #
SITE_ID = 9
SITE_NAME = "Security Shariant"
SITE_MESSAGE = "This is the security testing version of Shariant. Please avoid sharing real data in this environment."
SLACK['emoji'] = ':cop:'

# OIDC SETTINGS
OIDC_RP_CLIENT_ID = 'shariant-security'
OIDC_REQUIRED_GROUP = '/variantgrid/shariant_security'
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=https%3A%2F%2Ftest2.shariant.org.au'
_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'security_annotation_scratch')

# SEARCH_HGVS_GENE_SYMBOL_USE_MANE = True

SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_security_static")
STATICFILES_DIRS = (SHARIANT_TEST_STATIC_FILES_DIR, SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS
SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)

# Testing out Keycloak 26

# OIDC SETTINGS
USE_OIDC = True
LOGIN_URL = '/oidc_login/'
OIDC_RP_SIGN_ALGO = 'RS256'
# TODO, change this then move it to
OIDC_RP_CLIENT_SECRET = get_secret('OIDC.client_secret')

KEY_CLOAK_BASE = 'https://auth.shariant.org.au'
KEY_CLOAK_REALM = 'agha'
KEY_CLOAK_PROTOCOL_BASE = KEY_CLOAK_BASE + '/realms/' + KEY_CLOAK_REALM + '/protocol/openid-connect'
OIDC_OP_JWKS_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/certs'
OIDC_OP_AUTHORIZATION_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/auth'
OIDC_OP_TOKEN_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/token'
OIDC_OP_USER_ENDPOINT = KEY_CLOAK_PROTOCOL_BASE + '/userinfo'
OIDC_USER_SERVICES = KEY_CLOAK_BASE + '/realms/' + KEY_CLOAK_REALM + '/account'
OIDC_OP_LOGOUT_URL_METHOD = 'auth.backend.provider_logout'