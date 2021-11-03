"""
Shariant - https://shariant.org.au

See https://github.com/sacgf/variantgrid/wiki/Annotation%20Setup

"""

from variantgrid.settings.env.shariant import *  # pylint: disable=wildcard-import, unused-wildcard-import

# we import Shariant settings, all that's left is the overrides

# import all the base settings #
SITE_ID = 6

URLS_NAME_REGISTER.update({
    "condition_aliases": False,
    "clinvar_key_summary": True
})

DISCORDANCE_EMAIL = '"Shariant Updates" <updates@shariant.org.au>'
SEND_EMAILS = False
HEALTH_CHECK_ENABLED = False

HEARTBEAT_URL = 'https://heartbeat.uptimerobot.com/m785768698-f47e172d7c428063e2f452af783dadfcc07da004'

# OIDC SETTINGS
OIDC_RP_CLIENT_ID = 'shariant-test'
OIDC_REQUIRED_GROUP = None
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=https%3A%2F%2Ftest.shariant.org.au'

_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'test_annotation_scratch')

LIFTOVER_NCBI_REMAP_ENABLED = True

SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_test_static")
STATICFILES_DIRS = (SHARIANT_TEST_STATIC_FILES_DIR, SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
SITE_NAME = "Shariant Test"

VARIANT_VCF_DB_PREFIX = "stv"

SITE_MESSAGE = "This is the test version of Shariant. Please avoid sharing real data in this environment."
# DISCORDANCE_PAUSE_TEMP_VARIANT_MATCHING = True