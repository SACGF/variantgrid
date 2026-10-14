"""
Shariant Demo - https://demo.shariant.org.au
"""

from variantgrid.settings.env.shariantcommon import *  # pylint: disable=wildcard-import, unused-wildcard-import

# we import Shariant settings, all that's left is the overrides

# import all the base settings #
SITE_ID = 8
SITE_NAME = "Shariant Demo"
SITE_SHORT_NAME = "Shariant Demo"
SITE_MESSAGE = "This is the demo version of Shariant. Please avoid sharing real data in this environment."

# Discordance emails are logged, not sent, on this host (SEND_EMAILS is False in shariantcommon)
SILENCED_SYSTEM_CHECKS = [*SILENCED_SYSTEM_CHECKS, "variantgrid.W001"]

SLACK['emoji'] = ':tv:'

# OIDC SETTINGS
OIDC_RP_CLIENT_ID = 'shariant-demo'
OIDC_REQUIRED_GROUP = '/variantgrid/shariant_demo'
LOGOUT_REDIRECT_URL = "https://test2.shariant.org.au"

_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'demo_annotation_scratch')


CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL = "lab"

# SEARCH_HGVS_GENE_SYMBOL_USE_MANE = True

_SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_demo_static")
STATICFILES_DIRS = (_SHARIANT_TEST_STATIC_FILES_DIR,) + STATICFILES_DIRS
SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
