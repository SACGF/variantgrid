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
LOGOUT_REDIRECT_URL = "https://test2.shariant.org.au"

# SCRATCH
_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'security_annotation_scratch')

# SEARCH_HGVS_GENE_SYMBOL_USE_MANE = True

SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_security_static")
STATICFILES_DIRS = (SHARIANT_TEST_STATIC_FILES_DIR, SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS
SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)

