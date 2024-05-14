"""
Shariant - https://shariant.org.au

See https://github.com/sacgf/variantgrid/wiki/Annotation%20Setup

"""

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.env.shariantcommon import *  # pylint: disable=wildcard-import, unused-wildcard-import

# import all the base settings #
SITE_ID = 5  # shariant.org.au
SITE_NAME = "Shariant"

# HEARTBEAT_URL = 'https://heartbeat.uptimerobot.com/m788641874-4c58c98a716180f36670e551a0bd03fff47abfea'
SEND_EMAILS = True
OIDC_RP_CLIENT_ID = 'shariant'
OIDC_REQUIRED_GROUP = '/variantgrid/shariant_production'
LOGIN_URL = '/oidc_login/'
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=https%3A%2F%2Fshariant.org.au'


# Prod hasn't been updated yet, default now is 110
ANNOTATION_VEP_VERSION = "108"
ANNOTATION_VEP_VERSION_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_code", ANNOTATION_VEP_VERSION)
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "ensembl-vep")
ANNOTATION_VEP_PLUGINS_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "plugins")
