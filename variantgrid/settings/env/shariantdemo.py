"""
Shariant Demo - https://demo.shariant.org.au
"""

from variantgrid.settings.env.shariant import *  # pylint: disable=wildcard-import, unused-wildcard-import

# we import Shariant settings, all that's left is the overrides

# import all the base settings #
SITE_ID = 8

OIDC_REQUIRED_GROUP = '/variantgrid/shariant_demo'

SLACK['emoji'] = ':tv:'

URLS_NAME_REGISTER.update({
    "condition_aliases": False,
    "clinvar_key_summary": True
})

DISCORDANCE_EMAIL = '"Shariant Updates" <updates@shariant.org.au>'
SEND_EMAILS = False
HEALTH_CHECK_ENABLED = True

HEARTBEAT_URL = 'https://heartbeat.uptimerobot.com/m785768698-f47e172d7c428063e2f452af783dadfcc07da004'

# CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL = "lab"
CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL = "logged_in_users"
CLASSIFICATION_OMNI_IMPORTER_INCLUDE_SOURCE = True

# OIDC SETTINGS
OIDC_RP_CLIENT_ID = 'shariant-demo'
OIDC_REQUIRED_GROUP = None
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=https%3A%2F%2Fdemo.shariant.org.au'

_ANNOTATION_BASE_DIR = "/data/annotation"  # Set this to where you downloaded annotation (${ANNOTATION_BASE_DIR} from wiki)
ANNOTATION_VCF_DUMP_DIR = os.path.join(_ANNOTATION_BASE_DIR, 'demo_annotation_scratch')

ANNOTATION[BUILD_GRCH37]["columns_version"] = 3
ANNOTATION[BUILD_GRCH38]["columns_version"] = 3

ANNOTATION[BUILD_GRCH37]["vep_config"].update({
    "cosmic": "annotation_data/GRCh37/Cosmic_GenomeScreensMutant_v99_GRCh37.vcf.gz",
    "dbnsfp": "annotation_data/GRCh37/dbNSFP4.5a.grch37.stripped.gz",
    "mastermind": "annotation_data/GRCh37/mastermind_cited_variants_reference-2023.10.02-grch37.vcf.gz",
})
ANNOTATION[BUILD_GRCH38]["vep_config"].update({
    "cosmic": "annotation_data/GRCh38/Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz",
    "dbnsfp": "annotation_data/GRCh38/dbNSFP4.5a.grch38.stripped.gz",
    "mastermind": "annotation_data/GRCh38/mastermind_cited_variants_reference-2023.10.02-grch38.vcf.gz",
})

SEARCH_HGVS_GENE_SYMBOL_USE_MANE = True

SHARIANT_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_static")
SHARIANT_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "shariant_demo_static")
STATICFILES_DIRS = (SHARIANT_TEST_STATIC_FILES_DIR, SHARIANT_STATIC_FILES_DIR,) + STATICFILES_DIRS

SHARIANT_TEMPLATES_DIR = os.path.join(VARIANTGRID_APP_DIR, "templates/shariant_templates")
TEMPLATES[0]["DIRS"].insert(0, SHARIANT_TEMPLATES_DIR)
SITE_NAME = "Shariant Demo"

VARIANT_VCF_DB_PREFIX = "stv"

SITE_MESSAGE = "This is the demo version of Shariant. Please avoid sharing real data in this environment."

URLS_NAME_REGISTER.update({
    "vus": True
})

# VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS = True
