""" variantgrid.com """

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.annotation_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import


# SYNC_DETAILS = get_shariant_sync_secrets()

# import all the base settings #
SITE_ID = 3  # vg.com

SITE_DESCRIPTION = (
    "VariantGrid is an open source web platform for storing, annotating, analysing and "
    "classifying human variants."
)

WEB_HOSTNAME = 'variantgrid.com'
WEB_IP = '3.104.38.188'

DEBUG = False

PARTITION_ARCHIVE_DIR = "/data/database_dumps/partition_dumps/"
ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'

# This deployment used dbNSFP rankscores before raw scores - keep them visible (see annotation_settings.py)
ANNOTATION_SHOW_LEGACY_RANKSCORES = True

_BIG_DISK_BASE_DIR = "/data"
ANNOTATION_REFERENCE_BASE_DIR = os.path.join(_BIG_DISK_BASE_DIR, "annotation")
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.expanduser("~/perlbrew_runner.sh")

ANNOTATION[BUILD_GRCH38]["enabled"] = True
ANNOTATION[BUILD_T2TV2]["enabled"] = True
ANNOTATION_VCF_DUMP_DIR = os.path.join(_BIG_DISK_BASE_DIR, 'annotation_scratch')

GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1  # MedEx
DEFAULT_FROM_EMAIL = 'noreply@variantgrid.com'
SEND_EMAILS = True
ADMIN_EMAIL_NOTIFICATION = "admin@variantgrid.com"
CONTACT_US_ENABLED = True

# variantgrid.com static overrides (e.g. vc_settings.js enables the 'public' / "3rd Party
# Databases" share level - the gate for sharing outward to ClinVar / MatchMaker Exchange).
VGAWS_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "vgaws_static")
STATICFILES_DIRS = (VGAWS_STATIC_FILES_DIR,) + STATICFILES_DIRS

# MatchMaker Exchange - our own contact, sent as the `contact` block on every outbound patient.
MME_CONTACT = {
    "name": "Centre for Cancer Biology",
    "href": "mailto:david.lawrence@adelaide.edu.au",
    "institution": "Adelaide University",
}


# Needed in production (when debug=False)
ALLOWED_HOSTS = ['variantgrid.com', 'www.variantgrid.com', WEB_HOSTNAME, WEB_IP]
CSRF_TRUSTED_ORIGINS = [f"https://{WEB_HOSTNAME}"]

SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

HGVS_DEFAULT_METHOD = "biocommons_hgvs"

PEDIGREE_MADELINE2_COMMAND = "madeline2"

# Lock down menu - hide some VariantGrid urls / menu
URLS_NAME_REGISTER.update({"sequencing_data": False})

RECAPTCHA_PUBLIC_KEY = get_secret('RECAPTCHA.public_key')
RECAPTCHA_PRIVATE_KEY = get_secret('RECAPTCHA.private_key')
REGISTRATION_OPEN = True

SOMALIER["enabled"] = True
# Somalier reference/1k data + binary moved off slow EFS to local disk - the ancestry step
# reads ~2500 small 1kg-somalier/*.somalier files per run, which stalls badly on EFS.
SOMALIER["annotation_base_dir"] = "/opt/annotation/somalier"

UPLOAD_ENABLED = True

USER_CREATE_ORG_LABS = {
    "unknown": "unknown",
}

USER_CREATE_ORG_MESSAGE = {
    "unknown": "Users must belong to a lab and organisation to perform classifications, and share data with others. "
               "Because we don't know who you are yet, we've assinged you to 'Unknown'. If you want "
               "to properly set things up, please contact david.lawrence@"
               "sa.gov.au (using your institutional email). Thanks!",
}
