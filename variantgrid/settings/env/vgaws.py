""" variantgrid.com """

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import


# SYNC_DETAILS = get_shariant_sync_secrets()

# import all the base settings #
SITE_ID = 3  # vg.com

WEB_HOSTNAME = 'variantgrid.com'
WEB_IP = '54.253.175.206'

DEBUG = False
PYTHON_COMMAND = os.path.join(BASE_DIR, "env", "bin", "python3")
MANAGE_COMMAND = [PYTHON_COMMAND, os.path.join(BASE_DIR, "manage.py")]

ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'

_BIG_DISK_BASE_DIR = "/data"
ANNOTATION_REFERENCE_BASE_DIR = os.path.join(_BIG_DISK_BASE_DIR, "annotation")
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "VEP")
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "ensembl-vep")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.expanduser("~/perlbrew_runner.sh")
_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "fasta")

ANNOTATION[BUILD_GRCH38]["enabled"] = True
ANNOTATION_VCF_DUMP_DIR = os.path.join(_BIG_DISK_BASE_DIR, 'annotation_scratch')

GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1  # MedEx
DEFAULT_FROM_EMAIL = 'noreply@variantgrid.com'
SEND_EMAILS = True
ADMIN_EMAIL_NOTIFICATION = "admin@variantgrid.com"
CONTACT_US_ENABLED = True


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
SOMALIER["annotation_base_dir"] = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "somalier")

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
