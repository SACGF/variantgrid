""" variantgrid.com """

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

# SYNC_DETAILS = get_shariant_sync_secrets()

# import all the base settings #
SITE_ID = 3  # vg.com

WEB_HOSTNAME = 'test.variantgrid.com'
WEB_IP = '203.101.229.193'

DEBUG = False

ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'
SLACK['emoji'] = ':mouse:'

_BIG_DISK_BASE_DIR = "/data"
ANNOTATION_REFERENCE_BASE_DIR = os.path.join(_BIG_DISK_BASE_DIR, "annotation")
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "VEP")
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "ensembl-vep")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")
_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "fasta")


ANNOTATION[BUILD_GRCH37].update({
    "columns_version": 3,
    "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.25_GRCh37.p13_genomic.fna.gz"),
})

ANNOTATION[BUILD_GRCH38].update({
    "columns_version": 3,
    "enabled": True,
    "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
})


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

ANNOTATION_VCF_DUMP_DIR = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, 'scratch')

GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1  # MedEx
DEFAULT_FROM_EMAIL = 'noreply@variantgrid.com'
SEND_EMAILS = True
ADMIN_EMAIL_NOTIFICATION = "admin@variantgrid.com"
CONTACT_US_ENABLED = True

# Needed in production (when debug=False)
ALLOWED_HOSTS = [WEB_HOSTNAME, WEB_IP]

HEALTH_CHECK_ENABLED = False
HGVS_DEFAULT_METHOD = "biocommons_hgvs"

SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

PEDIGREE_MADELINE2_COMMAND = "madeline2"

RECAPTCHA_PUBLIC_KEY = get_secret('RECAPTCHA.public_key')
RECAPTCHA_PRIVATE_KEY = get_secret('RECAPTCHA.private_key')

#UPLOAD_ENABLED = False

# Lock down menu - hide some VariantGrid urls / menu
URLS_NAME_REGISTER.update({"sequencing_data": False})

VG_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "variantgrid_test_static")
STATICFILES_DIRS = (VG_TEST_STATIC_FILES_DIR,) + STATICFILES_DIRS

SOMALIER["enabled"] = True
SOMALIER["annotation_base_dir"] = os.path.join(ANNOTATION_REFERENCE_BASE_DIR, "somalier")

USER_CREATE_ORG_LABS = {
    "unknown": "unknown",
}

USER_CREATE_ORG_MESSAGE = {
    "unknown": "Users must belong to a lab and organisation to perform classifications, and share data with others. "
               "Because we don't know who you are yet, we've assinged you to 'Unknown'. If you want "
               "to properly set things up, please contact david.lawrence@"
               "sa.gov.au (using your institutional email). Thanks!",
}
