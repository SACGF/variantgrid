from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

# ANNOTATION_ENTREZ_EMAIL = 'your_email@yourdomain.com'

WEB_HOSTNAME = 'test.variantgrid.com'
WEB_IP = '203.101.231.225'

ALLOWED_HOSTS = ["localhost", WEB_HOSTNAME, WEB_IP]

# SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

# PEDIGREE_MADELINE2_COMMAND = "madeline2"
# Here is how you'd customise your VEP version and whether you want RefSeq etc:

ANNOTATION_VEP_VERSION = "112"
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
ANNOTATION_VEP_VERSION_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_code", ANNOTATION_VEP_VERSION)
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "ensembl-vep")
ANNOTATION_VEP_PLUGINS_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "plugins")

ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'
SLACK['emoji'] = ':mouse:'

ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

ANNOTATION_VCF_DUMP_DIR = os.path.join(ANNOTATION_BASE_DIR, 'scratch')

# GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1  # MedEx
HEALTH_CHECK_ENABLED = False
HGVS_DEFAULT_METHOD = "biocommons_hgvs"



ANNOTATION[BUILD_GRCH37].update({
    "annotation_consortium": "RefSeq",
})
ANNOTATION[BUILD_GRCH38].update({
    "annotation_consortium": "RefSeq",
    "enabled": True,
})

ANNOTATION[BUILD_T2TV2]["enabled"] = True

PEDIGREE_MADELINE2_COMMAND='madeline2'
LIFTOVER_BCFTOOLS_ENABLED = True

SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

RECAPTCHA_PUBLIC_KEY = get_secret('RECAPTCHA.public_key')
RECAPTCHA_PRIVATE_KEY = get_secret('RECAPTCHA.private_key')

# Stop Chrome warnings for non-https test site
SECURE_CROSS_ORIGIN_OPENER_POLICY = None

#UPLOAD_ENABLED = False

# Lock down menu - hide some VariantGrid urls / menu
URLS_NAME_REGISTER.update({"sequencing_data": False})

VG_TEST_STATIC_FILES_DIR = os.path.join(VARIANTGRID_APP_DIR, "static_files", "variantgrid_test_static")
STATICFILES_DIRS = (VG_TEST_STATIC_FILES_DIR,) + STATICFILES_DIRS

SOMALIER["enabled"] = True
SOMALIER["annotation_base_dir"] = os.path.join(ANNOTATION_BASE_DIR, "somalier")

USER_CREATE_ORG_LABS = {
    "unknown": "unknown",
}

USER_CREATE_ORG_MESSAGE = {
    "unknown": "Users must belong to a lab and organisation to perform classifications, and share data with others. "
               "Because we don't know who you are yet, we've assinged you to 'Unknown'. If you want "
               "to properly set things up, please contact david.lawrence@"
               "sa.gov.au (using your institutional email). Thanks!",
}
