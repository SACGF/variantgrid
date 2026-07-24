from variantgrid.settings.components.annotation_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

# ANNOTATION_ENTREZ_EMAIL = 'your_email@yourdomain.com'

WEB_HOSTNAME = 'test.variantgrid.com'
WEB_IP = '129.127.17.180'

ALLOWED_HOSTS = ["localhost", WEB_HOSTNAME, WEB_IP]

# SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

# PEDIGREE_MADELINE2_COMMAND = "madeline2"
# This box runs the latest annotation, which is now the package default (#1625) so needs no VEP overrides here.

ANNOTATION_ENTREZ_EMAIL = 'davmlaw@gmail.com'
SLACK['emoji'] = ':mouse:'

ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

ANNOTATION_VCF_DUMP_DIR = os.path.join(ANNOTATION_BASE_DIR, 'scratch')

# GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID = 1  # MedEx
HEALTH_CHECK_ENABLED = False
HGVS_DEFAULT_METHOD = "biocommons_hgvs"

# Beacon v2 genomic data-sharing endpoint (#1661) - enabled on test.variantgrid.com for
# manual testing. environment="test" is the spec's own signal that this is a non-production
# instance. (HEALTH_CHECK_ENABLED=False above, so the nightly Slack line won't post here -
# inspect BeaconInboundQuery rows via admin when testing.)
BEACON_ENABLED = True
BEACON_CONFIG = {
    **BEACON_CONFIG,
    "beacon_id": "com.variantgrid.test.beacon",
    "environment": "test",
    "organization": {
        "id": "variantgrid",
        "name": "VariantGrid (test)",
        "welcome_url": "https://test.variantgrid.com/",
        "contact_url": "mailto:admin@variantgrid.com",
    },
}

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

# Latest annotation (VEP 115, columns_version 4, dbNSFP 5.x, gnomAD 4.1, ...) is the package default now (#1625)
# ANNOTATION_ANNOTSV_ENABLED = True