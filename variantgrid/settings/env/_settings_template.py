from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

# ANNOTATION_ENTREZ_EMAIL = 'your_email@yourdomain.com'

WEB_HOSTNAME = 'yourdomain.com'
WEB_IP = '127.0.0.1'

ALLOWED_HOSTS = ["localhost", WEB_HOSTNAME, WEB_IP]

# SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

# PEDIGREE_MADELINE2_COMMAND = "madeline2"


# Here is how you'd customise your VEP version and whether you want RefSeq etc:

_different_vep_version = False
if _different_vep_version:
    ANNOTATION_VEP_VERSION = "112"
    ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
    ANNOTATION_VEP_VERSION_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_code", ANNOTATION_VEP_VERSION)
    ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "ensembl-vep")
    ANNOTATION_VEP_PLUGINS_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "plugins")

_use_grch38 = False
if _use_grch38:
    ANNOTATION[BUILD_GRCH38]["enabled"] = True

_use_refseq = True
if _use_refseq:
    ANNOTATION[BUILD_GRCH37].update({
        "annotation_consortium": "RefSeq",
    })
    ANNOTATION[BUILD_GRCH38].update({
        "annotation_consortium": "RefSeq",
    })

