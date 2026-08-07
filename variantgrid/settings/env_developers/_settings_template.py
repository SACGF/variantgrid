"""
Settings template for a development machine - copy to env_developers/<hostname>.py

Deployed servers copy env/_settings_template.py instead, which is served over HTTPS and has examples
for customising VEP version, annotation consortium and genome builds.
"""
from variantgrid.settings.components.annotation_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

DEBUG = True

WEB_HOSTNAME = 'localhost'
WEB_IP = '127.0.0.1'

ALLOWED_HOSTS = ["localhost", WEB_HOSTNAME, WEB_IP]

# ANNOTATION_ENTREZ_EMAIL = 'your_email@yourdomain.com'

# PEDIGREE_MADELINE2_COMMAND = "madeline2"

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
