from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

# ANNOTATION_ENTREZ_EMAIL = 'your_email@yourdomain.com'

WEB_HOSTNAME = 'yourdomain.com'
WEB_IP = '127.0.0.1'

ALLOWED_HOSTS = [WEB_HOSTNAME, WEB_IP]

# SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTOCOL', 'https')

# PEDIGREE_MADELINE2_COMMAND = "madeline2"
