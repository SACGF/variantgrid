"""
Settings for running the unit test suite on GitHub Actions (.github/workflows/django-tests.yml).

The workflow sets DJANGO_SETTINGS_MODULE=variantgrid.settings.env.github_actions explicitly,
as hostname-based settings detection can't work on ephemeral CI runners.

Database credentials deliberately match the defaults in secret_settings.py (snpdb/snpdb@localhost),
which is what the workflow's postgres service container provides - so no settings_config.json
is needed on the runner.
"""
from variantgrid.settings.components.annotation_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.celery_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.default_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import
from variantgrid.settings.components.seqauto_settings import *  # pylint: disable=wildcard-import, unused-wildcard-import

WEB_HOSTNAME = 'localhost'
WEB_IP = '127.0.0.1'

ALLOWED_HOSTS = ["localhost", WEB_HOSTNAME, WEB_IP]

ANNOTATION[BUILD_GRCH37].update({
    "annotation_consortium": "RefSeq",
})
ANNOTATION[BUILD_GRCH38].update({
    "annotation_consortium": "RefSeq",
    "enabled": True,
})

ANNOTATION[BUILD_T2TV2]["enabled"] = True
