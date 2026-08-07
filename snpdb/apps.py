import logging
import mimetypes
import sys

from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_save


class SnpdbConfig(AppConfig):
    name = 'snpdb'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from django.contrib.auth.models import Group, User

        from seqauto.signals.signals_list import backend_vcf_import_success_signal
        from snpdb.models import Trio
        from snpdb.signals.signal_handlers import (
            backend_vcf_import_success_handler,
            group_post_save_handler,
            trio_post_save_handler,
            user_post_save_handler,
        )
        # pylint: enable=import-outside-toplevel,unused-import

        backend_vcf_import_success_signal.connect(backend_vcf_import_success_handler)

        if not settings.UNIT_TEST:
            # Add newly created users to public group
            post_save.connect(user_post_save_handler, sender=User)

        if settings.REQUESTS_DISABLE_IPV6:
            import requests  # pylint: disable=import-outside-toplevel
            requests.packages.urllib3.util.connection.HAS_IPV6 = False

        # Make global settings share read only with this group by default
        post_save.connect(group_post_save_handler, sender=Group)

        # Add newly created users to public group
        post_save.connect(trio_post_save_handler, sender=Trio)

        # Disable annoying matplotlib findfont messages
        logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

        if sys.version_info < (3, 10):
            raise SystemExit("VariantGrid requires Python 3.10 or later.")

        # So static serve of MEDIA_ROOT just prompts to download VCF (used to rename to vcf.vcf)
        mimetypes.add_type("application/octet-stream", ".vcf", strict=True)
