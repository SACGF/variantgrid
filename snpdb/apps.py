import logging

from django.apps import AppConfig
from django.db.models.signals import post_save


class SnpdbConfig(AppConfig):
    name = 'snpdb'

    def ready(self):
        # pylint: disable=import-outside-toplevel
        from snpdb.models import Trio
        from django.contrib.auth.models import User, Group
        from seqauto.signals import backend_vcf_import_success_signal
        from snpdb.signals.signal_handlers import backend_vcf_import_success_handler, trio_post_save_handler, \
            user_post_save_handler, group_post_save_handler
        from snpdb.signals import vcf_health_check  # pylint: disable=unused-import
        from snpdb.signals import disk_usage_health_check  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel

        backend_vcf_import_success_signal.connect(backend_vcf_import_success_handler)

        # Add newly created users to public group
        post_save.connect(user_post_save_handler, sender=User)

        # Make global settings share read only with this group by default
        post_save.connect(group_post_save_handler, sender=Group)

        # Add newly created users to public group
        post_save.connect(trio_post_save_handler, sender=Trio)

        # Disable annoying matplotlib findfont messages
        logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
