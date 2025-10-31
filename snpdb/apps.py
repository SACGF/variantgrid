import logging
import mimetypes
import sys
from typing import AnyStr, Any

from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_save


class SnpdbConfig(AppConfig):
    name = 'snpdb'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from snpdb.models import Trio
        from django.contrib.auth.models import User, Group
        from seqauto.signals.signals_list import backend_vcf_import_success_signal
        from snpdb.signals.signal_handlers import backend_vcf_import_success_handler, trio_post_save_handler, \
            user_post_save_handler, group_post_save_handler
        from snpdb.signals import vcf_health_check
        from snpdb.signals import disk_usage_health_check
        from snpdb.signals import lab_search
        from snpdb.signals import organization_search
        from snpdb.signals import user_search
        from snpdb.signals import cohort_search
        from snpdb.signals import sample_search
        from snpdb.signals import vcf_search
        from snpdb.signals import variant_search
        from snpdb.signals import variant_zygosity_preview_extra
        from snpdb.signals import clinvar_export_search
        from snpdb.signals import scv_search
        from snpdb.signals import genomics_search
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

        if settings.DEBUG:
            # override print statement to include line number
            import traceback

            def dprint(*args, **kwargs):
                '''Pre-pends the filename and linenumber to the print
                statement'''
                stack = traceback.extract_stack()[:-1]
                last = stack[-1]
                filename: str
                line_no: Any
                # Handle different versions of the traceback module
                if hasattr(last, 'filename'):
                    filename = last.filename
                    line_no = last.lineno
                else:
                    filename = last[0]
                    line_no = last[1]

                out_str = ""
                if "VariantGrid" in filename:
                    out_str = "{}:{} - ".format(filename, line_no)

                # Prepend the filename and linenumber
                __builtins__['oldprint'](out_str, *args, **kwargs)


            if 'oldprint' not in __builtins__:
                __builtins__['oldprint'] = __builtins__['print']
            __builtins__['print'] = dprint