import logging

from django.core.management.base import BaseCommand

from annotation.annotation_files_check import annotation_data_exists
from variantgrid.deployment_validation.library_version_checks import check_library_versions
from variantgrid.deployment_validation.tool_version_checks import check_tool_versions


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--quiet', action='store_true', help="Suppress info message")
        parser.add_argument('--die-if-invalid', action='store_true')

    def handle(self, *args, **options):
        quiet = options["quiet"]
        die_if_invalid = options["die_if_invalid"]

        checks = {
            "Annotation data exists": annotation_data_exists(flat=True),
            "Library versions": check_library_versions(),
            "Tool versions": check_tool_versions(),
        }

        for check_type, check in checks.items():
            for k, data in check.items():
                if data.get("valid"):
                    if not quiet:
                        logging.info(f"%s, %s: OK", check_type, k)
                else:
                    fix = data.get("fix")
                    msg = f"{check_type=} {k} INVALID. Fix: {fix}"
                    if die_if_invalid:
                        raise ValueError(msg)
                    else:
                        logging.info(msg)
