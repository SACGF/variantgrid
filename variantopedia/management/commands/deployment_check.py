import logging

from django.conf import settings
from django.core.management.base import BaseCommand

from variantgrid.deployment_validation.annotation_files_check import annotation_data_exists, check_cdot_data
from variantgrid.deployment_validation.annotation_status_checks import check_annotation_versions, check_variant_annotation_runs_status
from variantgrid.deployment_validation.celery_checks import check_celery_tasks
from variantgrid.deployment_validation.column_check import check_variantgrid_columns
from variantgrid.deployment_validation.library_version_checks import check_library_versions
from variantgrid.deployment_validation.somalier_check import check_somalier
from variantgrid.deployment_validation.tool_version_checks import check_tool_versions
from variantgrid.deployment_validation.vep_check import check_vep


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--quiet', action='store_true', help="Suppress info message")
        parser.add_argument('--die-if-invalid', action='store_true')

    def handle(self, *args, **options):
        quiet = options["quiet"]
        die_if_invalid = options["die_if_invalid"]

        checks = {
            "Annotation data exists": annotation_data_exists(flat=True),
            "Annotation Versions": check_annotation_versions(),
            "Variant Annotation status": check_variant_annotation_runs_status(),
            "Library versions": check_library_versions(),
            "Tool versions": check_tool_versions(),
            "cdot data": check_cdot_data(),
            "Celery Tasks": check_celery_tasks(),
            "Columns": check_variantgrid_columns(),
            "VEP": check_vep(),
        }
        if settings.SOMALIER.get("enabled"):
            checks["somalier"] = check_somalier()

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
                        logging.error(msg)

                if warning := data.get("warning"):
                    logging.warning(warning)
