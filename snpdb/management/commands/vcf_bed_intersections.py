from django.core.management.base import BaseCommand
import logging

from snpdb.models import VCFBedIntersection
from snpdb.models.models_enums import ProcessingStatus
from snpdb.tasks.vcf_bed_file_task import create_backend_vcf_bed_intersections
from upload.models import BackendVCF


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--clear-errors', action='store_true')
        parser.add_argument('--clear-non-success', action='store_true')
        parser.add_argument('--clear-all', action='store_true')

    def handle(self, *args, **options):
        clear_all = options.get("clear_all")
        clear_errors = options.get("clear_errors")
        clear_non_success = options.get("clear_non_success")

        if clear_all:
            VCFBedIntersection.objects.all().delete()
        elif clear_errors:
            VCFBedIntersection.objects.filter(status=ProcessingStatus.ERROR).delete()
        elif clear_non_success:
            VCFBedIntersection.objects.exclude(status=ProcessingStatus.SUCCESS).delete()

        if not clear_non_success:
            non_success_pbis = VCFBedIntersection.objects.exclude(status=ProcessingStatus.SUCCESS)
            if non_success_pbis.exists():
                logging.warning("There are %d VCFBedIntersections that have status != SUCCESS", non_success_pbis.count())
                for pbi in non_success_pbis:
                    logging.warning(str(pbi))

        # This will create a lot of warnings for duplicates if already run
        for backend_vcf in BackendVCF.objects.all():  # filter(combo_vcf__sample_sheet__sequencing_sample__enrichment_kit__in=enrichment_kit_ids):
            create_backend_vcf_bed_intersections(backend_vcf)
