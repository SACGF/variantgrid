from django.core.management.base import BaseCommand

from annotation.models import AnnotationRun
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from annotation.vcf_files.import_vcf_annotations import handle_vep_warnings


class Command(BaseCommand):
    def handle(self, *args, **options):
        for ar in AnnotationRun.objects.filter(vep_skipped_count__isnull=True):
            bulk_inserter = BulkVEPVCFAnnotationInserter(ar, validate_columns=False)
            bulk_inserter.batch_id = 99999  # So we don't have to worry about overwriting old temp data
            handle_vep_warnings(ar, bulk_inserter)
            ar.save()
