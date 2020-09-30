from django.core.management.base import BaseCommand
import sys

from library.file_utils import open_handle_gzip
from upload.models import UploadStep, UploadPipeline
from upload.vcf.import_vcf_extract_unknown_variants_step import import_vcf_extract_unknown_variants_and_split_file


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--upload-step-pk', help='Upload Step Primary Key', required=True)

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]
        upload_step_pk = options["upload_step_pk"]

        if vcf_filename == '-':
            vcf_file = sys.stdin
        else:
            vcf_file = open_handle_gzip(vcf_filename, mode='rt')

        upload_step = UploadStep.objects.get(pk=upload_step_pk)
        items_processed = import_vcf_extract_unknown_variants_and_split_file(upload_step, vcf_file)

        # This step is the "variant count" (post normalisation) so use it for both upload_step and total pipeline items
        UploadStep.objects.filter(pk=upload_step.pk).update(items_processed=items_processed,
                                                            items_to_process=items_processed)

        upload_pipeline_qs = UploadPipeline.objects.filter(pk=upload_step.upload_pipeline.pk)
        upload_pipeline_qs.update(items_processed=items_processed)
