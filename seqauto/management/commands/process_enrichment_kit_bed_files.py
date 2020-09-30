from django.core.management.base import BaseCommand
import logging

from library.guardian_utils import admin_bot
from seqauto.models import EnrichmentKit
from snpdb.models.models_enums import ImportSource, ProcessingStatus
from upload.models import UploadedFile, UploadedFileTypes, UploadedBed
from upload.upload_processing import process_uploaded_file


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        user = admin_bot()

        if options.get("clear"):
            logging.info("Clearing existing enrichment_kit intervals")
            for enrichment_kit in EnrichmentKit.objects.filter(genomic_intervals__isnull=False):
                enrichment_kit.genomic_intervals.delete()

        unprocessed_enrichment_kits = EnrichmentKit.objects.filter(bed_file__isnull=False, genomic_intervals__isnull=True)
        for enrichment_kit in unprocessed_enrichment_kits:
            name = f"EnrichmentKit {enrichment_kit} bed file"
            logging.info("Processing %s", name)
            uploaded_file = UploadedFile.objects.create(path=enrichment_kit.bed_file,
                                                        import_source=ImportSource.COMMAND_LINE,
                                                        name=name,
                                                        user=user,
                                                        file_type=UploadedFileTypes.BED)

            process_uploaded_file(uploaded_file, run_async=False)
            upload_pipeline = uploaded_file.uploadpipeline

            if upload_pipeline.status == ProcessingStatus.SUCCESS:
                uploaded_bed = UploadedBed.objects.get(uploaded_file=uploaded_file)
                enrichment_kit.genomic_intervals = uploaded_bed.genomic_intervals_collection
                enrichment_kit.save()
            else:
                logging.error("%s processing_status: %s", enrichment_kit.bed_file, upload_pipeline.get_status_display())
