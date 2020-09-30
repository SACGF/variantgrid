from django.core.management.base import BaseCommand, CommandError
import logging
import os

from upload.models import UploadedFileTypes, UploadPipeline
from upload.uploaded_file_type import retry_upload_pipeline


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--uploaded_file_type', required=True)

    def handle(self, *args, **options):
        uploaded_file_type = options['uploaded_file_type']

        uft_dict = dict(UploadedFileTypes.CHOICES)
        uft_description = uft_dict.get(uploaded_file_type)
        if uft_description is None:
            script_name = os.path.basename(__file__)
            valid_ufts = ','.join(sorted(uft_dict))
            msg = f"Usage: {script_name} --uploaded_file_type=X (where X is one of {valid_ufts})"
            raise CommandError(msg)

        logging.info("Reloading UFPPs of type '%s'", uft_description)
        for upload_pipeline in UploadPipeline.objects.filter(uploaded_file__file_type=uploaded_file_type):
            retry_upload_pipeline(upload_pipeline)
