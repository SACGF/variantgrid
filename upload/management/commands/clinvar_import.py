#!/usr/bin/env python3

from django.core.management.base import BaseCommand

from library.guardian_utils import admin_bot
from library.log_utils import console_logger
from snpdb.models.models_enums import ImportSource
from upload.models import UploadedFileTypes
from upload.upload_processing import process_vcf_file


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('vcf')

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]

        logger = console_logger()
        name = 'ClinVar import'
        user = admin_bot()

        _, result = process_vcf_file(vcf_filename, name, user,
                                     import_source=ImportSource.COMMAND_LINE,
                                     file_type=UploadedFileTypes.CLINVAR)
        if result:
            logger.info(f"Result {result} status = {result.status}")
