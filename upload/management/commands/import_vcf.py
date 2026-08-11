#!/usr/bin/env python3

import logging

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from snpdb.models.models_enums import ImportSource
from upload.upload_metadata import GENOME_BUILD, SOURCE
from upload.upload_processing import process_vcf_file


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--name', required=True)
        parser.add_argument('--user', required=True)
        # Upload metadata - facts the VCF header may not carry. @see upload.upload_metadata
        parser.add_argument('--genome-build', help='Declared genome build, used where the header has no contigs')
        parser.add_argument('--source', help='Declared variant caller, matched against VCFSourceSettings')
        parser.add_argument('vcf')

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]
        username = options["user"]
        name = options["name"]
        user = User.objects.get(username=username)

        metadata = {}
        for key, option in {GENOME_BUILD: "genome_build", SOURCE: "source"}.items():
            if value := options.get(option):
                metadata[key] = value

        (_, result) = process_vcf_file(vcf_filename, name, user, import_source=ImportSource.COMMAND_LINE,
                                       metadata=metadata)
        if result:
            logging.info("Result %s status = %s", result, result.status)
