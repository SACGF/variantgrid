#!/usr/bin/env python3

from django.core.management.base import BaseCommand
from django.db.utils import DatabaseError
import logging

from annotation.annotation_versions import get_variant_annotation_version
from annotation.tasks.annotation_scheduler_task import annotation_scheduler
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):

    def handle(self, *args, **options):
        logging.info("Checking/Creating VariantAnnotationVersion...")
        for genome_build in GenomeBuild.builds_with_annotation():
            vav = get_variant_annotation_version(genome_build)
            try:
                # Some upgrade migrations caused partitions to be deleted
                vav.create_partition()
            except DatabaseError:
                pass

        logging.info("Scheduling annotation...")
        annotation_scheduler()
