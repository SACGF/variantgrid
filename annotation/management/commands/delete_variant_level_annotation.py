#!/usr/bin/env python3

import logging

from django.core.management.base import BaseCommand

from annotation.models.models import VariantAnnotationVersion, AnnotationRangeLock
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--variant_annotation_version_id', type=int, help='Delete a particular version')
        group.add_argument('--genome-build', help='Delete all for a genome build')
        group.add_argument('--all', action='store_true', help='Delete *everything*')

    def handle(self, *args, **options):
        variant_annotation_version_id = options["variant_annotation_version_id"]
        build_name = options["genome_build"]

        qs = VariantAnnotationVersion.objects.all()  # --all
        if variant_annotation_version_id:
            qs = qs.filter(pk=variant_annotation_version_id)
        elif build_name:
            genome_build = GenomeBuild.get_name_or_alias(build_name)
            qs = qs.filter(genome_build=genome_build)

        for variant_annotation_version in qs:
            logging.info("Deleting all VariantAnnotation for %s", variant_annotation_version)
            variant_annotation_version.truncate_related_objects()

            # Will cascade delete AnnotationRun then VariantAnnotation
            AnnotationRangeLock.objects.filter(version=variant_annotation_version).delete()
