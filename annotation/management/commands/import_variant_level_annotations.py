#!/usr/bin/env python3
"""
    This is used when an Annotation import crashes, and we want to import the VCF without having to re-run the
    annotation which will require re-annotating the VCF
"""

from django.core.management.base import BaseCommand
from django.utils import timezone
import logging
import os

from annotation.models.models import AnnotationRun, VariantAnnotation, VariantTranscriptAnnotation
from annotation.signals import annotation_run_complete_signal
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--annotation_run_id', type=int, required=True)
        parser.add_argument('--test-no-insert', action='store_true')

    def handle(self, *args, **options):
        annotation_run_id = options["annotation_run_id"]
        insert_variants = not options["test_no_insert"]

        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        annotation_run.error_exception = None
        annotation_run.save()

        logging.info("Deleting Variant/Transcript annotation within range lock")
        variant_annotation_version = annotation_run.annotation_range_lock.version
        range_kwargs = {"version": variant_annotation_version,
                        "variant_id__gte": annotation_run.annotation_range_lock.min_variant,
                        "variant_id__lte": annotation_run.annotation_range_lock.max_variant}
        VariantAnnotation.objects.filter(**range_kwargs)
        VariantTranscriptAnnotation.objects.filter(**range_kwargs)

        logging.info("Starting import VCF annotations")
        import_vcf_annotations(annotation_run, insert_variants=insert_variants)
        logging.info("Finished import, setting complete")

        annotation_run.upload_end = timezone.now()
        annotation_run.save()

        annotation_run_complete_signal.send(sender=os.path.basename(__file__),
                                            variant_annotation_version=variant_annotation_version)
