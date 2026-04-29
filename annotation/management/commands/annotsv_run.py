#!/usr/bin/env python3
"""
Re-run the AnnotSV stage for a single SV AnnotationRun.

Used to recover from a failed best-effort AnnotSV stage without re-running VEP,
or to backfill AnnotSV annotations for existing SV runs after enabling AnnotSV
on a deployment.
"""

import logging
import os

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from annotation.annotsv_annotation import run_annotsv
from annotation.models.models import AnnotationRun
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.vcf_files.bulk_annotsv_tsv_inserter import import_annotsv_tsv


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--annotation-run", type=int, required=True,
                            help="AnnotationRun PK")
        parser.add_argument("--skip-run", action="store_true",
                            help="Skip running AnnotSV; only import the existing TSV")

    def handle(self, *args, **options):
        if not settings.ANNOTATION_ANNOTSV_ENABLED:
            raise CommandError("ANNOTATION_ANNOTSV_ENABLED is False - enable it in settings first")

        run_id = options["annotation_run"]
        annotation_run = AnnotationRun.objects.get(pk=run_id)

        if annotation_run.pipeline_type != VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
            raise CommandError(f"AnnotationRun {run_id} is not a STRUCTURAL_VARIANT run")

        if not options["skip_run"]:
            if not annotation_run.vcf_dump_filename:
                raise CommandError(f"AnnotationRun {run_id} has no vcf_dump_filename - run VEP first")

            annotsv_dir = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR,
                                       f"annotsv_{annotation_run.pk}")
            tsv, rc, _stdout, stderr = run_annotsv(
                annotation_run.vcf_dump_filename, annotsv_dir,
                annotation_run.genome_build, annotation_run.annotation_consortium,
            )
            if rc == 0 and os.path.exists(tsv):
                annotation_run.annotsv_tsv_filename = tsv
                annotation_run.annotsv_error = None
            else:
                annotation_run.annotsv_error = (stderr or "")[:100_000]
                annotation_run.save()
                raise CommandError(f"AnnotSV run failed (rc={rc}): {stderr[:1000]}")
            annotation_run.save()

        if not annotation_run.annotsv_tsv_filename:
            raise CommandError(f"AnnotationRun {run_id} has no annotsv_tsv_filename")

        updated = import_annotsv_tsv(annotation_run)
        logging.info("Updated %d VariantAnnotation rows", updated)
        self.stdout.write(self.style.SUCCESS(f"AnnotSV import complete: {updated} rows updated"))
