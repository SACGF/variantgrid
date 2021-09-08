import logging
import os

from django.core.management.base import BaseCommand

from annotation.annotation_version_querysets import get_unannotated_variants_qs
from annotation.models.models import AnnotationVersion
from annotation.signals import annotation_run_complete_signal
from snpdb.models.models_genome import GenomeBuild


def set_annotation_version_complete(annotation_version, force):
    logging.info("Setting %s to complete", annotation_version)
    # This is a brute force check - as annotation range lock can sometimes miss some variants
    # due to them being skipped by VEP
    unannotated_variants_qs = get_unannotated_variants_qs(annotation_version=annotation_version)
    unannotated_count = unannotated_variants_qs.count()
    if unannotated_count:
        msg = f"There are {unannotated_count} unannotated variants! "
        msg += "You need to call 'annotate_unannotated_variants'"
        if force:
            print(f"Warning: {msg}")
        else:
            raise ValueError(msg)

    annotation_run_complete_signal.send(sender=os.path.basename(__file__),
                                        variant_annotation_version=annotation_version.variant_annotation_version)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--annotation-version', type=int, help='Annotation Version (defaults to latest)')
        parser.add_argument('--force', action='store_true', help='Force even if unannotated variants detected')

    def handle(self, *args, **options):
        annotation_version_id = options.get("annotation_version")
        force = options.get("force")
        if annotation_version_id:
            annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
            set_annotation_version_complete(annotation_version, force=force)
        else:
            for genome_build in GenomeBuild.builds_with_annotation():
                annotation_version = AnnotationVersion.latest(genome_build)
                set_annotation_version_complete(annotation_version, force=force)
