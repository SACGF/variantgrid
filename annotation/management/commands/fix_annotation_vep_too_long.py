"""
Backfill vep_skipped_reason=TOO_LONG rows for SVs over VEP's size cap (#2104).

The dump leaves those SVs out, and before #2104 their rows were only written when VEP happened to skip
something else in the same run - so most finished ranges hold long SVs with no VariantAnnotation, which the
variant page reports as a failed annotation. Works one finished STRUCTURAL_VARIANT run (= one range lock's
pk range) at a time, probing snpdb_variant for a long SV in range before touching the annotation partition.
Rows are written through annotation/vcf_files/import_vcf_annotations.py:insert_vep_too_long_skipped, gene
overlaps included, so fix_annotation_sv_overlaps has nothing to add for them.
"""
import logging

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.annotation_version_querysets import filter_vep_sv_max_size
from annotation.models.models import AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.vcf_files.import_vcf_annotations import insert_vep_too_long_skipped
from snpdb.archive import DataArchivedError
from snpdb.models import Variant

# Import processing files are named by batch - keep clear of any the run's original import left behind
BACKFILL_BATCH_ID = 2104


def fix_annotation_vep_too_long():
    for vav in VariantAnnotationVersion.objects.order_by("pk"):
        sv_max_size = vav.sv_max_size or settings.ANNOTATION_VEP_SV_MAX_SIZE
        if not sv_max_size or vav.get_any_annotation_version() is None:
            continue
        try:
            _fix_variant_annotation_version(vav, sv_max_size)
        except DataArchivedError as e:
            logging.info("Skipping %s: %s", vav, e)


def _fix_variant_annotation_version(vav: VariantAnnotationVersion, sv_max_size: int):
    sv_runs_qs = AnnotationRun.objects.filter(
        annotation_range_lock__version=vav,
        pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
        status=AnnotationStatus.FINISHED,
    ).select_related("annotation_range_lock").order_by("annotation_range_lock__min_variant_id")
    num_runs = 0
    num_rows = 0
    for annotation_run in sv_runs_qs:
        range_lock = annotation_run.annotation_range_lock
        long_sv_qs = Variant.objects.filter(pk__gte=range_lock.min_variant_id,
                                            pk__lte=range_lock.max_variant_id,
                                            svlen__isnull=False)
        # pk-ordered first() rather than exists() so the planner only reads the range
        if filter_vep_sv_max_size(long_sv_qs, sv_max_size, too_long=True).order_by("pk").first() is None:
            continue
        if inserted := insert_vep_too_long_skipped(annotation_run, sv_max_size=sv_max_size,
                                                   batch_id=BACKFILL_BATCH_ID):
            num_runs += 1
            num_rows += inserted
            logging.info("%s: wrote %d TOO_LONG rows", annotation_run, inserted)
    logging.info("%s: wrote %d TOO_LONG rows over %d runs", vav, num_rows, num_runs)


class Command(BaseCommand):
    """ Idempotent - only variants still without a VariantAnnotation row get one. """
    category = "one-off"

    def handle(self, *args, **options):
        fix_annotation_vep_too_long()
