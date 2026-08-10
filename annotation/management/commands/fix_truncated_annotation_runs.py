"""
Find AnnotationRuns whose VEP output was silently truncated, and mark them ERROR so they can be
retried (#1701).

VEP could exit 0 having written a short annotated VCF, and the import lane trusted that exit code: the
missing records became VariantAnnotation rows carrying vep_skipped_reason=UNKNOWN and the run reached
FINISHED. Those rows make the variants look annotated, so they are never offered for annotation again.

Detection is pure arithmetic over columns already on the row - records dumped versus records imported
plus records VEP said it skipped - so no VEP output needs to still be on disk.
"""
from django.core.management.base import BaseCommand

from annotation.models import AnnotationRun
from annotation.models.models_enums import AnnotationStatus
from annotation.vep_warnings import annotation_run_shortfall, parse_vep_warnings

ERROR_EXCEPTION_FORMAT = (
    "#1701: VEP output truncated - dump_count={dump_count}, annotated_count={annotated_count}, "
    "vep warned {warned_count} (shortfall {shortfall}). Full retry required."
)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--mark-error', action='store_true',
                            help="Write a descriptive error_exception on each affected run, moving it to "
                                 "ERROR so it can be retried from the admin UI")
        parser.add_argument('--min-shortfall', type=int, default=1,
                            help="Only report runs missing at least this many records (default: 1)")

    def handle(self, *args, **options):
        mark_error = options["mark_error"]
        min_shortfall = options["min_shortfall"]

        qs = AnnotationRun.objects.filter(status=AnnotationStatus.FINISHED,
                                          dump_count__isnull=False,
                                          annotated_count__isnull=False)
        qs = qs.select_related("annotation_range_lock__version").order_by("pk")

        header = ("run", "variant_annotation_version", "pipeline", "dump_count",
                  "annotated_count", "warned", "shortfall")
        self.stdout.write("\t".join(header))

        total_runs = 0
        total_variants = 0
        for annotation_run in qs:
            shortfall = annotation_run_shortfall(annotation_run.dump_count,
                                                 annotation_run.annotated_count,
                                                 annotation_run.vep_warnings)
            if shortfall < min_shortfall:
                continue
            warned_count = parse_vep_warnings(annotation_run.vep_warnings).skipped_count
            total_runs += 1
            total_variants += shortfall
            self.stdout.write("\t".join(str(f) for f in (
                annotation_run.pk,
                annotation_run.annotation_range_lock.version_id,
                annotation_run.pipeline_type,
                annotation_run.dump_count,
                annotation_run.annotated_count,
                warned_count,
                shortfall,
            )))

            if mark_error:
                annotation_run.error_exception = ERROR_EXCEPTION_FORMAT.format(
                    dump_count=annotation_run.dump_count,
                    annotated_count=annotation_run.annotated_count,
                    warned_count=warned_count,
                    shortfall=shortfall,
                )
                annotation_run.save()  # get_status() -> ERROR

        self.stdout.write(f"{total_runs} affected runs, {total_variants} variants missing annotation")
        if not total_runs:
            return

        if mark_error:
            self.stdout.write("Marked ERROR. Re-annotate from the admin UI - the per-run retry button, or "
                              "'Retry all failed runs' for a whole variant annotation version.")
            self.stdout.write("'Retry all failed runs' retries *every* ERROR run for that version, so any "
                              "unrelated failures go with them - use the per-run button where that matters.")
            self.stdout.write("Re-annotation is slow (roughly half an hour of VEP per run), and concurrent "
                              "VEP memory pressure is what caused this - retry in batches.")
        else:
            self.stdout.write("Re-run with --mark-error to move these to ERROR so they can be retried.")
