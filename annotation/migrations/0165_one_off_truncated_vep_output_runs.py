from django.db import migrations

from annotation.vep_warnings import annotation_run_shortfall
from manual.operations.manual_operations import ManualOperation

FINISHED = 'F'  # AnnotationStatus.FINISHED


def _has_truncated_annotation_runs(apps):
    """ Runs that finished having imported fewer records than were handed to VEP, without VEP reporting
        the difference as skipped - @see https://github.com/SACGF/variantgrid/issues/1701 """
    AnnotationRun = apps.get_model("annotation", "AnnotationRun")
    qs = AnnotationRun.objects.filter(status=FINISHED,
                                      dump_count__isnull=False,
                                      annotated_count__isnull=False)
    for dump_count, annotated_count, vep_warnings in qs.values_list("dump_count", "annotated_count",
                                                                    "vep_warnings").iterator():
        if annotation_run_shortfall(dump_count, annotated_count, vep_warnings) > 0:
            return True
    return False


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0164_annotationrun_vep_skipped_variants_filename"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_truncated_annotation_runs", "--mark-error"]),
                        note="Mark AnnotationRuns whose VEP output was silently truncated as ERROR so they "
                             "can be retried - see #1701",
                        test=_has_truncated_annotation_runs),
    ]
