from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_abandoned_annotation_runs(apps):
    """ Unfinished runs with no range lock, or on a HISTORICAL version - @see https://github.com/SACGF/variantgrid/issues/1940 """
    AnnotationRun = apps.get_model("annotation", "AnnotationRun")
    unfinished_qs = AnnotationRun.objects.exclude(status="F")  # AnnotationStatus.FINISHED
    return unfinished_qs.filter(annotation_range_lock__isnull=True).exists() or \
        unfinished_qs.filter(annotation_range_lock__version__status="HISTORICAL").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0198_cached_web_resource_to_snpdb"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["clean_up_abandoned_annotation_runs"]),
                        note="Delete annotation runs nothing will dispatch (no range lock, or a HISTORICAL "
                             "version) and their empty range locks (#1940)",
                        test=_has_abandoned_annotation_runs),
    ]
