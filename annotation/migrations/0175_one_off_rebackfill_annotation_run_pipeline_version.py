from django.db import migrations

ACTIVE = "ACTIVE"  # AnnotationPipelineVersion.Status.ACTIVE


def _backfill_run_pipeline_version(apps, _schema_editor):
    """ #720: 0174 stamped these once, then reset_for_retry rebuilt the row without carrying
        pipeline_version over, so retrying an errored run put it straight back to NULL and
        check_tool_version refused to run it. reset_for_retry now keeps the version - this puts the runs
        it already cleared back onto ACTIVE so they can be retried. """
    AnnotationRun = apps.get_model("annotation", "AnnotationRun")
    AnnotationPipelineVersion = apps.get_model("annotation", "AnnotationPipelineVersion")

    for pipeline_version in AnnotationPipelineVersion.objects.filter(status=ACTIVE):
        unversioned_qs = AnnotationRun.objects.filter(
            pipeline_type=pipeline_version.pipeline_type,
            pipeline_version__isnull=True,
            annotation_range_lock__version__genome_build=pipeline_version.genome_build_id,
        )
        # A lock that already has a run at this version keeps it - stamping a second one onto it would
        # break one_versioned_run_per_range_lock_type_and_version
        already_at_version = AnnotationRun.objects.filter(
            pipeline_type=pipeline_version.pipeline_type,
            pipeline_version=pipeline_version,
        ).values_list("annotation_range_lock_id", flat=True)
        unversioned_qs.exclude(annotation_range_lock_id__in=already_at_version).update(
            pipeline_version=pipeline_version)


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0174_one_off_backfill_annotation_run_pipeline_version"),
    ]

    operations = [
        migrations.RunPython(_backfill_run_pipeline_version, migrations.RunPython.noop),
    ]
