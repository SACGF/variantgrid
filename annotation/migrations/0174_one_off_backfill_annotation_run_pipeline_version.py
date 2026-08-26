from django.db import migrations

ACTIVE = "ACTIVE"  # AnnotationPipelineVersion.Status.ACTIVE


def _backfill_run_pipeline_version(apps, _schema_editor):
    """ #720: runs of a versioned pipeline created before AnnotationPipelineVersion existed have no
        version to check the installed tool against, so annotate() dies on the None. They were scheduled
        against whatever was installed at the time, which is what ACTIVE names now - stamp them with it
        so they run (and so the scheduler sees the lock as covered rather than making a second run). """
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
        ("annotation", "0173_annotationrun_dispatch_count"),
    ]

    operations = [
        migrations.RunPython(_backfill_run_pipeline_version, migrations.RunPython.noop),
    ]
