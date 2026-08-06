from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_columns_version_5_frameshift_rows(apps):
    """ Only deployments with a columns_version=5 partition carrying frameshift rows have
        anything to backfill - @see https://github.com/SACGF/variantgrid/issues/579 """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")
    vav_ids = VariantAnnotationVersion.objects.filter(columns_version=5).values_list("pk", flat=True)
    return VariantAnnotation.objects.filter(version_id__in=vav_ids,
                                            consequence__contains="frameshift_variant").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0162_ptc_annotation"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["backfill_ptc_annotation"]),
                        note="Populate the predicted stop codon position and PTC-aware NMD prediction "
                             "on existing columns_version=5 annotation (#579)",
                        test=_has_columns_version_5_frameshift_rows),
    ]
