from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_un_backfilled_columns_version_5(apps):
    """ Only deployments with a columns_version=5 partition still flagged as un-backfilled have
        anything to do - @see https://github.com/SACGF/variantgrid/issues/579 """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version=5, backfilled_ptc=False).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0162_ptc_annotation"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["backfill_ptc_annotation"]),
                        note="Populate the predicted stop codon position and PTC-aware NMD prediction "
                             "on existing columns_version=5 annotation (#579)",
                        test=_has_un_backfilled_columns_version_5),
    ]
