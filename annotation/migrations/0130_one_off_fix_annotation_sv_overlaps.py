from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _check_has_long_sv_annotation(apps):
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")

    VEP_SKIPPED_TOO_LONG = 'l'

    qs = VariantAnnotation.objects.filter(vep_skipped_reason=VEP_SKIPPED_TOO_LONG,
                                          overlapping_symbols__isnull=True)
    return qs.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0129_delete_columnvepfield'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_annotation_sv_overlaps"]),
                        test=_check_has_long_sv_annotation)
    ]
