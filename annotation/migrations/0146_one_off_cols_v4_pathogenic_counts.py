from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _check_has_columns_version4(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version=4).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0145_clinvarversion_data_archive_reason_and_more'),
    ]

    operations = [
        ManualOperation.operation_manage(["fix_columns_version4_damage_counts"],
                                         test=_check_has_columns_version4)
    ]
