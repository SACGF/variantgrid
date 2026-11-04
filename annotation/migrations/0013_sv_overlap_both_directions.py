from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_variant_annotation_versions(apps):
    return apps.get_model("annotation", "VariantAnnotationVersion").objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0012_vep_settings_snapshots'),
        ('manual', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='variantannotationversion',
            name='sv_overlap_reciprocal',
        ),
        ManualOperation.operation_manage(["one_off_sv_overlap_both_directions"],
                                         test=_test_has_variant_annotation_versions),
    ]
