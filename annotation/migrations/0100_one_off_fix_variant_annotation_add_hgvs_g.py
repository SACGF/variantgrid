# Generated by Django 4.2.10 on 2024-04-03 05:42

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _check_missing_hgvs_g(apps):
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")

    if VariantAnnotation.objects.filter(variant_class='SN', hgvs_g__isnull=True).exists():
        return True

    return False


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0099_variantannotation_hgvs_g'),
    ]

    operations = [
        ManualOperation.operation_manage(["fix_variant_annotation_add_hgvs_g"],
                                         test=_check_missing_hgvs_g)
    ]