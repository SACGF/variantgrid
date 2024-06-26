# Generated by Django 3.1 on 2020-11-23 02:24

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def test_has_sample_stats(apps):
    SampleVariantAnnotationStats = apps.get_model("annotation", "SampleVariantAnnotationStats")
    return SampleVariantAnnotationStats.objects.all().exists()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0009_vcfannotationstats'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["calculate_sample_stats", "--clear"]),
                        test=test_has_sample_stats)
    ]
