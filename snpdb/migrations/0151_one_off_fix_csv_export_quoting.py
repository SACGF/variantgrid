# Generated by Django 4.2.9 on 2024-10-08 06:29

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _delete_failed_csv_exports(apps, _schema_editor):
    CachedGeneratedFile = apps.get_model("snpdb", "CachedGeneratedFile")
    AFFECTED_GENERATORS = ["export_sample_to_downloadable_file", "export_cohort_to_downloadable_file"]
    CachedGeneratedFile.objects.filter(generator__in=AFFECTED_GENERATORS, exception__contains='NoneType').delete()


def _test_for_csv_export(apps):
    CachedGeneratedFile = apps.get_model("snpdb", "CachedGeneratedFile")

    AFFECTED_GENERATORS = ["export_sample_to_downloadable_file", "export_cohort_to_downloadable_file"]
    return CachedGeneratedFile.objects.filter(generator__in=AFFECTED_GENERATORS, filename__icontains='.csv').exists()


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0150_rename_cachedgeneratedfile2_cachedgeneratedfile"),
    ]

    operations = [
        migrations.RunPython(_delete_failed_csv_exports),
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_fix_csv_export_quoting"]),
                        test=_test_for_csv_export)
    ]
