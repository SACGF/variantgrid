# Generated by Django 4.0.7 on 2023-02-06 05:48

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _check_has_classifications(apps):
    Classification = apps.get_model("classification", "Classification")
    return Classification.objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0090_remove_classification_classification_import_and_more'),
    ]

    operations = [
        ManualOperation(
            # sort order previously put X after 9 and before 10 (instead of after 22)
            task_id=ManualOperation.task_id_manage(["fix_variant_matching", "--sort"]),
            test=_check_has_classifications
        )
    ]
