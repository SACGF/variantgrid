# Generated by Django 4.2.2 on 2023-08-14 23:42
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _override_blank_condition(apps):
    Classification = apps.get_model("classification", "Classification")
    condition_check = Classification.objects.filter(
            condition_resolution__isnull=False,
            evidence__condition__validation__icontains="error")
    return condition_check.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0108_discordancereporttriage'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["condition_matching_to_override_blank_condition"]),
                        test=_override_blank_condition)
    ]