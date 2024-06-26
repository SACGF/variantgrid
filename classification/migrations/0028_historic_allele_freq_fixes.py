# Generated by Django 3.1.6 on 2021-04-27 05:31

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _check_has_classifications_with_allele_frequency(apps):
    Classification = apps.get_model("classification", "Classification")
    return Classification.objects.filter(evidence__icontains='allele_frequency').exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0027_new_ekey_sample_date'),
    ]

    operations = [
        ManualOperation.operation_manage(["fix_allele_frequency_history"],
                                         test=_check_has_classifications_with_allele_frequency)
    ]
