# Generated by Django 3.2.4 on 2021-09-22 01:35

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_clingen_alleles(apps):
    ClinGenAllele = apps.get_model("snpdb", "ClinGenAllele")
    return ClinGenAllele.objects.all().exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0047_rename_behalf_org_id_clinvarkey_org_id'),
    ]

    operations = [
        ManualOperation(
            task_id=ManualOperation.task_id_manage(["clingen_allele_replace", "--indels"]),
            test=_test_has_clingen_alleles),
    ]
