# Generated by Django 3.1 on 2021-01-29 00:15
from datetime import datetime

from django.db import migrations
from django.utils import timezone

from manual.operations.manual_operations import ManualOperation


def test_need_hgnc(apps):
    GeneVersion = apps.get_model("genes", "GeneVersion")
    HGNCGeneNames = apps.get_model("genes", "HGNCGeneNames")
    HGNCGeneNamesImport = apps.get_model("genes", "HGNCGeneNamesImport")

    if not HGNCGeneNames.objects.exists():
        return False

    fixed_hgnc_date = timezone.make_aware(datetime(2021, 1, 29))
    no_new_import = not HGNCGeneNamesImport.objects.filter(created__gte=fixed_hgnc_date).exists()
    gv_but_no_fixes = GeneVersion.objects.exists() and not GeneVersion.objects.filter(hgnc__isnull=False).exists()
    return no_new_import or gv_but_no_fixes


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0014_auto_20210120_1723'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_hgnc"]), test=test_need_hgnc)
    ]