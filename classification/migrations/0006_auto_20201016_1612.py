# Generated by Django 3.1 on 2020-10-16 05:42

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0005_auto_20201001_1706'),
        ('manual', '0002_deployment'),
    ]

    operations = [
        # These flags are no longer used, instead they are handled by ImportedAlleleInfo
        # ManualOperation(task_id=ManualOperation.task_id_manage(["fix_matching_flags", "--apply"]),
        #                 test=_check_has_flags_of_interest)
    ]
