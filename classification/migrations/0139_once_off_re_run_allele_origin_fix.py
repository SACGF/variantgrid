# Generated by Django 4.2.10 on 2024-05-29 04:26

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0138_auto_20240527_1154'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(
            ["classification_set_allele_context"]
        ))
    ]
