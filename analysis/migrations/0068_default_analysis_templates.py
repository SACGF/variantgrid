# Generated by Django 4.1.3 on 2022-11-29 02:52

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_no_templates(apps):
    AnalysisTemplate = apps.get_model("analysis", "AnalysisTemplate")
    return not AnalysisTemplate.objects.all().exists()


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0067_alter_analysis_options'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["analysis_create_default_templates"]),
                        test=_test_no_templates)
    ]
