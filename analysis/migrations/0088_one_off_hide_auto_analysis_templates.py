# Generated by Django 4.2.11 on 2024-11-27 00:00
from django.conf import settings
from django.db import migrations


def _one_off_hide_auto_analysis_templates(apps, _schema_editor):
    # This is now done in analysis_create_default_templates - so can be deleted
    AnalysisTemplateVersion = apps.get_model("analysis", "AnalysisTemplateVersion")

    HIDDEN_TEMPLATES = [
        settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE,
        settings.ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT,
    ]

    AnalysisTemplateVersion.objects.filter(template__name__in=HIDDEN_TEMPLATES).update(appears_in_autocomplete=False)


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0087_alter_karyomappinggene_unique_together"),
    ]

    operations = [
        migrations.RunPython(_one_off_hide_auto_analysis_templates)
    ]
