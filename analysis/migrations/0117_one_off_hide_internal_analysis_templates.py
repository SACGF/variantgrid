from django.conf import settings
from django.db import migrations


def _hide_internal_analysis_templates(apps, _schema_editor):
    """ Internal auto-analysis templates (used for sample/cohort grid exports) shouldn't appear
        on the analysis templates page. AnalysisTemplate permissions delegate to its Analysis,
        so marking the analysis invisible hides the template from filter_for_user() """
    AnalysisTemplate = apps.get_model("analysis", "AnalysisTemplate")
    AnalysisTemplateVersion = apps.get_model("analysis", "AnalysisTemplateVersion")
    Analysis = apps.get_model("analysis", "Analysis")

    hidden_template_names = []
    for setting_name in ["ANALYSIS_TEMPLATES_AUTO_SAMPLE", "ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT"]:
        if template_name := getattr(settings, setting_name, None):
            hidden_template_names.append(template_name)

    hidden_templates = AnalysisTemplate.objects.filter(name__in=hidden_template_names)
    analysis_ids = hidden_templates.filter(analysis__isnull=False).values_list("analysis_id", flat=True)
    Analysis.objects.filter(pk__in=list(analysis_ids)).update(visible=False)
    AnalysisTemplateVersion.objects.filter(template__in=hidden_templates).update(appears_in_autocomplete=False)


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0116_classificationsnode_allele_origin_and_more'),
    ]

    operations = [
        migrations.RunPython(_hide_internal_analysis_templates, migrations.RunPython.noop),
    ]
