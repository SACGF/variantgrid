from django.db import migrations

from classification.models.classification_report_models import ReportNames
from classification.report.default_templates import generic_case_template, generic_json_template


def _create_default_case_report_template(apps, _schema_editor):
    """ The generic pair a deployment gets before it writes its own - get_or_create, so a lab that
        has already edited the row keeps what it wrote """
    ClassificationReportTemplate = apps.get_model("classification", "ClassificationReportTemplate")
    ClassificationReportTemplate.objects.get_or_create(
        name=ReportNames.DEFAULT_CASE_REPORT,
        defaults={
            "case_template": generic_case_template(),
            "json_template": generic_json_template(),
        })


def _delete_default_case_report_template(apps, _schema_editor):
    ClassificationReportTemplate = apps.get_model("classification", "ClassificationReportTemplate")
    ClassificationReportTemplate.objects.filter(name=ReportNames.DEFAULT_CASE_REPORT).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("classification", "0178_case_report"),
    ]

    operations = [
        migrations.RunPython(_create_default_case_report_template,
                             reverse_code=_delete_default_case_report_template),
    ]
