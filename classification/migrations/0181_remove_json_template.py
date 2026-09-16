"""

Drop ClassificationReportTemplate.json_template (#444).

A report's JSON is an interface to another system rather than a document a lab restyles, so it is
built in Python by the app that keeps it in step with that system - see
classification/report/renderers.py:render_json. Nothing ever had a json_template: the shipped default
was blank (the canonical context dump), which is still what a report gets when no app answers.

"""
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("classification", "0180_fold_change_ekey"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="classificationreporttemplate",
            name="json_template",
        ),
    ]
