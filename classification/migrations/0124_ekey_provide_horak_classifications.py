# Generated by Django 4.2.5 on 2024-01-05 02:00

from django.db import migrations

from classification.evidence_key_rename import OptionUpdator


def _update_clinical_significance_to_classification(apps, _schema):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')
    clin_sig = EvidenceKey.objects.get(key='clinical_significance')

    option_updator = OptionUpdator(clin_sig)
    option_updator.set_attributes("LP", namespace="acmg")
    option_updator.set_attributes("P", namespace="acmg")
    option_updator.set_attributes("R", namespace="germline")
    option_updator.set_attributes("D", namespace="germline")
    option_updator.ensure_option(
    {
             "vg": "4",
             "key": "LO",
             "index": 7,
             "label": "Likely Oncogenic",
             "bucket": 3,
             "clinvar": "Likely oncogenic",
             "namespace": "horak"
         }
    )
    option_updator.ensure_option(
    {
            "vg": "5",
            "key": "O",
            "index": 8,
            "label": "Oncogenic",
            "bucket": 3,
            "clinvar": "Oncogenic",
            "namespace": "horak"
        }
    )
    option_updator.preferred_order([
        "B", "LB", "VUS", "VUS_A", "VUS_B", "VUS_C", "LP", "P", "LO", "O", "D", "R", "A"
    ])
    option_updator.save()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0123_auto_20231205_1110'),
    ]

    operations = [
        migrations.RunPython(_update_clinical_significance_to_classification)
    ]
