# Generated by Django 4.2.9 on 2024-11-12 04:11
# This can be skipped if you bulk insert EKeys from DB copy

from django.db import migrations

def _ekey_allele_origin_germline_somatic_only(apps, _schema_editor):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')
    EvidenceKey.objects.filter(pk='allele_origin').update(
        mandatory=True,
        description="Indicate whether variant is Germline (heritable) or Somatic (non-heritable). To specify confidence, set field 'allele_origin_confirmation' to 'confirmed' or 'unconfirmed'",
        options = [
            {'key': 'germline',
            'index': 1,
            'clinvar': 'germline',
            'namespaces': ['germline']},
            {'key': 'somatic',
            'index': 2,
            'clinvar': 'somatic',
            'namespaces': ['somatic']},
        ]
    )



class Migration(migrations.Migration):

    dependencies = [
        ("classification", "0154_alter_importedalleleinfo_status"),
    ]

    operations = [
        migrations.RunPython(_ekey_allele_origin_germline_somatic_only),
    ]