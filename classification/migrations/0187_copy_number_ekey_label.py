from django.db import migrations


def _set_copy_number_label(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(pk="copy_number").update(label="Copy number")


def _reverse_copy_number_label(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(pk="copy_number").update(label="CN")


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0186_one_off_rematch_failed_gene_level_allele_infos'),
    ]

    operations = [
        migrations.RunPython(_set_copy_number_label, _reverse_copy_number_label),
    ]
