from django.db import migrations


def _set_lab_share_level(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key="curation_system_migration_index", max_share_level="private").update(max_share_level="lab")


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0228_importedalleleinfo_transcript_genome_ref_differ'),
    ]

    operations = [
        migrations.RunPython(_set_lab_share_level, migrations.RunPython.noop)
    ]
