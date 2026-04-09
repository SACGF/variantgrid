from django.db import migrations


def rename_institution_to_organisation(apps, schema_editor):
    Classification = apps.get_model('classification', 'Classification')
    ClassificationModification = apps.get_model('classification', 'ClassificationModification')

    Classification.objects.filter(share_level='institution').update(share_level='organisation')
    ClassificationModification.objects.filter(share_level='institution').update(share_level='organisation')


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0161_importedalleleinfo_hgvs_converter_data_version_and_more'),
    ]

    operations = [
        migrations.RunPython(rename_institution_to_organisation, migrations.RunPython.noop),
    ]
