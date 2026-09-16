from django.db import migrations

_DESCRIPTION = ('ClinVar and internal classifications, germline pair then somatic pair. Hover a '
                'chip for the per-record summary, expand the row for the full ClinVar record and '
                'the classifications behind it')
_PREVIOUS_DESCRIPTION = ('ClinVar and internal classifications, germline pair then somatic pair. Hover a '
                         'chip for the per-record summary, click it to scroll to the underlying column')


def _set_description(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="classifications").update(description=_DESCRIPTION)


def _restore_description(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="classifications").update(description=_PREVIOUS_DESCRIPTION)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0244_open_targets_backfill_note_link'),
    ]

    operations = [
        migrations.RunPython(_set_description, reverse_code=_restore_description),
    ]
