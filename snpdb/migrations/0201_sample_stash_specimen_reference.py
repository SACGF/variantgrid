"""
Specimen's TextField primary key is about to be replaced with a surrogate key (#1704), so drop the
FK pointing at it, keeping the reference_id text in a temporary column. 0202 uses that to re-attach
each sample to the Extraction created under its old specimen.
"""
from django.db import migrations, models
from django.db.models import F


def _flush_deferred_constraints(schema_editor):
    """ Postgres refuses ALTER TABLE while deferred FK trigger events from a data change are pending """
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("SET CONSTRAINTS ALL IMMEDIATE")


def _stash_specimen_reference(apps, schema_editor):
    Sample = apps.get_model("snpdb", "Sample")
    Sample.objects.filter(specimen__isnull=False).update(_specimen_reference_id=F("specimen_id"))
    _flush_deferred_constraints(schema_editor)


class Migration(migrations.Migration):

    dependencies = [
        ('patients', '0012_manual_patient_code_backfill_reminder'),
        ('snpdb', '0200_ptc_columns'),
    ]

    operations = [
        migrations.AddField(
            model_name='sample',
            name='_specimen_reference_id',
            field=models.TextField(null=True),
        ),
        migrations.RunPython(_stash_specimen_reference, migrations.RunPython.noop, elidable=True),
        migrations.RemoveField(
            model_name='sample',
            name='specimen',
        ),
    ]
