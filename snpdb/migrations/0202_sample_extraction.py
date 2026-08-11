"""
#1704 - Sample.extraction replaces Sample.specimen, specimen now being reached via
sample.extraction.specimen. patients.0013 made one Extraction per existing Specimen, so each sample
re-attaches to the extraction under the specimen it was pointing at.
"""
import django.db.models.deletion
from django.db import migrations, models


def _flush_deferred_constraints(schema_editor):
    """ Postgres refuses ALTER TABLE while deferred FK trigger events from a data change are pending """
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("SET CONSTRAINTS ALL IMMEDIATE")


def _restore_sample_extraction(apps, schema_editor):
    """ reference_id was globally unique while it was Specimen's primary key, so it identifies the row """
    Sample = apps.get_model("snpdb", "Sample")
    Extraction = apps.get_model("patients", "Extraction")
    extraction_by_reference_id = dict(Extraction.objects.values_list("specimen__reference_id", "pk"))

    samples = []
    for sample in Sample.objects.filter(_specimen_reference_id__isnull=False).iterator():
        sample.extraction_id = extraction_by_reference_id.get(sample._specimen_reference_id)
        samples.append(sample)
    Sample.objects.bulk_update(samples, ["extraction_id"], batch_size=2000)
    _flush_deferred_constraints(schema_editor)


class Migration(migrations.Migration):

    dependencies = [
        ('patients', '0013_extraction_specimen_surrogate_pk'),
        ('snpdb', '0201_sample_stash_specimen_reference'),
    ]

    operations = [
        migrations.AddField(
            model_name='sample',
            name='extraction',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL,
                                    to='patients.extraction'),
        ),
        migrations.RunPython(_restore_sample_extraction, migrations.RunPython.noop, elidable=True),
        migrations.RemoveField(
            model_name='sample',
            name='_specimen_reference_id',
        ),
    ]
