"""
#1704 - Specimen -> Extraction -> Sample.

Specimen swaps its TextField primary key for a surrogate one (reference_id becomes unique per
patient rather than globally), gains ExternallyManagedModel's external_pk, and hands
nucleic_acid_source down to the new Extraction. Every existing Specimen gets one Extraction
carrying whatever nucleic_acid_source it had.

FKs into Specimen are dropped and rebuilt around the key change - snpdb.0201 does that for Sample,
this does it for PatientRecord.
"""
import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models
from django.db.models import F
from django.utils import timezone


def _flush_deferred_constraints(schema_editor):
    """ Postgres refuses ALTER TABLE while deferred FK trigger events from a data change are pending """
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("SET CONSTRAINTS ALL IMMEDIATE")


def _stash_specimen_reference(apps, schema_editor):
    PatientRecord = apps.get_model("patients", "PatientRecord")
    PatientRecord.objects.filter(specimen__isnull=False).update(_specimen_reference_id=F("specimen_id"))
    _flush_deferred_constraints(schema_editor)


def _create_extractions(apps, schema_editor):
    Specimen = apps.get_model("patients", "Specimen")
    Extraction = apps.get_model("patients", "Extraction")
    extractions = [Extraction(specimen=specimen, nucleic_acid_source=specimen.nucleic_acid_source)
                   for specimen in Specimen.objects.all().iterator()]
    Extraction.objects.bulk_create(extractions, batch_size=2000)
    _flush_deferred_constraints(schema_editor)


def _restore_specimen_reference(apps, schema_editor):
    """ reference_id was globally unique while it was the primary key, so it identifies the row """
    PatientRecord = apps.get_model("patients", "PatientRecord")
    Specimen = apps.get_model("patients", "Specimen")
    specimen_by_reference_id = dict(Specimen.objects.values_list("reference_id", "pk"))

    patient_records = []
    for patient_record in PatientRecord.objects.filter(_specimen_reference_id__isnull=False).iterator():
        patient_record.specimen_id = specimen_by_reference_id.get(patient_record._specimen_reference_id)
        patient_records.append(patient_record)
    PatientRecord.objects.bulk_update(patient_records, ["specimen_id"], batch_size=2000)
    _flush_deferred_constraints(schema_editor)


class Migration(migrations.Migration):

    dependencies = [
        ('patients', '0012_manual_patient_code_backfill_reminder'),
        ('snpdb', '0201_sample_stash_specimen_reference'),
    ]

    operations = [
        migrations.AddField(
            model_name='patientrecord',
            name='_specimen_reference_id',
            field=models.TextField(null=True),
        ),
        migrations.RunPython(_stash_specimen_reference, migrations.RunPython.noop, elidable=True),
        migrations.RemoveField(
            model_name='patientrecord',
            name='specimen',
        ),

        # Nothing points at Specimen.reference_id now, so the key can be swapped
        migrations.AlterField(
            model_name='specimen',
            name='reference_id',
            field=models.TextField(),
        ),
        migrations.AddField(
            model_name='specimen',
            name='id',
            field=models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID'),
            preserve_default=False,
        ),
        migrations.AlterUniqueTogether(
            name='specimen',
            unique_together={('patient', 'reference_id')},
        ),
        migrations.AddField(
            model_name='specimen',
            name='created',
            field=django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, default=timezone.now,
                                                                    verbose_name='created'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='specimen',
            name='modified',
            field=django_extensions.db.fields.ModificationDateTimeField(auto_now=True, default=timezone.now,
                                                                        verbose_name='modified'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='specimen',
            name='external_pk',
            field=models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE,
                                       to='patients.externalpk'),
        ),

        migrations.CreateModel(
            name='Extraction',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True,
                                                                              verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True,
                                                                                   verbose_name='modified')),
                ('nucleic_acid_source', models.CharField(blank=True,
                                                         choices=[('D', 'DNA'), ('R', 'RNA')],
                                                         default='D', max_length=1, null=True)),
                ('external_pk', models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE,
                                                     to='patients.externalpk')),
                ('specimen', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE,
                                               to='patients.specimen')),
            ],
        ),
        migrations.RunPython(_create_extractions, migrations.RunPython.noop, elidable=True),
        migrations.RemoveField(
            model_name='specimen',
            name='nucleic_acid_source',
        ),

        migrations.AddField(
            model_name='patientrecord',
            name='specimen',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL,
                                    related_name='related_specimen', to='patients.specimen'),
        ),
        migrations.RunPython(_restore_specimen_reference, migrations.RunPython.noop, elidable=True),
        migrations.RemoveField(
            model_name='patientrecord',
            name='_specimen_reference_id',
        ),
    ]
