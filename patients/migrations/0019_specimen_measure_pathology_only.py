from django.db import migrations, models


def _delete_sequencing_measures(apps, _schema_editor):
    """ TMB, MSI, GIS, ploidy and the caller's tumour fraction are now the analysis' own record
        (seqauto.DragenTSO500CombinedVariantOutput), recreated by re-sending the CombinedVariantOutputs -
        the old uploads carry no 'sequencing_run', which keys that row. What stays is the pathologist's
        tumour content (Mocha) """
    SpecimenMeasure = apps.get_model("patients", "SpecimenMeasure")
    SpecimenMeasure.objects.exclude(measure_type="F").delete()
    SpecimenMeasure.objects.filter(method__startswith="DRAGEN").delete()


class Migration(migrations.Migration):

    dependencies = [
        ('patients', '0018_alter_patient_last_name'),
    ]

    operations = [
        migrations.RunPython(_delete_sequencing_measures, migrations.RunPython.noop),
        migrations.AlterField(
            model_name='specimenmeasure',
            name='measure_type',
            field=models.CharField(choices=[('F', 'Tumour content (pathology)')], max_length=1),
        ),
        # Only the DRAGEN import wrote these
        migrations.RemoveField(
            model_name='specimenmeasure',
            name='extraction',
        ),
        migrations.RemoveField(
            model_name='specimenmeasure',
            name='measured_date',
        ),
        migrations.RemoveField(
            model_name='specimenmeasure',
            name='threshold',
        ),
        migrations.RemoveField(
            model_name='specimenmeasure',
            name='threshold_source',
        ),
    ]
