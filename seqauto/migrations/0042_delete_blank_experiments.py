from django.db import migrations


def delete_blank_experiments(apps, schema_editor):
    """ The old ExperimentManager stripped "RPT" off names, so eg "RPT" was stored as "".
        SequencingRun.experiment is SET_NULL so nothing else is lost """
    Experiment = apps.get_model("seqauto", "Experiment")
    Experiment.objects.filter(name="").delete()


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0041_remove_seqautorecord_data_state_and_more"),
    ]

    operations = [
        migrations.RunPython(delete_blank_experiments, migrations.RunPython.noop),
    ]
