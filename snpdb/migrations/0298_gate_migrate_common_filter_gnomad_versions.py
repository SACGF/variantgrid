from django.db import migrations

from manual.operations.manual_operations import ManualOperation

GNOMAD_TASK_ID = ManualOperation.task_id_manage(["migrate_common_filter_gnomad_versions"])


def _gate_on_annotation(apps, schema_editor):
    """ The command fails without annotation: gated, the menu shows it blocked rather than it dying partway.
        Set here rather than in 0188, which deployments have already applied """
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationTask.objects.filter(pk=GNOMAD_TASK_ID).update(requires=["annotation-version"])


def _clear_gate(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationTask.objects.filter(pk=GNOMAD_TASK_ID).update(requires=[])


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0297_retire_somalier_download_manual_task"),
        ("manual", "0008_systemd_env_files"),
    ]

    operations = [
        migrations.RunPython(_gate_on_annotation, _clear_gate),
    ]
