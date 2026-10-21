from django.db import migrations

from manual.operations.manual_operations import ManualOperation

RETIRED_TASK_ID = ManualOperation.task_id_manage(["download_annotation_data"])


def _complete_retired_task(apps, schema_editor):
    """ deployment_check's somalier_allele_order fails until the hosted somalier build is installed (#183) """
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    if task := ManualMigrationTask.objects.filter(pk=RETIRED_TASK_ID).first():
        ManualMigrationAttempt.objects.create(
            task=task,
            requires_retry=False,
            note="Retired - deployment_check reports a somalier build that needs replacing",
        )


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0296_one_off_fix_symbolic_del_svlen_sign"),
        ("manual", "0008_systemd_env_files"),
    ]

    operations = [
        migrations.RunPython(_complete_retired_task, migrations.RunPython.noop),
    ]
