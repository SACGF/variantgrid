from django.db import migrations

from manual.operations.manual_operations import ManualOperation

RETIRED_TASK_ID = ManualOperation.task_id_manage(["one_off_recalculate_trio_karyomapping"])


def _delete_karyomapping_counts(apps, schema_editor):
    """ Counts from before reference calls were binned are near zero for joint-called trios. Deleted rather than
        recounted here - the trio page queues a recount for a trio with none (#1754) """
    ContigKaryomappingCounts = apps.get_model("analysis", "ContigKaryomappingCounts")
    GenomeKaryomappingCounts = apps.get_model("analysis", "GenomeKaryomappingCounts")
    ContigKaryomappingCounts.objects.all().delete()
    GenomeKaryomappingCounts.objects.all().delete()


def _complete_retired_task(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    if task := ManualMigrationTask.objects.filter(pk=RETIRED_TASK_ID).first():
        ManualMigrationAttempt.objects.create(
            task=task,
            requires_retry=False,
            note="Retired - trio karyomapping is recounted when the trio page is viewed",
        )


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0156_alter_allelefrequencynode_options_and_more"),
        ("manual", "0008_systemd_env_files"),
    ]

    operations = [
        migrations.RunPython(_delete_karyomapping_counts, migrations.RunPython.noop),
        migrations.RunPython(_complete_retired_task, migrations.RunPython.noop),
    ]
