from datetime import UTC, datetime

from django.db import migrations

from manual.operations.manual_operations import ManualOperation

PLUGIN_BACKFILL_TASK_ID = ManualOperation.task_id_other(
    "Backfill UTRAnnotator, AVI and SpliceVault columns by hand (#2426) - see note")
ALPHAGENOME_BACKFILL_TASK_ID = ManualOperation.task_id_other(
    "Backfill AlphaGenome splicing and the #2426 plugin columns by hand (#2433) - see note")


def _outstanding_required(apps, task_id: str):
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    ManualMigrationRequired = apps.get_model("manual", "ManualMigrationRequired")
    last_success = ManualMigrationAttempt.objects.filter(task_id=task_id, requires_retry=False).order_by("-created").first()
    since = last_success.created if last_success else datetime.min.replace(tzinfo=UTC)
    return ManualMigrationRequired.objects.filter(task_id=task_id, created__gte=since)


def _drop_superseded_plugin_backfill(apps, _schema_editor):
    """ 0007's note has the steps for every flagged version and all four sources, so the #2426 step registered
        alongside it (whose commands predate the AlphaGenome pin) is redundant - and following it alone would mark
        versions backfilled without AlphaGenome """
    if _outstanding_required(apps, ALPHAGENOME_BACKFILL_TASK_ID).exists():
        _outstanding_required(apps, PLUGIN_BACKFILL_TASK_ID).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0007_alphagenome_splicing"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        migrations.RunPython(_drop_superseded_plugin_backfill, migrations.RunPython.noop),
    ]
