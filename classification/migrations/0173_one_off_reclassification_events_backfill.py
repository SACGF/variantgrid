from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_classifications_without_timelines(apps):
    """ The reclassification analytics page reads a materialised timeline that has to be built from
        existing published modification history - @see https://github.com/SACGF/variantgrid/issues/1523 """
    Classification = apps.get_model("classification", "Classification")
    ReclassificationEvent = apps.get_model("classification", "ReclassificationEvent")
    return Classification.objects.exists() and not ReclassificationEvent.objects.exists()


class Migration(migrations.Migration):
    dependencies = [
        ("classification", "0172_reclassificationevent"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["reclassification_events_backfill"]),
                        note="Build germline clinical significance timelines from published modification "
                             "history, so the reclassification analytics page has data (#1523)",
                        test=_has_classifications_without_timelines),
    ]
