from django.db import migrations

from annotation.vep_plugin_backfill import BACKFILL_FLAG, backfill_note
from manual.operations.manual_operations import ManualOperation


def _has_versions_to_backfill(apps) -> bool:
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(**{BACKFILL_FLAG: False}).exists()


def _backfill_note(apps) -> str:
    """ Only registered where the plugins were configured before migrating. Enabling them later is
        `manage.py enable_vep_plugin_backfill`, which prints the steps """
    note = backfill_note(apps.get_model("annotation", "VariantAnnotationVersion"))
    return (note + "\n\nA plugin switched on after this migration: 'manage.py enable_vep_plugin_backfill' stamps "
            "its versions and prints the steps for them")


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0004_utr_annotator_avi_splicevault"),
        ("manual", "0001_initial"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Backfill UTRAnnotator, AVI and SpliceVault columns by hand (#2426) "
                                              "- see note"],
                                        note=_backfill_note, test=_has_versions_to_backfill),
    ]
