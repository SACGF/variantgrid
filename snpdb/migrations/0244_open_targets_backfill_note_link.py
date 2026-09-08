from django.db import migrations

from manual.operations.manual_operations import ManualOperation

# snpdb/0242 carried the whole dump/annotate/import recipe as the note, which is more than the
# migrator's task list can usefully show - it lives on the issue now, and the note points at it.
_ARGS = ["Backfill the open_targets_* columns by hand (#1822) - see note"]
_TASK_ID = ManualOperation.task_id_other(_ARGS)

_BACKFILL_NOTE = """Existing GRCh38 columns_version 5 annotation holds only the Open Targets
associations each variant leads, so a variant in someone else's GWAS credible set shows no L2G score.
Optional - the columns fill themselves at the next full re-annotation.
Recipe: https://github.com/SACGF/variantgrid/issues/1822#issuecomment-5519236046"""


def _drop_long_note(apps, _schema_editor):
    apps.get_model("manual", "ManualMigrationRequired").objects.filter(task_id=_TASK_ID).delete()


def _has_open_targets_annotation(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version__gte=5).exclude(open_targets__isnull=True) \
                                           .exclude(open_targets="").exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0243_merge_20260903_1122'),
    ]

    operations = [
        migrations.RunPython(_drop_long_note, migrations.RunPython.noop),
        ManualOperation.operation_other(args=_ARGS, note=_BACKFILL_NOTE, test=_has_open_targets_annotation),
    ]
