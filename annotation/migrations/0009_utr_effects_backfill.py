from datetime import UTC, datetime

from django.db import migrations, models

from annotation.vep_plugin_backfill import ACTIVE_STATUSES, UTR_EFFECTS_FLAG, backfill_note
from manual.operations.manual_operations import ManualOperation

# The #2426 / #2433 backfill steps registered by 0005 and 0007: this migration's note covers every version they
# did, in one pass with the UTR effects, so an outstanding copy of either would send the admin through VEP twice
SUPERSEDED_TASK_IDS = [
    ManualOperation.task_id_other("Backfill UTRAnnotator, AVI and SpliceVault columns by hand (#2426) - see note"),
    ManualOperation.task_id_other("Backfill AlphaGenome splicing and the #2426 plugin columns by hand (#2433) - see note"),
]
UTR_EFFECTS_TASK = "Re-annotate 5' UTR effects reported by Ensembl's UTRAnnotator by hand (#2499) - see note"


def _flag_versions_annotated_by_old_plugin(apps, _schema_editor):
    """ Every NEW/ACTIVE version that has run UTRAnnotator (the pin is set and it has locked variants) did so with
        Ensembl's plugin, whichever of a variant's effects it chose """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    AnnotationRangeLock = apps.get_model("annotation", "AnnotationRangeLock")
    annotated = AnnotationRangeLock.objects.filter(version_id=models.OuterRef("pk"))
    VariantAnnotationVersion.objects.filter(status__in=ACTIVE_STATUSES, utr_annotator__isnull=False) \
        .filter(models.Exists(annotated)).update(**{UTR_EFFECTS_FLAG: False})


def _has_versions_to_backfill(apps) -> bool:
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(**{UTR_EFFECTS_FLAG: False}).exists()


def _backfill_note(apps) -> str:
    return backfill_note(apps.get_model("annotation", "VariantAnnotationVersion"),
                         apps.get_model("annotation", "AnnotationRangeLock"))


def _outstanding_required(apps, task_id: str):
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    ManualMigrationRequired = apps.get_model("manual", "ManualMigrationRequired")
    last_success = ManualMigrationAttempt.objects.filter(task_id=task_id, requires_retry=False).order_by("-created").first()
    since = last_success.created if last_success else datetime.min.replace(tzinfo=UTC)
    return ManualMigrationRequired.objects.filter(task_id=task_id, created__gte=since)


def _drop_superseded_backfill_steps(apps, _schema_editor):
    if _outstanding_required(apps, ManualOperation.task_id_other(UTR_EFFECTS_TASK)).exists():
        for task_id in SUPERSEDED_TASK_IDS:
            _outstanding_required(apps, task_id).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0008_drop_superseded_plugin_backfill_step"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="backfilled_utr_effects",
            field=models.BooleanField(default=True),
        ),
        migrations.RunPython(_flag_versions_annotated_by_old_plugin, migrations.RunPython.noop),
        ManualOperation.operation_other(args=[UTR_EFFECTS_TASK], note=_backfill_note, test=_has_versions_to_backfill),
        migrations.RunPython(_drop_superseded_backfill_steps, migrations.RunPython.noop),
    ]
