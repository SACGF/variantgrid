from datetime import UTC, datetime

from django.db import migrations, models

from annotation.vep_plugin_backfill import ACTIVE_STATUSES, TRANSCRIPTS_FLAG, backfill_note
from manual.operations.manual_operations import ManualOperation

# The earlier #2426 / #2433 / #2499 backfill steps: this migration's note covers every version they did, in the one
# pass that fills the transcript rows, so an outstanding copy of any would send the admin through VEP twice
SUPERSEDED_TASK_IDS = [
    ManualOperation.task_id_other("Backfill UTRAnnotator, AVI and SpliceVault columns by hand (#2426) - see note"),
    ManualOperation.task_id_other("Backfill AlphaGenome splicing and the #2426 plugin columns by hand (#2433) - see note"),
    ManualOperation.task_id_other("Re-annotate 5' UTR effects reported by Ensembl's UTRAnnotator by hand (#2499) - see note"),
]
TRANSCRIPTS_TASK = "Backfill UTRAnnotator and SpliceVault transcript rows by hand (#2426) - see note"


def _flag_versions_that_ran_transcript_plugins(apps, _schema_editor):
    """ Every NEW/ACTIVE version that has run UTRAnnotator or SpliceVault (a pin set and locked variants) had its
        variants annotated before the switch-on backfilled into VariantAnnotation only """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    AnnotationRangeLock = apps.get_model("annotation", "AnnotationRangeLock")
    annotated = AnnotationRangeLock.objects.filter(version_id=models.OuterRef("pk"))
    ran_a_plugin = models.Q(utr_annotator__isnull=False) | models.Q(splicevault__isnull=False)
    VariantAnnotationVersion.objects.filter(status__in=ACTIVE_STATUSES).filter(ran_a_plugin) \
        .filter(models.Exists(annotated)).update(**{TRANSCRIPTS_FLAG: False})


def _has_versions_to_backfill(apps) -> bool:
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(**{TRANSCRIPTS_FLAG: False}).exists()


def _backfill_note(apps) -> str:
    return backfill_note(apps.get_model("annotation", "VariantAnnotationVersion"))


def _outstanding_required(apps, task_id: str):
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    ManualMigrationRequired = apps.get_model("manual", "ManualMigrationRequired")
    last_success = ManualMigrationAttempt.objects.filter(task_id=task_id, requires_retry=False).order_by("-created").first()
    since = last_success.created if last_success else datetime.min.replace(tzinfo=UTC)
    return ManualMigrationRequired.objects.filter(task_id=task_id, created__gte=since)


def _drop_superseded_backfill_steps(apps, _schema_editor):
    if _outstanding_required(apps, ManualOperation.task_id_other(TRANSCRIPTS_TASK)).exists():
        for task_id in SUPERSEDED_TASK_IDS:
            _outstanding_required(apps, task_id).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0009_utr_effects_backfill"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="backfilled_utr_splicevault_transcripts",
            field=models.BooleanField(default=True),
        ),
        migrations.RunPython(_flag_versions_that_ran_transcript_plugins, migrations.RunPython.noop),
        ManualOperation.operation_other(args=[TRANSCRIPTS_TASK], note=_backfill_note, test=_has_versions_to_backfill),
        migrations.RunPython(_drop_superseded_backfill_steps, migrations.RunPython.noop),
    ]
