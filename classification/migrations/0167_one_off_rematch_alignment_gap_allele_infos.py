from django.db import migrations

from manual.operations.manual_operations import ManualOperation

RETIRED_TASK_ID = ManualOperation.task_id_manage(["fix_legacy_classification_alignment_gap_hgvs_matching"])
_FVM = "manage*fix_variant_matching"


def _test_has_allele_infos(apps):
    """ The command works out for itself which records (if any) still hold matching made without alignment
        gap data - the converter version stamped on an ImportedAlleleInfo can't tell us, as it's re-stamped
        by coordinate/c.hgvs recalculation without the matched variant being re-derived """
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    return ImportedAlleleInfo.objects.exists()


def complete_retired_task(apps, schema_editor):
    """ fix_legacy_classification_alignment_gap_hgvs_matching has been deleted (it queried the pre-cdot
        transcript JSON schema and rematched via the pre-ImportedAlleleInfo API), so clear any outstanding
        row for it - the manage task registered below replaces it """
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    if task := ManualMigrationTask.objects.filter(pk=RETIRED_TASK_ID).first():
        ManualMigrationAttempt.objects.create(
            task=task,
            requires_retry=False,
            note="Retired - replaced by fix_legacy_imported_allele_info_alignment_gap",
        )


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0166_reminder_legacy_somatic_forwardport'),
        ('manual', '0005_backfill_manual_task_requires'),
    ]

    operations = [
        migrations.RunPython(complete_retired_task, migrations.RunPython.noop),
        # Runs against ImportedAlleleInfo, so it needs fix_variant_matching --extra to have linked
        # classifications to them first, and the same data gates that step uses.
        ManualOperation.operation_manage(
            ["fix_legacy_imported_allele_info_alignment_gap"],
            test=_test_has_allele_infos,
            requires=["cdot-current", "transcript-sequences-loaded", f"after:{_FVM} --extra"]),
    ]
