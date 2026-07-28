from django.db import migrations

# Manage-category commands that were deleted along with the feature they backfilled (audited: none
# were wrongly deleted). A DB that ran the old (pre-neutralised) migrations still holds an outstanding
# task row pointing at these now-missing commands; mark them complete so they drop out of the upgrade
# flow instead of failing when run. See the neutralised source migrations (e.g. snpdb/0100, genes/0027).
OBSOLETE_COMMANDS = {
    "one_off_calc_variant_end",                      # superseded by one_off_fix_variant_end (#414)
    "fix_panel_app_gene_list_permissions",           # feature removed (#690)
    "fix_historical_max_af",                         # max-af feature removed
    "remove_redundant_flags",                        # moved to a data migration (#1577)
    "fix_retrieve_transcript_version_sequence_info",
    "vep_download_fasta",
}


def complete_obsolete_tasks(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    for task in ManualMigrationTask.objects.all():
        category, _, line = task.id.partition("*")
        if category != "manage" or not line:
            continue
        command = line.split()[0]
        if command in OBSOLETE_COMMANDS:
            # A success attempt (requires_retry=False) dated now clears any earlier requirement, so the
            # task is no longer outstanding. Harmless if it was already complete.
            ManualMigrationAttempt.objects.create(
                task=task,
                requires_retry=False,
                note="Obsolete - command removed; auto-completed by manual/0004",
            )


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0003_manualgatesatisfied_manualmigrationtask_requires'),
    ]

    operations = [
        migrations.RunPython(complete_obsolete_tasks, migrations.RunPython.noop),
    ]
