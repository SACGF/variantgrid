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

# Obsolete steps whose command still EXISTS (so command-name matching would wrongly clear sibling
# variants) - matched by exact task id instead. These are historical in-place patches to EXISTING
# annotation. We are doing a full fresh annotation reload (old annotation kept only as legacy), so
# patching legacy annotation is moot - retire them here (matched by exact id so live commands and
# other flag variants are untouched). fix_annotation_link_transcripts is deliberately NOT listed:
# transcript links are heavily used to reach gene annotation and are quick, so that step still runs.
OBSOLETE_TASK_IDS = {
    # gene_annotation (gene-level) partial updates - all skipped for the fresh reload
    "manage*gene_annotation --missing",
    "manage*gene_annotation --add-new-to-existing",
    "manage*gene_annotation --add-dbnsfp-gene",
    "manage*gene_annotation --fix-bad-gene-annotation",
    "manage*gene_annotation --add-missing-omim",
    # variant-annotation column/field backfills onto existing annotation - skipped for the fresh reload
    "manage*fix_columns_version2_damage_counts",
    "manage*fix_columns_version3_gnomad_hemi_count",
    "manage*fix_columns_version4_damage_counts",
    "manage*fix_variant_annotation_add_hgvs_g",
    "manage*fix_annotation_sv_c_hgvs",
    "manage*fix_annotation_sv_overlaps",
    "manage*fix_historical_spliceai_max_ds",
}


def complete_obsolete_tasks(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    for task in ManualMigrationTask.objects.all():
        category, _, line = task.id.partition("*")
        if category != "manage" or not line:
            continue
        if line.split()[0] in OBSOLETE_COMMANDS or task.id in OBSOLETE_TASK_IDS:
            # A success attempt (requires_retry=False) dated now clears any earlier requirement, so the
            # task is no longer outstanding. Harmless if it was already complete.
            ManualMigrationAttempt.objects.create(
                task=task,
                requires_retry=False,
                note="Obsolete - auto-completed by manual/0004",
            )


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0003_manualgatesatisfied_manualmigrationtask_requires'),
        # The OBSOLETE_TASK_IDS skips are all created by annotation migrations (0022-0146). Depend on
        # the latest so this cleanup runs after every creator - otherwise a still-pending creator could
        # register the task after we clear it. (The OBSOLETE_COMMANDS orphans need no such dep: their
        # creators are already applied or neutralised.)
        ('annotation', '0146_one_off_cols_v4_pathogenic_counts'),
    ]

    operations = [
        migrations.RunPython(complete_obsolete_tasks, migrations.RunPython.noop),
    ]
