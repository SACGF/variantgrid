from django.core.management import get_commands
from django.db import migrations

# Any manage task whose command no longer exists can't be run, so it's completed generically below
# (a DB that ran the old, pre-neutralised migrations still holds such orphan rows). Audited examples,
# none wrongly deleted: one_off_calc_variant_end (#414), fix_panel_app_gene_list_permissions (#690),
# fix_historical_max_af, remove_redundant_flags, fix_retrieve_transcript_version_sequence_info,
# vep_download_fasta, classification_populate_sort_orders (column dropped again by classification/0155,
# sorting moved to Classification.summary - the command could only ever fail, so it was deleted).

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


# Obsolete 'other' (free-text) reminders. These were prerequisite/ordering notes now expressed as
# gates (import_transcript_fasta -> transcript-sequences-loaded; import gene annotation -> cdot-current)
# or superseded (deployment_check now verifies cdot + gnomAD files; one_off_fix_variant_end has a real
# manage task; dbNSFP skipped for the fresh reload; VEP110 advisory is moot). Retired by exact task id
# (single-arg operation_other, so the id is other*"<text>"). #1 has no current source - a DB orphan.
def _other_id(text: str) -> str:
    return f'other*"{text}"'


OBSOLETE_OTHER_TASK_IDS = {
    _other_id("VariantZygosityCount: Re-calculate and change config: "
              "https://github.com/SACGF/variantgrid/issues/1571#issuecomment-688651886"),
    _other_id("*** BEFORE rematching - import_transcript_fasta - see annotation page"),
    _other_id("** Before fix_legacy_classification_alignment_gap_hgvs_matching ** - Import gene "
              "annotation - see https://github.com/SACGF/variantgrid/issues/494#issuecomment-943977004"),
    _other_id("Download cdot data and run import_gene_annotation. "
              "See https://github.com/SACGF/variantgrid/issues/566#issuecomment-1097543081"),
    _other_id("Get gnomad AF>5 vcf.gz for upload pre-processing "
              "see https://github.com/SACGF/variantgrid/issues/521#issuecomment-1058830717"),
    _other_id("Import dbNSFP gene annotation (see annotation page)"),
    _other_id("Run 'python3 manage.py one_off_fix_variant_end' - this can be done outside of "
              "migrations (use screen)"),
    _other_id("You may want to update your annotation to VEP110 and latest annotation data "
              "See https://github.com/SACGF/variantgrid/wiki/Annotation-Column-Versions"),
}


def _is_obsolete(task_id: str, known_commands: set) -> bool:
    category, _, line = task_id.partition("*")
    if not line:
        return False
    if category == "manage":
        # Command deleted -> can't run, so obsolete; plus vetted skips whose command still exists.
        return line.split()[0] not in known_commands or task_id in OBSOLETE_TASK_IDS
    return task_id in OBSOLETE_OTHER_TASK_IDS


def complete_obsolete_tasks(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    known_commands = set(get_commands())
    for task in ManualMigrationTask.objects.all():
        if _is_obsolete(task.id, known_commands):
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
        # Run after every creator of a task we clear, otherwise a still-pending creator could register
        # the task after we clear it. OBSOLETE_TASK_IDS (annotation) + the obsolete 'other' reminders
        # span the annotation/genes/snpdb apps; depend on the latest creator in each. (The
        # OBSOLETE_COMMANDS orphans need no such dep: their creators are already applied or neutralised.)
        ('annotation', '0146_one_off_cols_v4_pathogenic_counts'),
        ('genes', '0054_one_off_import_gene_annotation_manual'),
        ('snpdb', '0118_reminder_one_off_fix_variant_end'),
    ]

    operations = [
        migrations.RunPython(complete_obsolete_tasks, migrations.RunPython.noop),
    ]
