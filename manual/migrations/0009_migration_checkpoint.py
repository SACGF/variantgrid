"""
The migration checkpoint (#736): the last migration of the old history. master squashed every app's migrations
into fresh initials without `replaces`, so a deployment reaches it through this branch: applying this migration
means every app is fully migrated and every manual task those migrations registered has been done. Once applied,
manual.checkpoint stops this branch running - switch to master / the release branch.
"""
from django.conf import settings
from django.db import migrations

from manual.upgrader import outstanding_tasks_json

FINAL_MIGRATIONS = [
    ('analysis', '0157_lazy_recount_trio_karyomapping'),
    ('annotation', '0201_annotationrangelock_unique_together_timestamps'),
    ('beacon', '0002_alter_beaconquerycache_exists'),
    ('classification', '0236_import_run_lab_record_ids_already_withdrawn'),
    ('email_manager', '0001_initial'),
    ('eventlog', '0008_api_request_log_notes'),
    ('expression', '0002_remove_cuffdiffrecord_cuff_diff_file_and_more'),
    ('flags', '0011_alter_flag_options_alter_flagcomment_options_and_more'),
    ('genes', '0105_alter_activesamplegenelist_options_and_more'),
    ('manual', '0008_systemd_env_files'),
    ('mme', '0002_inbound_matches_and_disclaimers'),
    ('ontology', '0028_backfill_phenotype_to_genes'),
    ('pathtests', '0006_alter_pathologytest_options_and_more'),
    ('patients', '0026_alter_clinician_options_and_more'),
    ('pedigree', '0004_alter_somalierpedigreerelate_created_and_more'),
    ('review', '0005_discordance_report_topic'),
    ('seqauto', '0055_alter_alignmentfile_options_alter_fastq_options_and_more'),
    ('snpdb', '0298_gate_migrate_common_filter_gnomad_versions'),
    ('sync', '0008_alter_classificationmodificationsyncrecord_options_and_more'),
    ('upload', '0058_alter_fileupload_options_alter_fileupload_created_and_more'),
    ('user_messages', '0002_migrate_legacy_django_messages'),
]
if "sapath" in {app.split(".")[0] for app in settings.INSTALLED_APPS}:
    FINAL_MIGRATIONS.append(('sapath', '0025_merge_tso500_kits_timestamps'))


def _refuse_while_manual_tasks_outstanding(_apps, _schema_editor):
    """ master no longer has the one-off commands these tasks run """
    if settings.UNIT_TEST:
        return  # A fresh test database registers the install-time tasks
    if outstanding := outstanding_tasks_json():
        task_lines = "\n".join(f"  {task['id']}" for task in outstanding)
        raise RuntimeError(f"{len(outstanding)} manual tasks are outstanding - master can't run them. Do them "
                           f"(scripts/upgrade.sh menu, or --auto-manage), mark any that don't apply complete, then "
                           f"migrate again:\n{task_lines}")


class Migration(migrations.Migration):

    dependencies = FINAL_MIGRATIONS

    operations = [
        migrations.RunPython(_refuse_while_manual_tasks_outstanding),
    ]
