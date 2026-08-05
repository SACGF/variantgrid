from django.db import migrations

# One-off "backfill from today": stamp dependency gates (see manual/gates.py) onto historical manage
# tasks. Going forward, new migrations declare these at the call site via
# ManualOperation.operation_manage(..., requires=[...]); this migration covers rows created before
# that field existed. It depends (below) on the migrations that create these tasks so it can't run
# before they exist and match zero rows.
_FVM = "manage*fix_variant_matching"
TASK_REQUIRES = {
    # After the ontology import (phenotype text is matched against ontology terms):
    "manage*match_patient_phenotypes --clear": ["ontology-imported"],
    # After variant annotation is upgraded (operator-confirmed manual gate):
    "manage*calculate_sample_stats": ["variant-annotation-current"],
    "manage*calculate_sample_stats --clear": ["variant-annotation-current"],
    # fix_variant_matching needs current cdot transcripts (the command itself hard-errors otherwise);
    # --extra also wants RefSeq transcript sequences loaded first, else c.hgvs resolution falls back
    # to slow per-transcript NCBI fetches (was the "BEFORE rematching - import_transcript_fasta"
    # reminder, now a gate). The four flags were inserted across successive migrations (classification
    # 0083 --extra, 0089 --revalidate_chgvs, 0091 --sort, 0101 --non-coding); chain them in that order.
    f"{_FVM} --extra": ["cdot-current", "transcript-sequences-loaded"],
    f"{_FVM} --revalidate_chgvs": ["cdot-current", f"after:{_FVM} --extra"],
    f"{_FVM} --sort": ["cdot-current", f"after:{_FVM} --revalidate_chgvs"],
    f"{_FVM} --non-coding": ["cdot-current", f"after:{_FVM} --sort"],
}


def backfill_requires(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    for task_id, requires in TASK_REQUIRES.items():
        ManualMigrationTask.objects.filter(pk=task_id).update(requires=requires)


def clear_requires(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    for task_id in TASK_REQUIRES:
        ManualMigrationTask.objects.filter(pk=task_id).update(requires=list())


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0004_complete_obsolete_manual_tasks'),
        # Must run after the migrations that create these tasks, else the update matches nothing.
        ('classification', '0101_one_off_fix_variant_matching_non_coding'),
        ('annotation', '0037_one_off_move_pheno_match_to_hgnc'),
        ('patients', '0004_one_off_fix_patient_imports_and_modifications'),
    ]

    operations = [
        migrations.RunPython(backfill_requires, clear_requires),
    ]
