from django.db import migrations

from manual.operations.manual_operations import ManualOperation

# snpdb/0239's note dumped every annotated variant, which on a full deployment is millions of records
# to put through VEP. The plugin writes all of its columns together, so only variants already holding
# Open Targets annotation can have a score - --not-null cuts the dump down to those.
_ARGS = ["Backfill open_targets_gwas_l2g_scores by hand (#1822) - see note"]
_TASK_ID = ManualOperation.task_id_other(_ARGS)

_BACKFILL_NOTE = """Only GRCh38 annotation made at columns_version 5 has Open Targets data to backfill.
The column fills itself at the next full re-annotation, so this is optional - skip it if one is coming up.

1. python3 manage.py annotation_backfill_columns --dump --genome-build GRCh38 \\
       --columns open_targets_gwas_l2g_scores --only-missing \\
       --not-null open_targets_variant_id --output ot_l2g.vcf.gz
2. Annotate ot_l2g.vcf.gz with VEP and the OpenTargets plugin, using the pipeline's flags for the
   build. The plugin matches on position and alleles, so it is the only one that needs loading,
   and --fields Feature,PICK,OpenTargets_gwasLocusToGeneScore keeps the output small
3. python3 manage.py annotation_backfill_columns --import --genome-build GRCh38 \\
       --columns open_targets_gwas_l2g_scores \\
       --csq open_targets_gwas_l2g_scores=OpenTargets_gwasLocusToGeneScore --vcf annotated.vcf.gz"""


def _already_backfilled(apps) -> bool:
    ManualMigrationAttempt = apps.get_model("manual", "ManualMigrationAttempt")
    return ManualMigrationAttempt.objects.filter(task_id=_TASK_ID, requires_retry=False).exists()


def _still_outstanding(apps) -> bool:
    """ Same test snpdb/0239 used, less any deployment that has since run its recipe - re-registering
        there would ask for a backfill that's already done """
    if _already_backfilled(apps):
        return False
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version__gte=5).exclude(open_targets__isnull=True) \
                                           .exclude(open_targets="").exists()


def _supersede_old_note(apps, _schema_editor):
    """ Drop the request carrying 0239's note, so an operator yet to run it sees only the cheaper recipe """
    if not _already_backfilled(apps):
        apps.get_model("manual", "ManualMigrationRequired").objects.filter(task_id=_TASK_ID).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0239_open_targets_l2g_scores_manual_backfill'),
    ]

    operations = [
        migrations.RunPython(_supersede_old_note, migrations.RunPython.noop),
        ManualOperation.operation_other(args=_ARGS, note=_BACKFILL_NOTE, test=_still_outstanding),
    ]
