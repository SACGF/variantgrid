from django.db import migrations

from manual.operations.manual_operations import ManualOperation

# snpdb/0236 registered this as a 'manage' task, so the migrator offers to run it - but
# annotation_backfill_columns can't run unattended: it needs a VEP run over the dumped VCF in
# between its --dump and --import halves, and the score is a plugin (CSQ) field rather than a top
# level INFO one, so the source has to be named on the command line.
_OLD_TASK_ID = ManualOperation.task_id_manage(["annotation_backfill_columns"])

_BACKFILL_NOTE = """Only GRCh38 annotation made at columns_version 5 has Open Targets data to backfill.
The column fills itself at the next full re-annotation, so this is optional - skip it if one is coming up.

1. python3 manage.py annotation_backfill_columns --dump --genome-build GRCh38 \\
       --columns open_targets_gwas_l2g_scores --only-missing --output ot_l2g.vcf.gz
2. Annotate ot_l2g.vcf.gz with VEP and the OpenTargets plugin, as the annotation pipeline does
3. python3 manage.py annotation_backfill_columns --import --genome-build GRCh38 \\
       --columns open_targets_gwas_l2g_scores \\
       --csq open_targets_gwas_l2g_scores=OpenTargets_gwasLocusToGeneScore --vcf annotated.vcf.gz"""


def _retire_manage_task(apps, _schema_editor):
    apps.get_model("manual", "ManualMigrationTask").objects.filter(pk=_OLD_TASK_ID).delete()


def _has_open_targets_annotation(apps):
    """ Same test snpdb/0236 used - only deployments holding Open Targets annotation have anything
        to backfill into the new per-record score column """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version__gte=5).exclude(open_targets__isnull=True) \
                                           .exclude(open_targets="").exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0238_composite_columns'),
        ('manual', '0005_backfill_manual_task_requires'),
    ]

    operations = [
        migrations.RunPython(_retire_manage_task, migrations.RunPython.noop),
        ManualOperation.operation_other(
            args=["Backfill open_targets_gwas_l2g_scores by hand (#1822) - see note"],
            note=_BACKFILL_NOTE,
            test=_has_open_targets_annotation),
    ]
