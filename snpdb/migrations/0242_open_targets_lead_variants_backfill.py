from django.db import migrations

from manual.operations.manual_operations import ManualOperation

# snpdb/0240's recipe re-ran VEP with the pipeline's flags, which took the OpenTargets plugin's
# lead_variants_only default of 1 - so it only ever saw the associations the variant leads. Now the
# pipeline passes lead_variants_only=0 the record sets differ, and open_targets_gwas_l2g_scores can
# only be imported alongside the sibling arrays it is zipped with.
_OLD_ARGS = ["Backfill open_targets_gwas_l2g_scores by hand (#1822) - see note"]
_ARGS = ["Backfill the open_targets_* columns by hand (#1822) - see note"]

_OPEN_TARGETS_COLUMNS = {
    "open_targets_gwas_l2g_score": "OpenTargets_gwasLocusToGeneScore",
    "open_targets_gwas_l2g_scores": "OpenTargets_gwasLocusToGeneScore",
    "open_targets_gwas_gene_id": "OpenTargets_gwasGeneId",
    "open_targets_gwas_diseases": "OpenTargets_gwasDiseases",
    "open_targets_study_type": "OpenTargets_studyType",
    "open_targets_study_id": "OpenTargets_studyId",
    "open_targets_is_lead": "OpenTargets_isLead",
    "open_targets_variant_id": "OpenTargets_variantId",
    "open_targets_qtl_gene_id": "OpenTargets_qtlGeneId",
    "open_targets_qtl_biosample": "OpenTargets_qtlBiosampleName",
}
_COLUMNS = ",".join(_OPEN_TARGETS_COLUMNS)
_CSQ = " \\\n       ".join(f"--csq {column}={source_field}"
                          for column, source_field in _OPEN_TARGETS_COLUMNS.items())
_FIELDS = ",".join(["Feature", "PICK"] + sorted(set(_OPEN_TARGETS_COLUMNS.values())))

_BACKFILL_NOTE = f"""Only GRCh38 annotation made at columns_version 5 has Open Targets data to backfill.
The columns fill themselves at the next full re-annotation, so this is optional - skip it if one is
coming up. Existing annotation holds only the associations each variant leads, so a variant sitting in
someone else's GWAS credible set shows no L2G score; this re-reads every association it belongs to,
open_targets_is_lead saying which of them the variant leads.
The open_targets_* arrays are parallel and displayed zipped together, so they are all rewritten from
the one VEP run. --not-null keeps the dump to variants Open Targets already matched - a variant whose
only associations are non-lead has to wait for the re-annotation.

1. python3 manage.py annotation_backfill_columns --dump --genome-build GRCh38 \\
       --columns {_COLUMNS} \\
       --not-null open_targets_variant_id --output ot.vcf.gz
2. Annotate ot.vcf.gz with VEP, using the pipeline's flags for the build. The plugin matches on
   position and alleles, so it is the only one that needs loading, and
   --plugin OpenTargets,file=<open_targets_NN.NN_vep.tsv.bgz>,cols=all,lead_variants_only=0
   --fields {_FIELDS}
   keeps the output to what the import reads
3. python3 manage.py annotation_backfill_columns --import --genome-build GRCh38 \\
       --columns {_COLUMNS} \\
       {_CSQ} \\
       --clear-missing --vcf annotated.vcf.gz"""


def _drop_old_note(apps, _schema_editor):
    """ The 0239/0240 recipe now writes arrays out of step with their siblings - retire it whether or
        not a deployment has run it, so the only request standing is the one below """
    task_id = ManualOperation.task_id_other(_OLD_ARGS)
    apps.get_model("manual", "ManualMigrationRequired").objects.filter(task_id=task_id).delete()


def _has_open_targets_annotation(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version__gte=5).exclude(open_targets__isnull=True) \
                                           .exclude(open_targets="").exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0241_open_targets_is_lead_column'),
    ]

    operations = [
        migrations.RunPython(_drop_old_note, migrations.RunPython.noop),
        ManualOperation.operation_other(args=_ARGS, note=_BACKFILL_NOTE, test=_has_open_targets_annotation),
    ]
