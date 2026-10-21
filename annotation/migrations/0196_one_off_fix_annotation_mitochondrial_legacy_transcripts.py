from django.db import migrations

from manual.operations.manual_operations import ManualOperation

IMPORT_CDOT_TASK_ID = ManualOperation.task_id_manage(["import_cdot_latest", "--annotation-consortium", "RefSeq",
                                                      "--force"])


def _has_refseq_variant_annotation(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(annotation_consortium="R").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0195_one_off_fix_gene_annotation_release_mitochondrial"),
    ]

    # Same task as 0192, requested again: on a long-lived database VEP's 'ND4' matched a gene-named pre-cdot
    # TranscriptVersion rather than nothing, so the first run found no unlinked MT rows (#2139)
    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_annotation_mitochondrial_transcripts"]),
                        note="Relink RefSeq mitochondrial annotation from legacy gene-named transcripts to cdot's "
                             "fake transcripts and delete the legacy rows (#2139)",
                        requires=[f"after:{IMPORT_CDOT_TASK_ID}"],
                        test=_has_refseq_variant_annotation),
    ]
