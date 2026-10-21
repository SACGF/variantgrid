from django.db import migrations

from manual.operations.manual_operations import ManualOperation

# 0192 imports cdot's 'fake-rna-<gene>' MT transcripts; this adds them to the releases imported before cdot had them
IMPORT_CDOT_TASK_ID = ManualOperation.task_id_manage(["import_cdot_latest", "--annotation-consortium", "RefSeq",
                                                      "--force"])


def _has_refseq_gene_annotation_release(apps):
    GeneAnnotationRelease = apps.get_model("genes", "GeneAnnotationRelease")
    return GeneAnnotationRelease.objects.filter(annotation_consortium="R").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0194_phenotype_models_to_patients"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_gene_annotation_release_mitochondrial"]),
                        note="Add the mitochondrial genes to RefSeq gene annotation releases and their "
                             "gene annotation (#2139)",
                        requires=[f"after:{IMPORT_CDOT_TASK_ID}"],
                        test=_has_refseq_gene_annotation_release),
    ]
