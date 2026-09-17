from django.db import migrations

from classification.models import ImportedAlleleInfoStatus
from manual.operations.manual_operations import ManualOperation
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME


def _has_failed_gene_level_allele_infos(apps):
    """ A splice target resolved before #1835 has the label in its seeded case, which never matched the
        upper-case Sequence the pipeline inserts, so the record was left Failed with no variant """
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    return ImportedAlleleInfo.objects.filter(
        status=ImportedAlleleInfoStatus.FAILED,
        variant_coordinate__startswith=f"{GENE_LEVEL_CONTIG_NAME}:").exists()


class Migration(migrations.Migration):
    dependencies = [
        # After the canonicalise, so the re-match reaches the relabelled variants rather than minting
        # a second set beside them
        ("classification", "0185_one_off_canonicalise_splice_labels"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["classification_rematch_stuck", "--status", "F",
                                                               "--gene-level", "--older-than-hours", "0"]),
                        note="Re-match gene-level classifications left Failed by a splice alt whose label case "
                             "did not match the stored Sequence (#1835)",
                        test=_has_failed_gene_level_allele_infos),
    ]
