from django.db import migrations

from manual.operations.manual_operations import ManualOperation
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME


def _has_gene_level_allele_infos(apps):
    """ Gene-level records validated before #1835 carry a fatal transcript_type_not_supported tag, so they
        are excluded from every export until they are re-validated """
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    return ImportedAlleleInfo.objects.filter(
        variant_coordinate__startswith=f"{GENE_LEVEL_CONTIG_NAME}:").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("classification", "0183_ekey_splicing_variant_options"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["revalidate_gene_level_allele_infos"]),
                        note="Re-validate gene-level classifications so they are no longer tagged as an "
                             "unsupported transcript type and excluded from exports (#1835)",
                        test=_has_gene_level_allele_infos),
    ]
