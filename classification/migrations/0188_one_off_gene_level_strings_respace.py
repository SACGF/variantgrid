from django.db import migrations

from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from manual.operations.manual_operations import ManualOperation

SPLICE_ALT_PREFIX = f"<{GeneLevelSymbolicAlt.SPLICE}:"


def _has_gene_level_strings_to_respace(apps):
    """ A copy number target imported before gene-level values kept their spaces reads
        'FGF4amplification'; a splice variant's annotation predates the 'splice' suffix """
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    Variant = apps.get_model("snpdb", "Variant")
    stripped_copy_number = ImportedAlleleInfo.objects.filter(
        imported_c_hgvs__iregex=r"^[A-Za-z0-9.\-]+(amplification|amp|gain|loss|deletion|del)$").exists()
    return stripped_copy_number or Variant.objects.filter(alt__seq__startswith=SPLICE_ALT_PREFIX).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("classification", "0187_copy_number_ekey_label"),
        ("annotation", "0186_alter_clinvar_version_alter_geneannotation_version_and_more"),
        ("genes", "0096_alter_genecoverage_gene_coverage_collection_and_more"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["gene_level_strings_respace"]),
                        note="Put the space back in copy number classification targets ('FGF4amplification' -> "
                             "'FGF4 amplification', so c.HGVS reads as g.HGVS does) and add the 'splice' "
                             "suffix to splice annotation and evidence (#1875). Run "
                             "'manage.py gene_level_strings_respace --dry-run' first",
                        test=_has_gene_level_strings_to_respace),
    ]
