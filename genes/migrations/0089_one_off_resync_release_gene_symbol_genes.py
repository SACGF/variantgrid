from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_alias_derived_release_gene_matches(apps):
    ReleaseGeneSymbolGene = apps.get_model("genes", "ReleaseGeneSymbolGene")
    return ReleaseGeneSymbolGene.objects.filter(match_info__contains="is an alias for").exists()


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0088_one_off_stamp_existing_pfam_domains_imported'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_rematch_release_symbols_to_genes"]),
                        note="Remove chained-alias gene matches (e.g. MT-TS2 -> RP8 -> PDCD2) and resync "
                             "ReleaseGeneSymbolGene to single-hop matching (#1669)",
                        test=_has_alias_derived_release_gene_matches),
    ]
