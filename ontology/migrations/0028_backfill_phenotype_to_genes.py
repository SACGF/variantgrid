from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _latest_phenotype_to_genes_import_has_no_relations(apps):
    """ The 2023 phenotype_to_genes.txt format imported as zero relations, so no HPO term reached a gene
        @see https://github.com/SACGF/variantgrid/issues/1862 """
    OntologyVersion = apps.get_model("ontology", "OntologyVersion")
    OntologyTermRelation = apps.get_model("ontology", "OntologyTermRelation")
    if ontology_version := OntologyVersion.objects.order_by("pk").last():
        return not OntologyTermRelation.objects.filter(
            from_import_id=ontology_version.hp_phenotype_to_genes_import_id).exists()
    return False


class Migration(migrations.Migration):

    dependencies = [
        ('ontology', '0027_alter_ontologyterm_ontology_service'),
        ('manual', '0003_manualgatesatisfied_manualmigrationtask_requires'),
    ]

    operations = [
        ManualOperation.operation_manage(["backfill_phenotype_to_genes"],
                                         note="Latest OntologyVersion has no HPO -> OMIM -> gene relations; downloads "
                                              "phenotype_to_genes.txt and loads it into the existing import (#1862)",
                                         test=_latest_phenotype_to_genes_import_has_no_relations),
    ]
