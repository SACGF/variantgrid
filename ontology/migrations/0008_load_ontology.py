# Generated by Django 3.1.6 on 2021-03-01 04:36

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_old_ontology_import(apps):
    OntologyImport = apps.get_model("ontology", "OntologyImport")
    return OntologyImport.objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('ontology', '0007_ontologyterm_aliases'),
        ('manual', '0002_deployment')
    ]

    operations = [
        ManualOperation(ManualOperation.operation_other("Download mondo.json and run ontology_import --mondo_json json_file_location",
                                                        note="Mondo parsing has changed significantly from first deploy",
                                                        test=_test_has_old_ontology_import)),
        ManualOperation(ManualOperation.operation_other("Download phenotype_to_genes.txt and run ontology_import --phenotype_to_genes phenotype_to_genes.txt",
                                                        note="Check annotation page for all ontology imports",
                                                        test=_test_has_old_ontology_import)),
    ]