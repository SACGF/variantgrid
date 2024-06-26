# Generated by Django 3.1 on 2021-01-12 04:55

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_old_ontology_import(apps):
    OntologyImport = apps.get_model("ontology", "OntologyImport")
    return OntologyImport.objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0002_deployment'),
        ('ontology', '0001_initial'),
    ]

    operations = [
        ManualOperation.operation_other("Download mondo.json and run ontology_import --mondo_json json_file_location",
                                        test=_test_has_old_ontology_import)
    ]
