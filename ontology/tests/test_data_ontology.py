from pathlib import Path

from django.conf import settings

from ontology.management.commands import ontology_import


def create_ontology_test_data():
    test_data_dir = Path(settings.BASE_DIR) / "ontology" / "tests" / "test_data"
    biomart_filename = test_data_dir / "biomart_omim.tsv"
    hpo_filename = test_data_dir / "small.owl"

    ontology_import.load_biomart(str(biomart_filename), True)
    # ontology_import.load_mondo(filename, False)
    ontology_import.load_hpo(str(hpo_filename), True)
    # ontology_import.load_hpo_disease(filename, False)
