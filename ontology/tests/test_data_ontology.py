import logging
import warnings
from pathlib import Path

from django.conf import settings
from django.utils import timezone

from ontology.management.commands import ontology_import
from ontology.models import OntologyVersion, OntologyImport


def create_ontology_test_data():
    # Suppress Pronto warnings, @see https://docs.python.org/3/library/warnings.html#temporarily-suppressing-warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        test_data_dir = Path(settings.BASE_DIR) / "ontology" / "tests" / "test_data"
        biomart_filename = test_data_dir / "biomart_omim.tsv"
        hpo_filename = test_data_dir / "small.owl"

        ontology_import.load_biomart(str(biomart_filename), True)
        # ontology_import.load_mondo(filename, False)
        ontology_import.load_hpo(str(hpo_filename), True)
        # ontology_import.load_hpo_disease(filename, False)


def create_test_ontology_version() -> OntologyVersion:
    kwargs = {}
    now = timezone.now()
    for field, (import_source, filename) in OntologyVersion.ONTOLOGY_IMPORTS.items():
        oi, _ = OntologyImport.objects.get_or_create(import_source=import_source, filename=filename,
                                                     defaults={"processed_date": now})
        kwargs[field] = oi
    ontology_version, _ = OntologyVersion.objects.get_or_create(**kwargs)
    return ontology_version
