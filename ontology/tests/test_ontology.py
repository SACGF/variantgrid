from django.test import TestCase

from ontology.tests.test_data_ontology import create_ontology_test_data


class Test(TestCase):
    def test_ontology_import_loaders_run(self):
        """ Smoke test - the fixture loaders build the whole test ontology without raising """
        create_ontology_test_data()
