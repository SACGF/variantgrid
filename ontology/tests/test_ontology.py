from django.test import TestCase

from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)


class Test(TestCase):
    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        cls.ontology_version = create_test_ontology_version()

    def test_ontology_import_loaders_run(self):
        """ Smoke test - the fixture loaders build the whole test ontology without raising """
        self.assertIsNotNone(self.ontology_version)

    def test_phenotype_to_genes_links_hpo_to_gene_symbols(self):
        """ Five column phenotype_to_genes.txt (2023 format) - HPO -> OMIM -> HGNC, ORPHA rows dropped (#1862) """
        def symbols_for(hpo_id):
            return set(self.ontology_version.gene_symbols_for_terms((hpo_id,)).values_list("symbol", flat=True))

        # CUBN is only ORPHA for this term, BCKDHB has no HGNC term
        self.assertEqual({"BCKDHA", "MEN1"}, symbols_for("HP:0001507"))
        # MEN1 is only ORPHA for this term, AMN has no HGNC term
        self.assertEqual({"BCKDHA", "CUBN"}, symbols_for("HP:0004323"))
