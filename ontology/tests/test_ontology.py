from django.test import TestCase

from ontology.models import OntologyService, OntologySnake, OntologyTerm
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

    def test_unknown_gene_symbol_has_no_hgnc_term(self):
        """ #999 - callers treat "no HGNC" as "no relationships" rather than an error """
        self.assertIsNone(OntologyTerm.get_gene_symbol_or_none("NOTAGENE123"))

    def test_unknown_gene_symbol_has_no_relationships(self):
        """ #999 - every OntologySnake gene entry point returns empty rather than raising.
            terms_for_gene_symbol matches MemoryOntologyTraverser, which already did """
        self.assertEqual([], list(OntologySnake.terms_for_gene_symbol("NOTAGENE123", OntologyService.MONDO)))
        self.assertEqual(set(), OntologySnake.mondo_terms_for_gene_symbol("NOTAGENE123"))
        self.assertEqual([], OntologySnake.direct_relationships_for_gene_symbol("NOTAGENE123"))
