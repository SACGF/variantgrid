from django.test import TestCase

from ontology.models import OntologyService, OntologyTerm
from ontology.ontology_traversal import DbOntologyTraverser, MemoryOntologyTraverser
from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)


class OntologyTraversalParityTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        cls.ov = create_test_ontology_version()

    def test_terms_for_gene_symbol_parity(self):
        db = DbOntologyTraverser(self.ov)
        mem = MemoryOntologyTraverser(self.ov)
        services = [OntologyService.OMIM, OntologyService.HPO, OntologyService.MONDO]
        for hgnc in OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC):
            for service in services:
                db_leaves = sorted(t.pk for t in db.terms_for_gene_symbol(hgnc.name, service, max_depth=1).leafs())
                mem_leaves = sorted(t.pk for t in mem.terms_for_gene_symbol(hgnc.name, service, max_depth=1).leafs())
                self.assertEqual(db_leaves, mem_leaves,
                                 msg=f"Mismatch for {hgnc.name} -> {service}")
