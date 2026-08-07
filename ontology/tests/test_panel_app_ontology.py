from django.test import TestCase
from django.utils import timezone

from genes.models import PanelAppServer
from genes.panel_app import get_hgnc_pk_from_api_record
from ontology.models import (
    OntologyImport,
    OntologyRelation,
    OntologyService,
    OntologyTerm,
    OntologyTermRelation,
    OntologyTermStatus,
)
from ontology.panel_app_ontology import _update_gene_relations, hgnc_ontology_terms_by_pk


class TestGetHgncPkFromApiRecord(TestCase):
    def test_parses_hgnc_id(self):
        self.assertEqual(10593, get_hgnc_pk_from_api_record({"gene_data": {"hgnc_id": "HGNC:10593"}}))

    def test_parses_bare_number(self):
        self.assertEqual(10593, get_hgnc_pk_from_api_record({"gene_data": {"hgnc_id": "10593"}}))

    def test_missing_hgnc_id(self):
        self.assertIsNone(get_hgnc_pk_from_api_record({"gene_data": {"gene_symbol": "BRCA1"}}))

    def test_missing_gene_data(self):
        self.assertIsNone(get_hgnc_pk_from_api_record({}))

    def test_null_gene_data(self):
        self.assertIsNone(get_hgnc_pk_from_api_record({"gene_data": None}))

    def test_unparsable_hgnc_id(self):
        self.assertIsNone(get_hgnc_pk_from_api_record({"gene_data": {"hgnc_id": "HGNC:junk"}}))


class TestPanelAppHgncResolution(TestCase):
    """ PanelApp records carry an HGNC ID, so gene relations resolve from that rather than matching on
        their gene symbol - which comes from a dated HGNC snapshot and can be a symbol we can't match. """

    def setUp(self):
        # PanelAppServer rows are created by a data migration
        self.server = PanelAppServer.australia_instance()
        self.o_import = OntologyImport.objects.create(import_source=OntologyService.HGNC,
                                                      filename="test", context="test", hash="N/A",
                                                      processor_version=1, completed=True,
                                                      processed_date=timezone.now())
        self.hgnc_term = OntologyTerm.objects.create(id="HGNC:32331", ontology_service=OntologyService.HGNC,
                                                     index=32331, name="CFAP276",
                                                     status=OntologyTermStatus.CONDITION,
                                                     from_import=self.o_import)
        self.omim_term = OntologyTerm.objects.create(id="OMIM:617000", ontology_service=OntologyService.OMIM,
                                                     index=617000, name="Some condition",
                                                     status=OntologyTermStatus.CONDITION,
                                                     from_import=self.o_import)

    def _panel_record(self):
        return [{
            "gene_data": {"gene_symbol": "UNMATCHABLE_SYMBOL", "hgnc_id": "HGNC:32331"},
            "confidence_level": "3",
            "evidence": ["Expert Review Green"],
            "phenotypes": ["Some condition OMIM:617000"],
        }]

    def test_hgnc_terms_by_pk(self):
        self.assertEqual({32331: self.hgnc_term}, hgnc_ontology_terms_by_pk([32331, 99999999]))

    def test_symbol_that_cannot_be_matched_raises_without_hgnc_term(self):
        """ Establishes the failure this bypasses - the symbol alone has no HGNC to match """
        with self.assertRaises(ValueError):
            OntologyTerm.get_gene_symbol("UNMATCHABLE_SYMBOL")

    def test_relation_created_via_hgnc_term(self):
        _update_gene_relations("UNMATCHABLE_SYMBOL",
                               results_fetcher=lambda _symbol: self._panel_record(),
                               hgnc_term=self.hgnc_term)

        relations = OntologyTermRelation.objects.filter(dest_term_id=self.hgnc_term.id,
                                                        relation=OntologyRelation.PANEL_APP_AU)
        self.assertEqual(1, relations.count())
        relation = relations.get()
        self.assertEqual("OMIM:617000", relation.source_term_id)
        self.assertEqual("Expert Review Green", relation.extra["strongest_classification"])
