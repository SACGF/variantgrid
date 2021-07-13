from unittest import TestCase
from ontology.models import OntologyTerm, OntologyService, OntologyRelation
from ontology.ontology_builder import OntologyBuilder


class TestClinVarExportModels(TestCase):

    def setup(self):
        ontology_builder = OntologyBuilder(filename="mock", context="mock", import_source="mock")
        ontology_builder.add_term(term_id="MONDO:0000001", name="Disease or Disorder")
        ontology_builder.add_term(term_id="MONDO:0000022", name="Digit Issues")
        ontology_builder.add_term(term_id="MONDO:0000033", name="Toe Issues")
        ontology_builder.add_term(term_id="MONDO:0000044", name="Big Toe Broken")
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000044", dest_term_id="MONDO:0000033", relation=OntologyRelation.IS_A)
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000033", dest_term_id="MONDO:0000022", relation=OntologyRelation.IS_A)
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000022", dest_term_id="MONDO:0000001", relation=OntologyRelation.IS_A)
        ontology_builder.complete()
