from unittest import TestCase

from classification.models import ConditionResolved
from classification.models.clinvar_export_prepare import ConditionGrouper, ConditionGroup
from ontology.models import OntologyTerm, OntologyRelation
from ontology.ontology_builder import OntologyBuilder


class TestClinVarExportModels(TestCase):

    def setup(self):
        ontology_builder = OntologyBuilder(filename="mock", context="mock", import_source="mock")
        ontology_builder.add_term(term_id="MONDO:0000001", name="Disease or Disorder")
        ontology_builder.add_term(term_id="MONDO:0000022", name="Digit Issues")
        ontology_builder.add_term(term_id="MONDO:0000033", name="Toe Issues")
        ontology_builder.add_term(term_id="MONDO:0000044", name="Big Toe Broken")
        ontology_builder.add_term(term_id="MONDO:0000500", name="Bad Lung")
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000044", dest_term_id="MONDO:0000033", relation=OntologyRelation.IS_A)
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000033", dest_term_id="MONDO:0000022", relation=OntologyRelation.IS_A)
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000022", dest_term_id="MONDO:0000001", relation=OntologyRelation.IS_A)
        ontology_builder.add_ontology_relation(source_term_id="MONDO:0000500", dest_term_id="MONDO:0000001", relation=OntologyRelation.IS_A)
        ontology_builder.complete()

    def test_grouping(self):
        def best_candidate(num1: int, num2: int):
            return max(num1, num2)

        condition_num_grouper: ConditionGrouper[int] = ConditionGrouper(best_candidate)
        condition_num_grouper.add_group(ConditionGroup(
            conditions=ConditionResolved(terms=[OntologyTerm.get_or_stub_cached("MONDO:0000033")])),
            candidate=1
        )
        condition_num_grouper.add_group(ConditionGroup(
            conditions=ConditionResolved(terms=[OntologyTerm.get_or_stub_cached("MONDO:0000022")])),
            candidate=2
        )
        condition_num_grouper.add_group(ConditionGroup(
            conditions=ConditionResolved(terms=[OntologyTerm.get_or_stub_cached("MONDO:0000500")])),
            candidate=3
        )