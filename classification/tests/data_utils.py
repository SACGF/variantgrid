from ontology.models import OntologyRelation
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException


class ConditionMock:

    MONDO_DISEASE_OR_DISORDER = "MONDO:0000001"
    MONDO_DIGIT_ISSUE = "MONDO:9000022"
    MONDO_TOE_ISSUE = "MONDO:9000033"
    MONDO_BIG_TOE_BROKEN = "MONDO:9000044"
    MONDO_LITTLE_TOE_BORKEN = "MONDO:9000045"
    MONDO_BAD_LUNG = "MONDO:9000500"
    MONDO_BAD_HEART = "MONDO:9000600"

    OMIM_BIG_TOE_BROKEN = "OMIM:9000044"

    @staticmethod
    def setUp():
        ontology_builder = OntologyBuilder(filename="mock", context="mock", import_source="mock")
        try:
            ontology_builder.ensure_hash_changed("3")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_DISEASE_OR_DISORDER, name="Disease or Disorder")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_DIGIT_ISSUE, name="Digit Issues")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_TOE_ISSUE, name="Toe Issues")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_BIG_TOE_BROKEN, name="Big Toe Broken")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_LITTLE_TOE_BORKEN, name="Little Toe Broken")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_BAD_LUNG, name="Bad Lung")
            ontology_builder.add_term(term_id=ConditionMock.MONDO_BAD_HEART, name="Bad Heart")
            ontology_builder.add_term(term_id=ConditionMock.OMIM_BIG_TOE_BROKEN, name="OMIM Big Toe")
            ontology_builder.add_ontology_relation(source_term_id=ConditionMock.MONDO_BIG_TOE_BROKEN, dest_term_id=ConditionMock.MONDO_TOE_ISSUE, relation=OntologyRelation.IS_A)
            ontology_builder.add_ontology_relation(source_term_id=ConditionMock.MONDO_TOE_ISSUE, dest_term_id=ConditionMock.MONDO_DIGIT_ISSUE, relation=OntologyRelation.IS_A)
            ontology_builder.add_ontology_relation(source_term_id=ConditionMock.MONDO_DIGIT_ISSUE, dest_term_id=ConditionMock.MONDO_DISEASE_OR_DISORDER, relation=OntologyRelation.IS_A)
            ontology_builder.add_ontology_relation(source_term_id=ConditionMock.MONDO_BAD_LUNG, dest_term_id=ConditionMock.MONDO_DISEASE_OR_DISORDER, relation=OntologyRelation.IS_A)
            ontology_builder.add_ontology_relation(source_term_id=ConditionMock.MONDO_BAD_HEART, dest_term_id=ConditionMock.MONDO_DISEASE_OR_DISORDER, relation=OntologyRelation.IS_A)
            ontology_builder.add_ontology_relation(source_term_id=ConditionMock.MONDO_BIG_TOE_BROKEN, dest_term_id=ConditionMock.OMIM_BIG_TOE_BROKEN, relation=OntologyRelation.EXACT)
            ontology_builder.complete()
        except OntologyBuilderDataUpToDateException:
            pass