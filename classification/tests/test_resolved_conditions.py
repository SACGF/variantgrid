from unittest import TestCase

from classification.models import ConditionResolved, MultiCondition
from classification.tests.data_utils import ConditionMock
from ontology.models import OntologyTerm, OntologySnake, OntologyTermRelation


class ResolvedConditionTest(TestCase):

    def setUp(self):
        ConditionMock.setUp()

    def test_ancestor_relationship(self):
        m_disease = OntologyTerm.get_or_stub(ConditionMock.MONDO_DISEASE_OR_DISORDER)
        m_digit = OntologyTerm.get_or_stub(ConditionMock.MONDO_DIGIT_ISSUE)
        m_toe = OntologyTerm.get_or_stub(ConditionMock.MONDO_TOE_ISSUE)
        m_big_toe = OntologyTerm.get_or_stub(ConditionMock.MONDO_BIG_TOE_BROKEN)
        m_bad_lung = OntologyTerm.get_or_stub(ConditionMock.MONDO_BAD_LUNG)
        o_big_toe = OntologyTerm.get_or_stub(ConditionMock.OMIM_BIG_TOE_BROKEN)

        snakes = OntologySnake.check_if_ancestor(m_big_toe, m_disease)
        self.assertEqual(len(snakes), 1)
        self.assertEqual([step.dest_term for step in snakes[0].show_steps()], [m_toe, m_digit, m_disease])

        snakes = OntologySnake.check_if_ancestor(m_big_toe, m_disease, max_levels=2)
        self.assertEqual(len(snakes), 0)

        snakes = OntologySnake.check_if_ancestor(m_bad_lung, m_digit, max_levels=10)
        self.assertEqual(len(snakes), 0)

        self.assertEqual(m_big_toe, OntologyTermRelation.as_mondo(o_big_toe))

    def test_resolved_condition(self):
        # m_disease = OntologyTerm.get_or_stub(ConditionMock.MONDO_DISEASE_OR_DISORDER)
        # m_digit = OntologyTerm.get_or_stub(ConditionMock.MONDO_DIGIT_ISSUE)
        m_toe = OntologyTerm.get_or_stub(ConditionMock.MONDO_TOE_ISSUE)
        m_big_toe = OntologyTerm.get_or_stub(ConditionMock.MONDO_BIG_TOE_BROKEN)
        m_bad_lung = OntologyTerm.get_or_stub(ConditionMock.MONDO_BAD_LUNG)
        o_big_toe = OntologyTerm.get_or_stub(ConditionMock.OMIM_BIG_TOE_BROKEN)

        simple_big_toe = ConditionResolved(terms=[m_big_toe])
        simple_toe = ConditionResolved(terms=[m_toe])
        simple_bad_lung = ConditionResolved(terms=[m_bad_lung])
        simple_o_big_toe = ConditionResolved(terms=[o_big_toe])

        self.assertFalse(simple_toe.is_multi_condition)
        self.assertEqual(m_big_toe, simple_big_toe.mondo_term)

        # when comparing 2 terms in same lineage, make sure more general term is found regardless of order
        self.assertEqual(simple_toe, ConditionResolved.more_general_term_if_related(simple_big_toe, simple_toe))
        self.assertEqual(simple_toe, ConditionResolved.more_general_term_if_related(simple_toe, simple_big_toe))
        # when comparing an item to itself, it should be considered the general term
        self.assertEqual(simple_toe, ConditionResolved.more_general_term_if_related(simple_toe, simple_toe))

        # two unrelated terms should return None
        self.assertIsNone(ConditionResolved.more_general_term_if_related(simple_big_toe, simple_bad_lung))

        complex_combo_1 = ConditionResolved(terms=[m_bad_lung, m_big_toe], join=MultiCondition.CO_OCCURRING)
        complex_combo_2 = ConditionResolved(terms=[m_big_toe, m_toe], join=MultiCondition.CO_OCCURRING)

        # ensure limitations of multi-condition are obeyed
        self.assertTrue(complex_combo_1.is_multi_condition)
        self.assertIsNone(complex_combo_1.mondo_term)  # multi-terms shouldn't return a value for hte single term
        self.assertEqual(complex_combo_1, ConditionResolved.more_general_term_if_related(complex_combo_1, complex_combo_1))  # multi-term should match itself
        self.assertIsNone(ConditionResolved.more_general_term_if_related(complex_combo_1, simple_toe))  # multi-terms can only match themselves
        self.assertIsNone(ConditionResolved.more_general_term_if_related(complex_combo_1, complex_combo_2))  # multi-terms can only match themselves

        # test cross ontology, should favor Mondo
        self.assertEqual(simple_big_toe, ConditionResolved.more_general_term_if_related(simple_big_toe, simple_o_big_toe))
        self.assertEqual(simple_big_toe, ConditionResolved.more_general_term_if_related(simple_o_big_toe, simple_big_toe))
