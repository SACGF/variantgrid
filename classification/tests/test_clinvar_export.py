from dataclasses import dataclass
from typing import List, Set, Optional
from unittest import TestCase

from classification.models import ConditionResolved
from classification.models.clinvar_export_prepare import ConsolidatingMerger
from classification.tests.data_utils import ConditionMock
from ontology.models import OntologyTerm


@dataclass
class MockEstablished:
    unique_id: int
    condition: ConditionResolved
    candidate: Optional[int]

    def __eq__(self, other):
        if isinstance(other, MockEstablished):
            return self.condition == other.condition and self.candidate == other.candidate and self.unique_id == other.unique_id
        return False

    def __hash__(self):
        return hash(self.unique_id)


@dataclass
class MockCandidate:
    condition: ConditionResolved
    candidate: int

    def __eq__(self, other):
        if isinstance(other, MockCandidate):
            return self.condition == other.condition and self.candidate == other.candidate
        return False


class ConditionGroupPrepareMergerTest(ConsolidatingMerger[MockEstablished, MockCandidate]):

    def __init__(self, established_candidates: Set[MockEstablished]):
        super().__init__()
        self._established_candidates = established_candidates
        self._new_established: List[MockEstablished] = list()
        self.next_id = 100

    def retrieve_established(self) -> Set[MockEstablished]:
        return set(self._established_candidates)

    def establish_new_candidate(self, new_candidate: MockCandidate) -> MockEstablished:
        me = MockEstablished(unique_id=self.next_id, condition=new_candidate.condition, candidate=new_candidate.candidate)
        self.next_id += 1
        self._new_established.append(me)
        return me

    def combine_candidates_if_possible(
            self,
            candidate_1: MockCandidate,
            candidate_2: MockCandidate) -> Optional[MockCandidate]:
        if general_condition := ConditionResolved.more_general_term_if_related(candidate_1.condition,
                                                                               candidate_2.condition):
            return MockCandidate(condition=general_condition, candidate=max(candidate_1.candidate, candidate_2.candidate))
        else:
            return None

    def merge_into_established_if_possible(self, established: MockEstablished, new_candidate: Optional[MockCandidate]) -> bool:
        if new_candidate is None:
            established.candidate = None
            return True
        elif new_candidate.condition.is_same_or_more_specific(established.condition):
            established.candidate = new_candidate.candidate
            return True
        else:
            return False


class TestClinVarExportModels(TestCase):

    def setUp(self):
        ConditionMock.setUp()

    def test_grouping(self):
        # m_disease = OntologyTerm.get_or_stub(ConditionMock.MONDO_DISEASE_OR_DISORDER)
        # m_digit = OntologyTerm.get_or_stub(ConditionMock.MONDO_DIGIT_ISSUE)
        m_toe = OntologyTerm.get_or_stub(ConditionMock.MONDO_TOE_ISSUE)
        m_big_toe = OntologyTerm.get_or_stub(ConditionMock.MONDO_BIG_TOE_BROKEN)
        m_bad_lung = OntologyTerm.get_or_stub(ConditionMock.MONDO_BAD_LUNG)
        m_bad_heart = OntologyTerm.get_or_stub(ConditionMock.MONDO_BAD_HEART)

        simple_big_toe = ConditionResolved(terms=[m_big_toe])
        simple_toe = ConditionResolved(terms=[m_toe])
        simple_bad_lung = ConditionResolved(terms=[m_bad_lung])
        simple_bad_heart = ConditionResolved(terms=[m_bad_heart])
        # simple_o_big_toe = ConditionResolved(terms=[o_big_toe])

        group_toe = MockCandidate(
            condition=simple_toe,
            candidate=1
        )
        group_big_toe = MockCandidate(
            condition=simple_big_toe,
            candidate=2
        )
        group_lung = MockCandidate(
            condition=simple_bad_lung,
            candidate=3
        )

        # test all new groups
        established = MockEstablished(unique_id=1, condition=simple_bad_heart, candidate=7)
        condition_num_grouper = ConditionGroupPrepareMergerTest(established_candidates={established, })

        # these two groups should merge
        condition_num_grouper.add_new_candidate(group_toe)
        self.assertEqual(len(condition_num_grouper.collapsed_new_candidates), 1)
        condition_num_grouper.add_new_candidate(group_big_toe)
        self.assertEqual(len(condition_num_grouper.collapsed_new_candidates), 1)

        # this group shouldn't merge
        condition_num_grouper.add_new_candidate(group_lung)
        self.assertEqual(len(condition_num_grouper.collapsed_new_candidates), 2)

        # no retrieve_established groups, so we should now have 2 retrieve_established groups
        condition_num_grouper.consolidate()
        new_groups = condition_num_grouper._new_established
        self.assertEqual(len(new_groups), 2)
        # takes the more general term from group_toe,

        self.assertEqual(new_groups[0].condition, simple_toe)
        self.assertEqual(new_groups[0].candidate, 2)

        self.assertEqual(new_groups[1].condition, simple_bad_lung)
        self.assertEqual(new_groups[1].candidate, 3)

        # no new candidate for bad lung, candidate value should be updated to 0
        self.assertIsNone(established.candidate)

        # TODO test merging and non merging
