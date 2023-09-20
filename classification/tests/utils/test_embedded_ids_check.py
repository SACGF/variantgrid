from unittest import TestCase

from classification.models import embedded_ids_check
from classification.tests.utils.data_utils import ConditionMock


class Embedded_Ids_Check(TestCase):

    def setUp(self):
        ConditionMock.setUp()

    def test_valid_ontology_id(self):
        ids = f'{ConditionMock.MONDO_TOE_ISSUE}{ConditionMock.MONDO_BIG_TOE_BROKEN};uncertain '
        self._run_test(ids, ['MONDO:9000033', 'MONDO:9000044'], 'U', 5)

    def test_multiple_ontology_ids(self):
        ids = f'{ConditionMock.MONDO_TOE_ISSUE}{ConditionMock.MONDO_DIGIT_ISSUE};co-occurring'
        self._run_test(ids, ['MONDO:9000033', 'MONDO:9000022'], 'C', 5)

    def test_single_ontology_id(self):
        ids = f'{ConditionMock.MONDO_BAD_HEART}'
        self._run_test(ids, ['MONDO:9000600'], 'N', 5)

    def test_multiple_ontology_ids_2(self):
        ids = f'{ConditionMock.MONDO_BAD_LUNG}{ConditionMock.MONDO_BIG_TOE_BROKEN};cooccurring'
        self._run_test(ids, ['MONDO:9000500', 'MONDO:9000044'], 'C', 5)

    def test_multiple_ontology_ids_with_text(self):
        ids = 'MONDO:0010726 Rett syndrome, MONDO:0013249 autosomal recessive nonsyndromic hgjghj jhjhjds srtghfdv  hearing loss 84A; UNCERTAIN'
        self._run_test(ids, ['MONDO:0010726', 'MONDO:0013249'], 'U', 7)

    def _run_test(self, ids, expected_ids, expected_condition, expected_messages_length):
        result = embedded_ids_check(ids)
        self.assertTrue(result.ids_found_in_text)
        self.assertTrue(result.validated)

        if len(result.terms) == len(expected_ids):
            self.assertEqual({term.id for term in result.terms}, set(expected_ids))
            self.assertEqual(result.condition_multi_operation, expected_condition)
            self.assertEqual(len(result.messages), expected_messages_length)
        else:
            self.assertEqual(result.terms[0].id, expected_ids[0])
            self.assertEqual(result.condition_multi_operation, expected_condition)
