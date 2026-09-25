from django.test import TestCase

from classification.enums import SpecialEKeys
from classification.models import Classification, ClassificationModification, EvidenceKeyMap


class ClinicalSignificanceCountsSummaryTestCase(TestCase):
    """ The "Compare with other classifications for this variant" counts on the classification page """

    @staticmethod
    def _record(evidence: dict) -> ClassificationModification:
        return ClassificationModification(published_evidence={key: {"value": value} for key, value in evidence.items()})

    def test_somatic_tiers_are_counted_alongside_germline(self):
        records = [
            self._record({SpecialEKeys.CLINICAL_SIGNIFICANCE: "P"}),
            self._record({SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE: "tier_1"}),
            self._record({SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE: "tier_1"}),
            self._record({SpecialEKeys.CLINICAL_SIGNIFICANCE: "VUS", SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE: "tier_3"}),
            self._record({SpecialEKeys.CLINICAL_SIGNIFICANCE: None}),
        ]
        germline = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value
        somatic = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value
        expected = f"{somatic('tier_1')} x2, {germline('P')} x1, {germline('VUS')} x1, {somatic('tier_3')} x1, Unclassified x1"
        self.assertEqual(Classification.clinical_significance_counts_summary(records), expected)
