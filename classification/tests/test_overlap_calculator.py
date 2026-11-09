"""
OverlapCalculatorBase.calculate_entries: shared non-comparable values (e.g. Drug Response) contribute to an Overlap
without counting towards its status, judged on the effective value so a triage amending to or from one counts.
"""
from typing import Optional

from django.test import TestCase

from classification.enums.classification_enums import OverlapStatus, TestingContextBucket
from classification.enums.overlaps_enums import (
    ClassificationResultValue,
    OverlapContributionStatus,
    OverlapEntrySourceTextChoices,
    TriageDecision,
    TriageStatus,
)
from classification.models.overlaps_model import OverlapContribution
from classification.services.overlap_calculator import OverlapCalculatorOncPath


def _contribution(value: str, amend_value: Optional[str] = None) -> OverlapContribution:
    triage = TriageDecision(status=TriageStatus.REVIEWED_WILL_FIX, amend_value=amend_value) if amend_value else TriageDecision()
    return OverlapContribution(
        source=OverlapEntrySourceTextChoices.CLASSIFICATION,
        value_type=ClassificationResultValue.ONC_PATH,
        value=value,
        contribution_status=OverlapContributionStatus.CONTRIBUTING,
        testing_context_bucket=TestingContextBucket.GERMLINE,
        triage_state=triage.to_dict(),
    )


class OverlapCalculatorNonComparableTest(TestCase):

    def _status(self, *contributions: OverlapContribution) -> OverlapStatus:
        return OverlapCalculatorOncPath.calculate_entries(contributions).status

    def test_non_comparable_value_ignored(self):
        self.assertEqual(self._status(_contribution("P"), _contribution("drug_response")), OverlapStatus.SINGLE_SUBMITTER)
        self.assertEqual(self._status(_contribution("drug_response")), OverlapStatus.NO_COUNTING_CONTRIBUTIONS)

    def test_amend_to_non_comparable_value(self):
        snapshot = OverlapCalculatorOncPath.calculate_entries([_contribution("P"), _contribution("B", amend_value="drug_response")])
        self.assertEqual(snapshot.status, OverlapStatus.SINGLE_SUBMITTER)
        self.assertTrue(snapshot.has_pending_values)

    def test_amend_from_non_comparable_value(self):
        self.assertEqual(self._status(_contribution("P"), _contribution("drug_response", amend_value="B")),
                         OverlapStatus.MEDICALLY_SIGNIFICANT)
