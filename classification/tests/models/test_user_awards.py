from datetime import timedelta

from django.test import TestCase
from django.utils import timezone

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import Classification, ClassificationModification
from classification.tests.models.test_utils import ClassificationTestUtils
from classification.user_awards import _cold_cases


class ColdCaseAwardTest(TestCase):
    """ cold_case counts a published modification whose previous modification is > 6 months older - #1819 """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _classification_with_gap(self, gap: timedelta) -> Classification:
        classification = Classification.create(user=self.user, lab=self.lab,
                                               data={SpecialEKeys.C_HGVS: {'value': 'c.301A>C'}},
                                               save=True, source=SubmissionSource.API)
        classification.publish_latest(user=self.user)
        classification.patch_value(patch={SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'P'}},
                                   clear_all_fields=False, user=self.user, source=SubmissionSource.API, save=True)
        classification.publish_latest(user=self.user)

        modifications = list(ClassificationModification.objects.filter(classification=classification).order_by("created"))
        self.assertGreaterEqual(len(modifications), 2)
        ClassificationModification.objects.filter(pk=modifications[0].pk).update(created=timezone.now() - gap)
        return classification

    def test_year_gap_counts(self):
        self._classification_with_gap(timedelta(days=365))
        self.assertEqual(_cold_cases(None), {None: {self.user.pk: 1}})

    def test_month_gap_ignored(self):
        self._classification_with_gap(timedelta(days=30))
        self.assertEqual(_cold_cases(None), {})
