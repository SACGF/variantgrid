from typing import Optional

from django.test import TestCase

from classification.enums import AlleleOriginBucket, ShareLevel, SubmissionSource
from classification.models import (
    AlleleGrouping,
    AlleleOriginGrouping,
    Classification,
    ClassificationGrouping,
    ClassificationModification,
)
from classification.tests.models.test_utils import ClassificationTestUtils
from library.utils import strip_json
from snpdb.models import Allele


class ClassificationGroupingCountsTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.record_count = 0

    def _grouping(self, clinical_significance: Optional[str]) -> ClassificationGrouping:
        self.record_count += 1
        allele_origin_grouping = AlleleOriginGrouping.objects.create(
            allele_grouping=AlleleGrouping.objects.create(allele=Allele.objects.create()),
            allele_origin_bucket=AlleleOriginBucket.GERMLINE
        )
        classification = Classification.objects.create(
            lab=self.lab,
            user=self.user,
            lab_record_id=f"test_{self.record_count}",
            summary=strip_json({"pathogenicity": {"classification": clinical_significance}})
        )
        modification = ClassificationModification.objects.create(
            classification=classification,
            user=self.user,
            source=SubmissionSource.API,
            is_last_published=True
        )
        return ClassificationGrouping.objects.create(
            allele_origin_grouping=allele_origin_grouping,
            lab=self.lab,
            allele_origin_bucket=AlleleOriginBucket.GERMLINE,
            share_level=ShareLevel.ALL_USERS,
            latest_classification_modification=modification
        )

    def test_vus_sub_levels_merge_and_no_data_sorts_last(self):
        for clinical_significance in ["B", "VUS", "VUS_A", "VUS_B", "P", None]:
            self._grouping(clinical_significance)

        counts = ClassificationGrouping.clinical_significance_counts(ClassificationGrouping.objects.all())
        self.assertEqual([(count.clinical_significance, count.count) for count in counts],
                         [("P", 1), ("VUS", 3), ("B", 1), (None, 1)])
        self.assertEqual([count.css_class for count in counts],
                         ["cs-p", "cs-vus", "cs-b", "cs-none"])
        self.assertEqual(counts[-1].label, "No Data")

    def test_oncogenic_sorts_with_germline_equivalent(self):
        for clinical_significance in ["LB", "O", "LO", "LP", "P"]:
            self._grouping(clinical_significance)

        counts = ClassificationGrouping.clinical_significance_counts(ClassificationGrouping.objects.all())
        self.assertEqual([count.clinical_significance for count in counts], ["P", "O", "LP", "LO", "LB"])

    def test_clinical_significance_q_matches_its_count(self):
        for clinical_significance in ["VUS", "VUS_A", "VUS_C", "P", None]:
            self._grouping(clinical_significance)

        for clinical_significance, expected in [("VUS", 3), ("P", 1),
                                                (ClassificationGrouping.NO_CLINICAL_SIGNIFICANCE, 1)]:
            filtered = ClassificationGrouping.objects.filter(
                ClassificationGrouping.clinical_significance_q(clinical_significance))
            self.assertEqual(filtered.count(), expected, clinical_significance)
