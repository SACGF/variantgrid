from datetime import date
from typing import Optional

from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

from classification.enums import ShareLevel, SpecialEKeys, SubmissionSource
from classification.models import (
    AlleleOriginGrouping,
    Classification,
    ClassificationGrouping,
    ClassificationGroupingEntry,
    ClassificationModification,
)
from classification.models.classification_grouping import ClassificationGroupingPathogenicDifference
from classification.tests.models.test_utils import ClassificationTestUtils
from library.utils import strip_json
from snpdb.models import Allele
from snpdb.models.models_enums import AlleleOriginBucket


class ClassificationGroupingCountsTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.record_count = 0

    def _grouping(self, clinical_significance: Optional[str]) -> ClassificationGrouping:
        self.record_count += 1
        allele_origin_grouping = AlleleOriginGrouping.objects.create(
            allele=Allele.objects.create(),
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


class ClassificationGroupingSortColumnsTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _grouping(self, summary: dict) -> ClassificationGrouping:
        classification = Classification.objects.create(lab=self.lab, user=self.user, lab_record_id="test_1",
                                                       summary=summary)
        ClassificationModification.objects.create(classification=classification, user=self.user,
                                                  source=SubmissionSource.API, is_last_published=True)
        grouping = ClassificationGrouping.objects.create(
            allele_origin_grouping=AlleleOriginGrouping.objects.create(allele=Allele.objects.create(),
                                                                       allele_origin_bucket=AlleleOriginBucket.SOMATIC),
            lab=self.lab,
            share_level=ShareLevel.ALL_USERS
        )
        ClassificationGroupingEntry.objects.create(classification=classification, grouping=grouping)
        grouping.update()
        return grouping

    def test_update_copies_sort_columns_and_clears_them_when_empty(self):
        grouping = self._grouping({"pathogenicity": {"classification": "P", "sort": 5},
                                   "somatic": {"clinical_significance": "tier_1", "sort": 3},
                                   "date": {"date": "2024-03-07", "type": "curation_date"}})
        self.assertEqual((grouping.latest_pathogenicity_sort, grouping.latest_somatic_sort, grouping.latest_curated_date),
                         (5, 3, date(2024, 3, 7)))

        grouping.classificationgroupingentry_set.all().delete()
        del grouping.classification_modifications
        grouping.update()
        grouping.refresh_from_db()
        self.assertEqual((grouping.latest_pathogenicity_sort, grouping.latest_somatic_sort, grouping.latest_curated_date),
                         (None, None, None))

    def test_malformed_curated_date_is_null(self):
        for date_value in ["", "not a date", None]:
            ClassificationGrouping.objects.all().delete()
            Classification.objects.all().delete()
            grouping = self._grouping({"date": {"date": date_value}})
            self.assertIsNone(grouping.latest_curated_date, date_value)


class ClassificationGroupingUpdateQueriesTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _grouping(self, record_count: int) -> ClassificationGrouping:
        grouping = ClassificationGrouping.objects.create(
            allele_origin_grouping=AlleleOriginGrouping.objects.create(allele=Allele.objects.create(),
                                                                       allele_origin_bucket=AlleleOriginBucket.GERMLINE),
            lab=self.lab,
            share_level=ShareLevel.ALL_USERS
        )
        for i in range(record_count):
            classification = Classification.objects.create(
                lab=self.lab, user=self.user, lab_record_id=f"test_{grouping.pk}_{i}",
                summary={"pathogenicity": {"classification": ["P", "LP", "VUS"][i % 3], "sort": i}})
            ClassificationModification.objects.create(
                classification=classification, user=self.user, source=SubmissionSource.API, is_last_published=True,
                published_evidence={SpecialEKeys.ZYGOSITY: {"value": ["heterozygous", "homozygous"][i % 2]},
                                    SpecialEKeys.CONDITION: {"value": f"condition {i % 2}"}})
            ClassificationGroupingEntry.objects.create(classification=classification, grouping=grouping)
        return grouping

    @staticmethod
    def _update_queries(grouping: ClassificationGrouping) -> list[str]:
        grouping = ClassificationGrouping.objects.get(pk=grouping.pk)
        with CaptureQueriesContext(connection) as ctx:
            grouping.update()
        return [q["sql"] for q in ctx.captured_queries]

    def test_update_queries_do_not_grow_with_classifications(self):
        self._update_queries(self._grouping(1))  # Fill the evidence key cache
        single_queries = self._update_queries(self._grouping(1))
        several = self._grouping(6)
        several_queries = self._update_queries(several)
        self.assertEqual(len(single_queries), len(several_queries))
        # SAVEPOINT/RELEASE, the last published modifications (one query, select_related), search terms
        # (allele info + resolved variant info, ClinVar SCVs, DELETE), the grouping UPDATE, and the post_save
        # refresh of its OverlapContribution (contributions + allele origin grouping)
        self.assertLessEqual(len(several_queries), 10)

        several.refresh_from_db()
        self.assertEqual(several.classification_count, 6)
        self.assertEqual(several.pathogenic_difference, ClassificationGroupingPathogenicDifference.CLIN_SIG_DIFFS)
        self.assertEqual(several.zygosity_values, ["heterozygous", "homozygous"])
