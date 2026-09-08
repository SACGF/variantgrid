from django.test import TestCase

from classification.enums import AlleleOriginBucket, ShareLevel, SubmissionSource
from classification.models import (
    AlleleGrouping,
    AlleleOriginGrouping,
    Classification,
    ClassificationGrouping,
    ClassificationGroupingEntry,
    ClassificationModification,
    ClinicalContext,
)
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.models import Allele, AlleleConversionTool


class AlleleMergeRehomingTestCase(TestCase):
    """ Allele.merge() moves Classification.allele in bulk - the clinical context and grouping it is filed
        under are derived records that stay behind unless the merge signal re-homes them (#1361) """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _classification_on(self, allele: Allele) -> Classification:
        classification = Classification.objects.create(
            lab=self.lab,
            user=self.user,
            lab_record_id="merge_rehome_test",
            allele=allele,
            allele_origin_bucket=AlleleOriginBucket.GERMLINE,
            share_level=ShareLevel.ALL_USERS,
        )
        modification = ClassificationModification.objects.create(
            classification=classification,
            user=self.user,
            source=SubmissionSource.API,
            is_last_published=True
        )
        classification.clinical_context = ClinicalContext.objects.create(
            allele=allele,
            allele_origin_bucket=AlleleOriginBucket.GERMLINE,
            name=ClinicalContext.default_name
        )
        classification.save()

        allele_origin_grouping = AlleleOriginGrouping.objects.create(
            allele_grouping=AlleleGrouping.objects.create(allele=allele),
            allele_origin_bucket=AlleleOriginBucket.GERMLINE
        )
        grouping = ClassificationGrouping.objects.create(
            allele_origin_grouping=allele_origin_grouping,
            lab=self.lab,
            allele_origin_bucket=AlleleOriginBucket.GERMLINE,
            share_level=ShareLevel.ALL_USERS,
            latest_classification_modification=modification
        )
        ClassificationGroupingEntry.objects.create(classification=classification, grouping=grouping)
        return classification

    def test_allele_merge_rehomes_classification(self):
        allele_a = Allele.objects.create()
        allele_b = Allele.objects.create()
        classification = self._classification_on(allele_b)

        self.assertTrue(allele_a.merge(AlleleConversionTool.BCFTOOLS_LIFTOVER, allele_b))

        classification.refresh_from_db()
        entry = ClassificationGroupingEntry.objects.get(classification=classification)
        self.assertEqual(classification.allele, allele_a)
        self.assertEqual(classification.clinical_context.allele, allele_a)
        self.assertEqual(entry.grouping.allele_origin_grouping.allele_grouping.allele, allele_a)
