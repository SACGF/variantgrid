from django.test import TestCase

from classification.enums import AlleleOriginBucket, SubmissionSource
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationModification,
)
from classification.tests.models.test_utils import ClassificationTestUtils


class ClassificationConsensusPatchTestCase(TestCase):
    """ What copy_scope and copy_allele_origin let through when copying from a previous classification """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _consensus_patch(self, allele_origin_bucket: AlleleOriginBucket) -> dict:
        classification = Classification.objects.create(
            lab=self.lab,
            user=self.user,
            lab_record_id=f"consensus_{allele_origin_bucket}",
            allele_origin_bucket=allele_origin_bucket
        )
        modification = ClassificationModification.objects.create(
            classification=classification,
            user=self.user,
            source=SubmissionSource.API,
            is_last_published=True,
            published_evidence={
                "condition": {"value": "MONDO:0007947"},  # ALLELE
                "h_summary": {"value": "gene summary"},  # GENE
                "segregation": {"value": "co-segregates"},  # GERMLINE
                "somatic:tmb_value": {"value": 12},  # NONE - this patient's tumour
            }
        )
        return ClassificationConsensus(modification=modification).consensus_patch

    def test_germline_source_copies_everything_but_the_tumour_measurement(self):
        patch = self._consensus_patch(AlleleOriginBucket.GERMLINE)
        self.assertEqual(patch["condition"], {"value": "MONDO:0007947"})
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertEqual(patch["segregation"], {"value": "co-segregates"})
        self.assertNotIn("somatic:tmb_value", patch)

    def test_somatic_source_leaves_germline_only_keys_behind(self):
        patch = self._consensus_patch(AlleleOriginBucket.SOMATIC)
        self.assertEqual(patch["condition"], {"value": "MONDO:0007947"})
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertNotIn("segregation", patch)
        self.assertNotIn("somatic:tmb_value", patch)
