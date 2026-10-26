from django.test import TestCase

from annotation.enums import CitationSource
from annotation.models.models_citations import CitationIdNormalized
from classification.enums import ShareLevel, SubmissionSource
from classification.models.classification import Classification
from classification.tasks.classification_citation_tasks import classification_citation_ids
from classification.tests.models.test_utils import ClassificationTestUtils


class ClassificationCitationsTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def _create(self, literature: str, publish: bool) -> Classification:
        lab, user = ClassificationTestUtils.lab_and_user()
        classification = Classification.create(user=user, lab=lab, lab_record_id=None,
                                               data={'literature': literature},
                                               save=True, source=SubmissionSource.API, make_fields_immutable=False)
        if publish:
            classification.publish_latest(user=user, share_level=ShareLevel.LAB)
        return classification

    def test_classification_citation_ids(self):
        republished = self._create("PMID: 11111111", publish=True)
        self._create("PMCID: PMC2222222 and OMIM:300000", publish=True)
        self._create("PMID: 33333333", publish=False)

        # only the last published version is exported
        _lab, user = ClassificationTestUtils.lab_and_user()
        republished.patch_value({'literature': "PMID: 44444444"}, user=user, source=SubmissionSource.API, save=True)
        republished.publish_latest(user=user, share_level=ShareLevel.LAB)

        self.assertEqual(classification_citation_ids(), {
            CitationIdNormalized(CitationSource.PUBMED, "44444444"),
            CitationIdNormalized(CitationSource.PUBMED_CENTRAL, "PMC2222222"),
        })
