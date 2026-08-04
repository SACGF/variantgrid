from django.contrib.auth.models import User
from django.test import override_settings

from classification.tests.models.test_utils import ClassificationTestUtils
from library.django_utils.unittest_utils import URLTestCase
from mme.models import MMESubmission, MMEInboundQuery, MMEInboundMatch
from mme.tests.fakes import make_classification

NODES = {"testnode": {"base_url": "https://node.test", "token": "secret-token", "api_version": "1.1",
                      "disclaimer": "baseline disclaimer"}}


@override_settings(MME_ENABLED=True, MME_NODES=NODES,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class Test(URLTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        ClassificationTestUtils.setUp()
        cls.lab, cls.user = ClassificationTestUtils.lab_and_user()
        cls.lab.contact_email = "curator@lab.org"
        cls.lab.mme_enabled = True
        cls.lab.save()
        cls.user_admin = User.objects.create_superuser("mme_url_admin", "foo@bar", "password")

        cls.classification = make_classification(cls.lab, cls.user)
        cls.submission = MMESubmission.objects.create(
            classification=cls.classification, node_id="testnode", external_patient_id="vg:1")
        inbound_query = MMEInboundQuery.objects.create(
            peer_node_id="testpeer", request_json={}, num_results=1)
        cls.inbound_match = MMEInboundMatch.objects.create(
            inbound_query=inbound_query, classification=cls.classification, score=0.5,
            remote_patient_id="remote-1", query_patient_json={})

        cls.PRIVATE_URL_NAMES_AND_KWARGS = [
            ("mme_classification_panel", {"classification_id": cls.classification.pk}, 200),
            ("mme_view_submission", {"submission_id": cls.submission.pk}, 200),
            ("mme_view_inbound_match", {"inbound_match_id": cls.inbound_match.pk}, 200),
        ]
        # Service requirements #9 and #7 oblige us to publish these without a login.
        cls.PUBLIC_URL_NAMES_AND_KWARGS = [
            ("mme_public_metrics", {}, 200),
            ("mme_public_disclaimers", {}, 200),
        ]

    def test_urls_for_owner(self):
        self._test_urls(self.PRIVATE_URL_NAMES_AND_KWARGS, self.user)

    def test_urls_for_admin(self):
        self._test_urls(self.PRIVATE_URL_NAMES_AND_KWARGS, self.user_admin)

    def test_public_urls_anonymous(self):
        self._test_urls(self.PUBLIC_URL_NAMES_AND_KWARGS)
