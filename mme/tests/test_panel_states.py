""" The MatchMaker Exchange card renders for every classification on an MME deployment and
names which layer of the opt-in is unmet - a card that vanishes on an ineligible record
teaches a curator nothing. """
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from classification.enums.classification_enums import ClinicalSignificance, ShareLevel
from classification.tests.models.test_utils import ClassificationTestUtils
from library.django_utils.unittest_utils import URLTestCase
from mme.models import MMEInboundMatch, MMEInboundQuery
from mme.tests.fakes import make_classification

NODES = {"testnode": {"base_url": "https://node.test", "token": "secret-token", "api_version": "1.1"}}


@override_settings(MME_ENABLED=True, MME_NODES=NODES,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class PanelStateTestCase(URLTestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()
        self.client.force_login(self.user)

    def _get(self, classification, **params):
        url = reverse("mme_classification_panel",
                      kwargs={"classification_id": classification.pk})
        return self.client.get(url, params)

    def test_eligible_shows_the_rule_on_the_submit_button(self):
        response = self._get(make_classification(self.lab, self.user))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Share with MatchMaker Exchange (VUS and above)")
        self.assertContains(response, "testnode")

    def test_benign_names_clinvar_as_the_right_channel(self):
        response = self._get(make_classification(
            self.lab, self.user, clinical_significance=ClinicalSignificance.BENIGN))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "candidate genes")
        self.assertContains(response, "Benign")
        self.assertContains(response, "ClinVar")
        self.assertNotContains(response, "Share with MatchMaker Exchange (VUS and above)")

    def test_unclassified_record_is_also_not_a_candidate(self):
        response = self._get(make_classification(
            self.lab, self.user, clinical_significance=None))
        self.assertContains(response, "ClinVar")

    def test_not_shared_publicly_names_the_share_level(self):
        response = self._get(make_classification(
            self.lab, self.user, share_level=ShareLevel.ALL_USERS.value))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Share this classification at")
        self.assertContains(response, "3rd Party Databases")

    def test_lab_not_opted_in_names_the_lab_admin(self):
        classification = make_classification(self.lab, self.user)
        self.lab.mme_enabled = False
        self.lab.save()
        response = self._get(classification)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "lab admin")

    def test_fragment_and_standalone_render_the_same_card(self):
        classification = make_classification(self.lab, self.user)
        standalone = self._get(classification).content.decode()
        fragment = self._get(classification, fragment=1).content.decode()

        # The card's own markup is identical; only the page chrome differs.
        for marker in ("Share with MatchMaker Exchange (VUS and above)",
                       "Genomic features", "Phenotype features (HPO)"):
            self.assertIn(marker, fragment)
            self.assertIn(marker, standalone)

        # ?fragment=1 skips the base template - it drops into the existing MME card body
        self.assertNotIn("menu", fragment)
        self.assertLess(len(fragment), len(standalone))

    def test_inbound_matches_are_listed_and_link_to_the_match(self):
        classification = make_classification(self.lab, self.user)
        inbound_query = MMEInboundQuery.objects.create(
            peer_node_id="genematcher", request_json={}, num_results=1)
        match = MMEInboundMatch.objects.create(
            inbound_query=inbound_query, classification=classification, score=0.5,
            remote_patient_id="remote-1", query_patient_json={})

        response = self._get(classification)
        self.assertContains(response, "genematcher")
        self.assertContains(response, "remote-1")
        self.assertContains(response,
                            reverse("mme_view_inbound_match", kwargs={"inbound_match_id": match.pk}))

    def test_user_without_view_permission_is_refused(self):
        # Not shared, and the stranger is in no group - no write, no readable modification
        classification = make_classification(self.lab, self.user,
                                             share_level=ShareLevel.LAB.value)
        stranger = User.objects.create_user("mme_stranger", "stranger@nowhere.org", "pw")
        self.client.logout()
        self.client.force_login(stranger)
        response = self._get(classification)
        self.assertIn(response.status_code, (302, 403, 404))


@override_settings(MME_ENABLED=True, MME_NODES=NODES,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class InboundMatchPageTestCase(URLTestCase):
    """ The notification's direct link shows everything the peer gave us. """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()
        self.client.force_login(self.user)

        classification = make_classification(self.lab, self.user)
        inbound_query = MMEInboundQuery.objects.create(
            peer_node_id="genematcher", request_json={}, num_results=1)
        self.match = MMEInboundMatch.objects.create(
            inbound_query=inbound_query, classification=classification, score=0.5,
            remote_patient_id="remote-1",
            query_patient_json={"id": "remote-1",
                                "contact": {"name": "Dr Peer", "href": "mailto:peer@peer.org"}})

    def test_shows_the_peer_their_contact_and_their_patient(self):
        response = self.client.get(
            reverse("mme_view_inbound_match", kwargs={"inbound_match_id": self.match.pk}))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "genematcher")
        self.assertContains(response, "Dr Peer")
        self.assertContains(response, "peer@peer.org")
