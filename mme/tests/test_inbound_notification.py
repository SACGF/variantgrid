""" The depositor - us - is notified when a peer's query matches one of our records, not
just the querier. Peers re-query on a schedule, so repeats must stay silent. """
from unittest.mock import patch

from django.test import TestCase, override_settings

from classification.tests.models.test_utils import ClassificationTestUtils
from mme.models import MMEInboundQuery, MMEInboundMatch
from mme.notifications import record_inbound_matches, notify_depositors
from mme.tests.fakes import make_classification
from user_messages.models import Message

QUERY_PATIENT = {"id": "remote-1", "contact": {"name": "Dr Peer", "href": "mailto:peer@peer.org"},
                 "genomicFeatures": [{"gene": {"id": "BRCA1"}}]}


class _Match:
    """ Stand-in for mme.matching.MMEMatch - record_inbound_matches only reads these. """

    def __init__(self, classification, score=0.75):
        self.classification = classification
        self.score = score
        self.patient = {}


@override_settings(MME_ENABLED=True, MME_FROM_EMAIL="mme@example.org", SEND_EMAILS=True,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class InboundNotificationTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.email = "lab@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()
        self.classification = make_classification(self.lab, self.user)

    def _serve_query(self, patient=None) -> MMEInboundQuery:
        patient = patient or QUERY_PATIENT
        inbound_query = MMEInboundQuery.objects.create(
            peer_node_id="genematcher", request_json={"patient": patient}, num_results=1)
        record_inbound_matches(inbound_query, [_Match(self.classification)], patient)
        return inbound_query

    def test_records_which_of_our_classifications_we_handed_over(self):
        inbound_query = self._serve_query()
        match = MMEInboundMatch.objects.get(inbound_query=inbound_query)
        self.assertEqual(match.classification, self.classification)
        self.assertEqual(match.remote_patient_id, "remote-1")
        self.assertEqual(match.query_patient_json, QUERY_PATIENT)
        self.assertIsNone(match.notified)

    def test_first_match_messages_curator_emails_lab_and_stamps_notified(self):
        inbound_query = self._serve_query()
        with patch("mme.notifications.EmailLog.send_mail") as mock_send:
            notify_depositors(inbound_query.pk)

        msg = Message.objects.filter(recipient=self.user).first()
        self.assertIsNotNone(msg)
        self.assertIn("genematcher", msg.subject)
        self.assertIn("Dr Peer", msg.body)              # their contact
        self.assertIn("/mme/inbound_match/", msg.body)  # a direct link to the rest

        _args, kwargs = mock_send.call_args
        self.assertEqual(kwargs["recipient_list"], ["lab@lab.org"])

        self.assertIsNotNone(MMEInboundMatch.objects.get(inbound_query=inbound_query).notified)

    def test_repeat_query_for_the_same_remote_patient_records_but_stays_silent(self):
        with patch("mme.notifications.EmailLog.send_mail"):
            notify_depositors(self._serve_query().pk)
            second_query = self._serve_query()
            notify_depositors(second_query.pk)

        self.assertEqual(MMEInboundMatch.objects.count(), 2)     # audit trail keeps both
        self.assertEqual(Message.objects.filter(recipient=self.user).count(), 1)
        self.assertIsNone(MMEInboundMatch.objects.get(inbound_query=second_query).notified)

    def test_a_different_remote_patient_is_a_new_case_and_does_notify(self):
        with patch("mme.notifications.EmailLog.send_mail"):
            notify_depositors(self._serve_query().pk)
            notify_depositors(self._serve_query(dict(QUERY_PATIENT, id="remote-2")).pk)

        self.assertEqual(Message.objects.filter(recipient=self.user).count(), 2)

    def test_notified_only_stamped_when_a_channel_accepted_it(self):
        inbound_query = self._serve_query()
        self.lab.email = ""
        self.lab.group_name = None
        self.lab.save()
        with patch("mme.notifications.Message.objects.create", side_effect=RuntimeError("db down")), \
                patch("mme.notifications.AdminNotificationBuilder"):
            notify_depositors(inbound_query.pk)

        self.assertIsNone(MMEInboundMatch.objects.get(inbound_query=inbound_query).notified)

    def test_unattributed_query_still_notifies(self):
        # Rows written before per-peer tokens existed carry no peer_node_id.
        inbound_query = MMEInboundQuery.objects.create(request_json={}, num_results=1)
        record_inbound_matches(inbound_query, [_Match(self.classification)], QUERY_PATIENT)
        with patch("mme.notifications.EmailLog.send_mail"):
            notify_depositors(inbound_query.pk)
        self.assertIn("an MME peer", Message.objects.filter(recipient=self.user).first().subject)
