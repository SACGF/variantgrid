""" `notify_match` is the recipient policy for both directions - the authoring curator AND
the lab, because a case can sit unmatched for years and the author may be long gone. """
from unittest.mock import MagicMock, patch

from django.contrib.auth.models import User
from django.core.exceptions import ImproperlyConfigured
from django.test import TestCase, override_settings

from classification.tests.models.test_utils import ClassificationTestUtils
from mme.apps import MMEConfig
from mme.client import submit
from mme.models import MMESubmission
from mme.notifications import lab_notification_emails, notify_match
from mme.tests.fakes import make_classification
from user_messages.models import Message

NODES = {"testnode": {"base_url": "https://node.test", "token": "secret-token", "api_version": "1.1"}}


@override_settings(MME_FROM_EMAIL="mme@example.org", SEND_EMAILS=True)
class LabNotificationEmailsTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def test_lab_email_wins_when_set(self):
        self.lab.email = "lab@lab.org"
        self.lab.save()
        self.assertEqual(lab_notification_emails(self.lab), ["lab@lab.org"])

    def test_falls_back_to_active_lab_members(self):
        # No lab inbox set: every active member hears, so a match still lands.
        self.assertEqual(sorted(lab_notification_emails(self.lab)), ["joe@joe.com", "joe@joe2.com"])

    def test_inactive_and_blank_email_members_excluded(self):
        User.objects.filter(username="joejoe2").update(is_active=False)
        User.objects.filter(username="joejoe").update(email="")
        self.assertEqual(lab_notification_emails(self.lab), [])

    def test_lab_with_no_group_and_no_email_yields_no_address(self):
        self.lab.group_name = None
        self.lab.email = ""
        self.lab.save()
        self.assertEqual(lab_notification_emails(self.lab), [])


@override_settings(MME_FROM_EMAIL="mme@example.org", SEND_EMAILS=True, MME_ENABLED=True,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class NotifyMatchTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.email = "lab@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()
        self.classification = make_classification(self.lab, self.user)

    def test_notifies_curator_and_lab(self):
        with patch("mme.notifications.EmailLog.send_mail") as mock_send:
            self.assertTrue(notify_match(self.classification, "subj", "body"))

        self.assertTrue(Message.objects.filter(recipient=self.user, subject="subj").exists())
        _args, kwargs = mock_send.call_args
        self.assertEqual(kwargs["recipient_list"], ["lab@lab.org"])
        self.assertEqual(kwargs["from_email"], "mme@example.org")

    def test_inactive_curator_drops_to_the_lab_alone(self):
        User.objects.filter(pk=self.user.pk).update(is_active=False)
        self.classification.refresh_from_db()
        with patch("mme.notifications.EmailLog.send_mail") as mock_send:
            self.assertTrue(notify_match(self.classification, "subj", "body"))

        self.assertFalse(Message.objects.filter(recipient=self.user).exists())
        mock_send.assert_called_once()

    def test_one_channel_failing_still_counts_as_sent(self):
        # In-app inbox down, lab email fine: the match still reaches somebody.
        with patch("mme.notifications.Message.objects.create", side_effect=RuntimeError("db down")), \
                patch("mme.notifications.EmailLog.send_mail") as mock_send:
            self.assertTrue(notify_match(self.classification, "subj", "body"))
        mock_send.assert_called_once()

    def test_every_channel_failing_raises_an_admin_alert(self):
        self.lab.email = ""
        self.lab.group_name = None
        self.lab.save()
        self.classification.refresh_from_db()
        with patch("mme.notifications.Message.objects.create", side_effect=RuntimeError("db down")), \
                patch("mme.notifications.AdminNotificationBuilder") as mock_nb:
            self.assertFalse(notify_match(self.classification, "subj", "body"))
        mock_nb.assert_called_once()
        mock_nb.return_value.send.assert_called_once()


@override_settings(MME_NODES=NODES, MME_ENABLED=True, MME_FROM_EMAIL="mme@example.org",
                   SEND_EMAILS=True, MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class OutboundNotificationTestCase(TestCase):
    """ The outbound path notifies through the same helper. """

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, self.user = ClassificationTestUtils.lab_and_user()
        lab.contact_email = "curator@lab.org"
        lab.email = "lab@lab.org"
        lab.mme_enabled = True
        lab.save()
        self.classification = make_classification(lab, self.user)
        self.submission = MMESubmission.objects.create(
            classification=self.classification, node_id="testnode", external_patient_id="vg:1")

    def test_results_notify_curator_and_lab(self):
        response = MagicMock()
        response.raise_for_status.return_value = None
        response.json.return_value = {"results": [
            {"score": {"patient": 0.9}, "patient": {"id": "remote-1", "contact": {}}}]}
        with patch("mme.client.requests.post", return_value=response), \
                patch("mme.notifications.EmailLog.send_mail") as mock_send:
            submit(self.submission)

        self.assertTrue(Message.objects.filter(recipient=self.user).exists())
        _args, kwargs = mock_send.call_args
        self.assertEqual(kwargs["recipient_list"], ["lab@lab.org"])


class StartupCheckTestCase(TestCase):
    """ Email sending is silently a no-op when from_email or SEND_EMAILS is unset. """

    def _ready(self):
        MMEConfig.ready(MMEConfig("mme", __import__("mme")))

    @override_settings(MME_ENABLED=True, MME_FROM_EMAIL=None, SEND_EMAILS=True)
    def test_enabled_without_from_email_rejected(self):
        with self.assertRaises(ImproperlyConfigured):
            self._ready()

    @override_settings(MME_ENABLED=True, MME_FROM_EMAIL="mme@example.org", SEND_EMAILS=False)
    def test_enabled_without_send_emails_rejected(self):
        with self.assertRaises(ImproperlyConfigured):
            self._ready()

    @override_settings(MME_ENABLED=True, MME_FROM_EMAIL="mme@example.org", SEND_EMAILS=True)
    def test_fully_configured_accepted(self):
        self._ready()

    @override_settings(MME_ENABLED=False, MME_FROM_EMAIL=None, SEND_EMAILS=False)
    def test_disabled_deployment_needs_nothing(self):
        self._ready()
