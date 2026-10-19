"""
OverlapServices.send_prepared_discordance_notifications: messages go out after the commit, each send isolated,
a failed lab message retried for that lab alone (#2279).
"""
from unittest.mock import MagicMock, patch

from django.test import TestCase

from classification.enums.classification_enums import (
    AlleleOriginBucket,
    OverlapStatus,
    ShareLevel,
    TestingContextBucket,
)
from classification.enums.overlaps_enums import (
    ClassificationResultValue,
    OverlapContributionStatus,
    OverlapEntrySourceTextChoices,
    OverlapSnapshot,
    OverlapType,
)
from classification.models.classification_grouping import (
    AlleleOriginGrouping,
    ClassificationGrouping,
)
from classification.models.overlaps_model import (
    Overlap,
    OverlapContribution,
    OverlapContributionNextStep,
    OverlapDiscordanceNotification,
)
from classification.services.overlaps_services import OverlapServices
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.models.models import Lab
from snpdb.models.models_variant import Allele


class DiscordanceNotificationSendTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab_a, _ = ClassificationTestUtils.lab_and_user()
        self.lab_b = Lab.objects.create(name="Labby B", organization=self.lab_a.organization, city="CityA",
                                        country=self.lab_a.country, group_name="instx/labby_b")
        self.notifications = [self._discordant_notification([self.lab_a, self.lab_b]) for _ in range(2)]

    def _discordant_notification(self, labs: list[Lab]) -> OverlapDiscordanceNotification:
        allele = Allele.objects.create()
        overlap = Overlap.objects.create(overlap_type=OverlapType.SINGLE_CONTEXT,
                                         value_type=ClassificationResultValue.ONC_PATH,
                                         allele=allele, testing_context_bucket=TestingContextBucket.GERMLINE)
        allele_origin_grouping = AlleleOriginGrouping.objects.create(allele=allele,
                                                                     allele_origin_bucket=AlleleOriginBucket.GERMLINE)
        for lab in labs:
            grouping = ClassificationGrouping.objects.create(allele_origin_grouping=allele_origin_grouping, lab=lab,
                                                             share_level=ShareLevel.ALL_USERS)
            contribution = OverlapContribution.objects.create(
                source=OverlapEntrySourceTextChoices.CLASSIFICATION, allele=allele, classification_grouping=grouping,
                value_type=ClassificationResultValue.ONC_PATH, contribution_status=OverlapContributionStatus.CONTRIBUTING,
                testing_context_bucket=TestingContextBucket.GERMLINE)
            OverlapContributionNextStep.objects.create(overlap=overlap, contribution=contribution)

        return OverlapDiscordanceNotification.objects.create(
            overlap=overlap,
            old_state=OverlapSnapshot(status=OverlapStatus.EXACT_AGREEMENT).to_dict(),
            new_state=OverlapSnapshot(status=OverlapStatus.MEDICALLY_SIGNIFICANT,
                                      lab_groups=[lab.group_name for lab in labs]).to_dict())

    def _send(self, failing_lab: Lab | None = None) -> tuple[list[MagicMock], dict[Lab, MagicMock]]:
        """ Runs the send with the messages mocked, returns the admin and lab messages that were sent """
        admin_messages = []
        lab_messages = {}

        def admin_message(notification, _labs):
            message = MagicMock(message=f"admin {notification.pk}")
            admin_messages.append(message)
            return message

        def lab_message(lab, _notifications):
            message = MagicMock(message=f"lab {lab}", lab=lab)
            if lab == failing_lab:
                message.send.side_effect = RuntimeError("Slack is down")
            lab_messages[lab] = message
            return message

        with patch.object(OverlapServices, "_admin_message", side_effect=admin_message), \
                patch.object(OverlapServices, "_lab_message", side_effect=lab_message):
            with self.captureOnCommitCallbacks() as callbacks:
                OverlapServices.send_prepared_discordance_notifications()
            self.assertFalse(any(message.send.called for message in admin_messages + list(lab_messages.values())))
            for callback in callbacks:
                callback()
        return admin_messages, lab_messages

    def test_failed_lab_send_is_retried_for_that_lab_only(self):
        with patch("classification.services.overlaps_services.report_exc_info") as report_exc_info:
            admin_messages, lab_messages = self._send(failing_lab=self.lab_a)
        self.assertEqual(len(admin_messages), 2)
        self.assertEqual(set(lab_messages), {self.lab_a, self.lab_b})
        for message in admin_messages + list(lab_messages.values()):
            message.send.assert_called_once()
        report_exc_info.assert_called_once()

        unsent = OverlapDiscordanceNotification.objects.filter(notification_sent_date__isnull=True)
        self.assertEqual({(n.overlap_id, n.lab) for n in unsent},
                         {(n.overlap_id, self.lab_a) for n in self.notifications})

        # the retry goes to lab A alone, without repeating the admin messages
        admin_messages, lab_messages = self._send()
        self.assertEqual(admin_messages, [])
        self.assertEqual(set(lab_messages), {self.lab_a})
        self.assertFalse(OverlapDiscordanceNotification.objects.filter(notification_sent_date__isnull=True).exists())

    def test_failed_retry_stays_unsent(self):
        with patch("classification.services.overlaps_services.report_exc_info"):
            self._send(failing_lab=self.lab_a)
            retry_pks = set(OverlapDiscordanceNotification.objects.filter(
                notification_sent_date__isnull=True).values_list("pk", flat=True))
            self._send(failing_lab=self.lab_a)
        unsent_pks = set(OverlapDiscordanceNotification.objects.filter(
            notification_sent_date__isnull=True).values_list("pk", flat=True))
        self.assertEqual(unsent_pks, retry_pks)
