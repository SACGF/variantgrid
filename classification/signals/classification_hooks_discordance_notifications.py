from typing import Optional, List, Set
from django.conf import settings
from django.db import transaction
from django.db.models import QuerySet
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone
from classification.enums import SpecialEKeys
from classification.models import DiscordanceReport, discordance_change_signal, EvidenceKeyMap, \
    ClassificationLabSummary
from classification.models.clinical_context_models import DiscordanceNotification, ClinicalContextChangeData, \
    ClinicalContextRecalcTrigger
from classification.models.discordance_models_utils import DiscordanceReportRowData
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab
from snpdb.utils import LabNotificationBuilder

"""
Responsible for emailing/slacking users when a discordance is detected.
In future we'd like this to occur in bulk (as a single import can create or solve many discordances).
Discordance is actually detected by classification_hooks_discordance_status.py
"""


@receiver(discordance_change_signal, sender=DiscordanceReport)
def notify_discordance_change(discordance_report: DiscordanceReport, clinical_context_change_data: ClinicalContextChangeData, **kwargs):
    if settings.DISCORDANCE_ENABLED and clinical_context_change_data.notify_worthy:
        prepare_discordance_notification(discordance_report=discordance_report, cause=clinical_context_change_data.cause_text)
        if clinical_context_change_data.cause_code != ClinicalContextRecalcTrigger.DELAYED:
            send_prepared_discordance_notifications()


def prepare_discordance_notification(discordance_report: DiscordanceReport, cause: Optional[str] = None):
    all_labs = discordance_report.involved_labs.keys()
    for lab in all_labs:
        DiscordanceNotification.objects.create(
            lab=lab,
            discordance_report=discordance_report,
            cause=cause
        )


def send_prepared_discordance_notifications(outstanding_notifications: Optional[QuerySet[DiscordanceNotification]] = None):
    if outstanding_notifications is None:
        outstanding_notifications = DiscordanceNotification.objects.filter(notification_sent_date__isnull=True).order_by('pk')

    with transaction.atomic():
        labs_to_notify = outstanding_notifications.order_by('lab').values_list('lab', flat=True).distinct()
        outstanding_notifications = outstanding_notifications.select_for_update()

        for lab_id in labs_to_notify:
            lab = Lab.objects.get(pk=lab_id)
            dr_ids: Set[int] = set()

            outstanding_lab_discordance_notifications: List[DiscordanceNotification] = list(outstanding_notifications.filter(lab=lab))
            combined_notifications: List[LabNotificationBuilder] = []

            for outstanding_notification in outstanding_lab_discordance_notifications:
                dr_id = outstanding_notification.discordance_report.pk

                if dr_id in dr_ids:
                    # don't notify about the same discordance report twice if somehow we have several pending notifications
                    # for the same discordance
                    continue

                dr_ids.add(dr_id)

                report_url = get_url_from_view_path(
                    reverse('discordance_report', kwargs={'discordance_report_id': dr_id}),
                )
                clin_sig_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

                lab_notification = LabNotificationBuilder(
                    lab=lab,
                    message=f"Discordance Update (DR_{dr_id})"
                )

                user_perspective = LabPickerData.for_lab(lab=outstanding_notification.lab)
                report_summary = DiscordanceReportRowData(discordance_report=outstanding_notification.discordance_report,
                                                          perspective=user_perspective)
                if resolution_text := outstanding_notification.discordance_report.resolution_text:
                    lab_notification.add_markdown(f"The below overlap is now marked as *{resolution_text}*")
                # notification.add_markdown(f"The labs {all_lab_names} are involved in the following discordance:")

                lab_notification.add_field(label="Discordance Detected On", value=report_summary.date_detected_str)

                c_hgvs_str = "\n".join((str(chgvs) for chgvs in report_summary.c_hgvses))
                lab_notification.add_field(label="c.HGVS", value=c_hgvs_str)

                sig_lab: ClassificationLabSummary
                for sig_lab in report_summary.lab_significances:
                    count_text = ""
                    if sig_lab.count > 1:
                        count_text = f" x {sig_lab.count}"
                    pending_text = ""
                    if sig_lab.pending:
                        pending_text = " (PENDING)"

                    if sig_lab.changed:
                        lab_notification.add_field(
                            f"{sig_lab.lab}{count_text}",
                            f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)} -> {clin_sig_key.pretty_value(sig_lab.clinical_significance_to)}{pending_text}")
                    else:
                        lab_notification.add_field(
                            f"{sig_lab.lab}{count_text}",
                            f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)}")

                # don't want to include notes in email as the text might be too sensitive
                lab_notification.add_markdown(f"Full details of the overlap can be seen here : <{report_url}|(DR_{dr_id})>")
                combined_notifications.append(lab_notification)

            # end loop
            # we now have combined_notifications to send
            dr_ids = list(sorted(dr_ids))
            current_date = timezone.now()
            subject: str

            # Admin Notification
            admin_notification = NotificationBuilder("Discordance notifications")
            admin_notification.add_markdown(":email: Sending Discordance Notifications")
            for dr_id in dr_ids:
                admin_notification.add_field("Lab", str(lab))
                link_text = []
                for dr_id in dr_ids:
                    report_url = get_url_from_view_path(reverse('discordance_report', kwargs={'discordance_report_id': dr_id}))
                    link_text.append(f"<{report_url}|(DR_{dr_id})>")

            admin_notification.add_field("Discordances", ", ".join(link_text))
            admin_notification.send()

            dr_count = len(dr_ids)
            if dr_count > 6:
                subject = f"Discordance Update for {dr_count} Discordances"
            else:
                subject = ", ".join([f"DR_{dr_id}" for dr_id in dr_ids])

            LabNotificationBuilder(
                lab=lab,
                message=subject
            ).merge(*combined_notifications).send()

            for outstanding_notification in outstanding_lab_discordance_notifications:
                outstanding_notification.notification_sent_date = current_date
            DiscordanceNotification.objects.bulk_update(outstanding_lab_discordance_notifications, fields=['notification_sent_date'])
