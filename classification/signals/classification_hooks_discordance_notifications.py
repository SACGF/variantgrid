from collections import defaultdict
from typing import Optional

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
        labs_to_notify = outstanding_notifications.values('lab').distinct()
        outstanding_notifications = outstanding_notifications.select_for_update()
        lab_notifications = defaultdict(list)

        for lab_info in labs_to_notify:
            lab_id = lab_info['lab']
            lab_obj = Lab.objects.get(pk=lab_id)

            discordance_entries = outstanding_notifications.filter(lab=lab_obj)

            for entry in discordance_entries:
                all_labs = entry.discordance_report.involved_labs.keys()
                # all_lab_names = ", ".join(lab.name for lab in all_labs)
                report_url = get_url_from_view_path(
                    reverse('discordance_report', kwargs={'discordance_report_id': entry.discordance_report.id}),
                )
                clin_sig_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

                labs_str = "\n".join(f"* {lab}" for lab in sorted(all_labs))
                nb = NotificationBuilder("Discordance notifications")
                nb.add_markdown(
                    f":fire_engine: :email: Sending Discordance Report <{report_url}|(DR_{entry.discordance_report.pk})> notification to {entry.lab}")
                # nb.add_field("Labs", labs_str)
                # nb.add_field("Trigger for notification", entry.cause)
                nb.send()

                notification = LabNotificationBuilder(lab=entry.lab,
                                                      message=f"Discordance Update (DR_{entry.discordance_report.id})")

                user_perspective = LabPickerData.for_lab(lab=entry.lab)
                report_summary = DiscordanceReportRowData(discordance_report=entry.discordance_report,
                                                          perspective=user_perspective)
                if resolution_text := entry.discordance_report.resolution_text:
                    notification.add_markdown(f"The below overlap is now marked as *{resolution_text}*")
                # notification.add_markdown(f"The labs {all_lab_names} are involved in the following discordance:")

                notification.add_field(label="Discordance Detected On", value=report_summary.date_detected_str)

                c_hgvs_str = "\n".join((str(chgvs) for chgvs in report_summary.c_hgvses))
                notification.add_field(label="c.HGVS", value=c_hgvs_str)

                sig_lab: ClassificationLabSummary
                for sig_lab in report_summary.lab_significances:
                    count_text = ""
                    if sig_lab.count > 1:
                        count_text = f" x {sig_lab.count}"
                    pending_text = ""
                    if sig_lab.pending:
                        pending_text = " (PENDING)"

                    if sig_lab.changed:
                        notification.add_field(f"{sig_lab.lab}{count_text}",
                                               f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)} -> {clin_sig_key.pretty_value(sig_lab.clinical_significance_to)}{pending_text}")
                    else:
                        notification.add_field(f"{sig_lab.lab}{count_text}",
                                               f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)}")
                # don't want to include notes in email as the text might be too sensitive
                notification.add_markdown(f"Full details of the overlap can be seen here : <{report_url}>")
                lab_notifications[lab_obj].append(notification)

    for lab, notifications in lab_notifications.items():
        report_ids = [f"DR_{notification.message.split('_')[1][:-1]}" for notification in notifications]
        num_report_ids = len(report_ids)  # Count of report IDs

        if num_report_ids > 6:
            combined_message = f"{num_report_ids} Discordance Updates"
        else:
            combined_message = f"Discordance Update for ({', '.join(report_ids)})"

        combined_notification = LabNotificationBuilder(
            lab=lab,
            message=combined_message
        )

        for notification in notifications:
            notification.sent = True
            combined_notification.merge(notification)
        combined_notification.send()
        current_date = timezone.now()
        # Updating the notification_sent_date for the DiscordanceNotification records
        DiscordanceNotification.objects.filter(lab=lab, notification_sent_date=None).update(
            notification_sent_date=current_date)
