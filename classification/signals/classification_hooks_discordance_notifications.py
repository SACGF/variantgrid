from collections import defaultdict
from typing import Optional

from django.conf import settings
from django.db import transaction
from django.db.models import QuerySet
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone
from classification.models import DiscordanceReport, discordance_change_signal, \
    OverlapDiscordanceNotification, Overlap, OverlapContributionStatus, ClassificationImportRun
from classification.models.clinical_context_models import DiscordanceNotification, ClinicalContextChangeData, \
    ClinicalContextRecalcTrigger
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder
from snpdb.models import Lab
from snpdb.utils import LabNotificationBuilder

"""
Responsible for emailing/slacking users when a discordance is detected.
In future we'd like this to occur in bulk (as a single import can create or solve many discordances).
Discordance is actually detected by classification_hooks_discordance_status.py
"""


@receiver(discordance_change_signal, sender=DiscordanceReport)
def notify_discordance_change(discordance_report: DiscordanceReport, clinical_context_change_data: ClinicalContextChangeData, **kwargs):
    # if settings.DISCORDANCE_ENABLED and clinical_context_change_data.notify_worthy:
    #     prepare_discordance_notification(discordance_report=discordance_report, cause=clinical_context_change_data.cause_text)
    #     if clinical_context_change_data.cause_code != ClinicalContextRecalcTrigger.DELAYED:
    #         NotificationBuilder("Send Notifications Triggered")\
    #             .add_markdown(f"Detected discordance change outside of an import for DR_{discordance_report.pk} casued by *{clinical_context_change_data.cause_code.value}* - will send out notifications now.")\
    #             .send()
    #         send_prepared_discordance_notifications()
    pass


def prepare_discordance_notification(discordance_report: DiscordanceReport, cause: Optional[str] = None):
    pass
    # all_labs = discordance_report.involved_labs.keys()
    # for lab in all_labs:
    #     DiscordanceNotification.objects.create(
    #         lab=lab,
    #         discordance_report=discordance_report,
    #         cause=cause
    #     )


def _report_url_for_id(overlap: Overlap, lab: Lab):
    return get_url_from_view_path(overlap.get_absolute_url())


def send_prepared_discordance_notifications(outstanding_notifications: Optional[QuerySet[OverlapDiscordanceNotification]] = None):
    if ClassificationImportRun.ongoing_imports():
        # don't send notifications while an import is ongoing
        print("There's an ongoing import, wont send notifications")
        return

    if outstanding_notifications is None:
        outstanding_notifications = OverlapDiscordanceNotification.objects.filter(notification_sent_date__isnull=True).order_by('pk')

    outstanding_notifications = outstanding_notifications.select_related("overlap")
    with transaction.atomic():
        outstanding_notifications = outstanding_notifications.select_for_update()

        current_date = timezone.now()
        relevant_notifications: list[OverlapDiscordanceNotification] = []
        notifications_by_lab: dict[Lab, list[OverlapDiscordanceNotification]] = defaultdict(list)

        for notification in outstanding_notifications:
            if not notification.is_still_relevant:
                notification.delete()
            else:
                # send one notification per discordance so we don't try to send a single message too big for slack
                overall_admin_notification = NotificationBuilder("Discordance notification").add_markdown(
                    ":email: Sending Discordance Notification")
                notification.notification_sent_date = current_date
                relevant_notifications.append(notification)
                # TODO detect which labs were once a part but then withdrew
                # though currently not possible to tell if they were part of the allele from years back or from minutes ago and just withdrew
                labs: set[Lab] = set()
                for contribution in notification.overlap.contributions.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING):
                    if lab := contribution.lab:
                        labs.add(lab)
                        notifications_by_lab[lab].append(notification)

                sorted_lab_str = "\n".join(str(lab) for lab in sorted(labs))

                overlap = notification.overlap
                discordance_status_icon = ":no_good:" if notification.new_status.is_discordant else ":handshake:"
                overlap_description = f"{overlap.scope_description} {overlap.value_type_label}"
                overlap_change = f"{notification.old_status.label} -> {notification.new_status.label} {discordance_status_icon}"

                overall_admin_notification.add_field("Overlap", f"<{_report_url_for_id(notification.overlap, None)}|Overlap_{notification.overlap_id}> {overlap_description}")
                overall_admin_notification.add_field("Involved Labs", sorted_lab_str)
                overall_admin_notification.add_field("Change", overlap_change)
                overall_admin_notification.send()

        # FIXME make an admin notification

        for lab, notifications in notifications_by_lab.items():
            notifications_list = list(sorted(notifications))
            notification_count = len(notifications)
            if notification_count > 6:
                subject = f"Discordance Update for {notification_count} Discordances"
            else:
                subject = "Discordance Update for (" + ", ".join([f"OL_{notification.overlap_id}" for notification in notifications_list]) + ")"

            lab_notification = LabNotificationBuilder(lab=lab, message=subject)

            for index, notification in enumerate(notifications):
                if not index == 0:
                    lab_notification.add_divider()

                report_url = _report_url_for_id(notification.overlap, lab)
                # report_summary = DiscordanceReportRowData(discordance_report=dr, perspective=user_perspective)
                lab_notification.add_markdown(f"The below overlap is now marked as *{notification.new_status.label}*")
                # notification.add_markdown(f"The labs {all_lab_names} are involved in the following discordance:")

                lab_notification.add_field(label="Discordance Report", value=f"<{report_url}|OL_{notification.overlap_id, lab}>")
                # lab_notification.add_field(label="Discordance Detected On", value=report_summary.date_detected_str)
                # FIXME add many more details about the Overlap

            lab_notification.send()
        OverlapDiscordanceNotification.objects.bulk_update(relevant_notifications, fields=['notification_sent_date'])


                # FIXME get c.HGVSs back
                # c_hgvs_str = "\n".join((str(chgvs) for chgvs in report_summary.c_hgvses))
                # lab_notification.add_field(label="c.HGVS", value=c_hgvs_str)

                # FIXME have to go through Overlap to get individual rows, and maybe Audit Log to see what was the trigger
                # sig_lab: ClassificationLabSummary
                # for sig_lab in report_summary.lab_significances:
                #     count_text = ""
                #     if sig_lab.count > 1:
                #         count_text = f" x {sig_lab.count}"
                #     pending_text = ""
                #     if sig_lab.pending:
                #         pending_text = " (PENDING)"
                #
                #     if sig_lab.changed:
                #         lab_notification.add_field(
                #             f"{sig_lab.lab}{count_text}",
                #             f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)} -> {clin_sig_key.pretty_value(sig_lab.clinical_significance_to)}{pending_text}")
                #     else:
                #         lab_notification.add_field(
                #             f"{sig_lab.lab}{count_text}",
                #             f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)}")


# def send_prepared_discordance_notifications(outstanding_notifications: Optional[QuerySet[DiscordanceNotification]] = None):
#     if outstanding_notifications is None:
#         outstanding_notifications = DiscordanceNotification.objects.filter(notification_sent_date__isnull=True).order_by('pk')
#
#     with transaction.atomic():
#         labs_to_notify = outstanding_notifications.order_by('lab').values_list('lab', flat=True).distinct()
#         outstanding_notifications = outstanding_notifications.select_for_update()
#
#         for lab_id in labs_to_notify:
#             lab = Lab.objects.get(pk=lab_id)
#             user_perspective = LabPickerData.for_lab(lab=lab)
#
#             outstanding_lab_discordance_notifications: list[DiscordanceNotification] = list(outstanding_notifications.filter(lab=lab))
#             current_date = timezone.now()
#
#             unique_ids: set[int] = set()
#             for outstanding_notification in outstanding_lab_discordance_notifications:
#                 outstanding_notification.notification_sent_date = current_date
#                 dr_id = outstanding_notification.discordance_report_id
#                 unique_ids.add(dr_id)
#
#             dr_ids: list[int] = list(sorted(unique_ids))
#             drs: list[DiscordanceReport] = DiscordanceReport.objects.filter(pk__in=dr_ids).order_by('pk')
#
#             dr_count = len(dr_ids)
#             if dr_count > 6:
#                 subject = f"Discordance Update for {dr_count} Discordances"
#             else:
#                 subject = "Discordance Update for (" + ", ".join([f"DR_{dr_id}" for dr_id in dr_ids]) + ")"
#
#             lab_notification = LabNotificationBuilder(lab=lab, message=subject)
#
#             def report_url_for_id(the_id):
#                 return get_url_from_view_path(
#                     reverse('discordance_report', kwargs={'discordance_report_id': the_id}),
#                 )
#
#             # Admin Notification
#             admin_notification = NotificationBuilder("Discordance notifications")\
#                 .add_markdown(":email: Sending Discordance Notifications")\
#                 .add_field("Lab", str(lab))\
#                 .add_field("Discordance IDs", ", ".join(f"<{report_url_for_id(dr_id)}|DR_{dr_id}>" for dr_id in dr_ids))
#
#             is_first = True
#             for dr in drs:
#                 if is_first:
#                     is_first = False
#                 else:
#                     lab_notification.add_divider()
#
#                 dr_id = dr.pk
#                 report_url = report_url_for_id(dr_id)
#                 clin_sig_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
#
#                 report_summary = DiscordanceReportRowData(discordance_report=dr,
#                                                           perspective=user_perspective)
#                 if resolution_text := dr.resolution_text:
#                     lab_notification.add_markdown(f"The below overlap is now marked as *{resolution_text}*")
#                 # notification.add_markdown(f"The labs {all_lab_names} are involved in the following discordance:")
#
#                 lab_notification.add_field(label="Discordance Report", value=f"<{report_url}|DR_{dr_id}>")
#
#                 lab_notification.add_field(label="Discordance Detected On", value=report_summary.date_detected_str)
#
#                 c_hgvs_str = "\n".join((str(chgvs) for chgvs in report_summary.c_hgvses))
#                 lab_notification.add_field(label="c.HGVS", value=c_hgvs_str)
#
#                 sig_lab: ClassificationLabSummary
#                 for sig_lab in report_summary.lab_significances:
#                     count_text = ""
#                     if sig_lab.count > 1:
#                         count_text = f" x {sig_lab.count}"
#                     pending_text = ""
#                     if sig_lab.pending:
#                         pending_text = " (PENDING)"
#
#                     if sig_lab.changed:
#                         lab_notification.add_field(
#                             f"{sig_lab.lab}{count_text}",
#                             f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)} -> {clin_sig_key.pretty_value(sig_lab.clinical_significance_to)}{pending_text}")
#                     else:
#                         lab_notification.add_field(
#                             f"{sig_lab.lab}{count_text}",
#                             f"{clin_sig_key.pretty_value(sig_lab.clinical_significance_from)}")
#
#             admin_notification.send()
#             lab_notification.send()
#             DiscordanceNotification.objects.bulk_update(outstanding_lab_discordance_notifications, fields=['notification_sent_date'])
