from django.conf import settings
from django.dispatch import receiver
from django.urls import reverse

from classification.enums import SpecialEKeys
from classification.models import DiscordanceReport, discordance_change_signal, EvidenceKeyMap
from classification.models.classification_groups import ClassificationGroups
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder
from snpdb.utils import LabNotificationBuilder


@receiver(discordance_change_signal, sender=DiscordanceReport)
def notify_discordance_change(discordance_report: DiscordanceReport, **kwargs):
    if settings.DISCORDANCE_ENABLED:
        send_discordance_notification(discordance_report=discordance_report)


def send_discordance_notification(discordance_report: DiscordanceReport):
    all_labs = discordance_report.all_actively_involved_labs
    all_lab_names = ", ".join(lab.name for lab in all_labs)
    groups = ClassificationGroups(discordance_report.all_classification_modifications)
    report_url = get_url_from_view_path(
        reverse('discordance_report', kwargs={'discordance_report_id': discordance_report.id}),
    )
    clin_sig_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    for lab in all_labs:
        notification = LabNotificationBuilder(lab=lab, message=f"Discordance Update (DR_{discordance_report.id})")
        if resolution_text := discordance_report.resolution_text:
            notification.add_markdown(f"The below overlap is now marked as *{resolution_text}*")
        notification.add_markdown(f"The labs {all_lab_names} are involved in the following discordance:")
        listing = ""
        for group in groups:
            listing += f"- {group.lab} `{group.most_recent.c_parts}` {clin_sig_key.pretty_value(group.clinical_significance)}"
            if group.count() > 1:
                listing += f" x {group.count()} records"
            if group.is_withdrawn:
                listing += " *WITHDRAWN*"
            listing += "\n"
        notification.add_markdown(listing)
        notification.add_markdown(f"Full details of the discordance can be seen here : <{report_url}>")
        notification.send()

    labs_notified = ", ".join(sorted([lab.name for lab in all_labs]))
    NotificationBuilder(message=f"Discordance Notification <{report_url}> sent to {labs_notified}")\
        .add_markdown(f":fire_engine: :email: Discordance Notification <{report_url}> sent to {labs_notified}").send()
