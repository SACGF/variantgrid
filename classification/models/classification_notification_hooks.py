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
    all_labs = discordance_report.all_actively_involved_labs()
    all_lab_names = ", ".join(lab.name for lab in all_labs)
    groups = ClassificationGroups(discordance_report.all_classification_modifications())
    report_url = get_url_from_view_path(
        reverse('discordance_report', kwargs={'report_id': discordance_report.id}),
    )
    clin_sig_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    for lab in all_labs:
        notification = LabNotificationBuilder(lab=lab, message=f"Discordance Update ({discordance_report.id})")
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

        NotificationBuilder(message=f"Discordance notification re Discordance Report <{report_url}> sent to {lab.name}", emoji=":email:").send()


"""
@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):  # pylint: disable=unused-argument
    if settings.DISCORDANCE_ENABLED:
        resolution = flag_comment.resolution
        if resolution:
            flag: Flag = flag_comment.flag
            flag_type = flag.flag_type
            if flag_type == flag_types.classification_flag_types.discordant:
                vc = Classification.objects.filter(flag_collection=flag.collection.id).select_related('lab').get()
                email_discordance_for_classification(vc)


def email_discordance_for_report(report: DiscordanceReport):
    for vcm in report.all_classification_modifications():
        vc = vcm.classification
        email_discordance_for_classification(vc)


def email_discordance_for_classification(vc: Classification) -> bool:
    discordant_flag = vc.flag_collection.get_flag_of_type(flag_type=flag_types.classification_flag_types.discordant, open_only=False)
    if discordant_flag:
        is_open = discordant_flag.resolution.status == FlagStatus.OPEN
        lab = vc.lab

        subject = f"Discordance detected ({vc.id})"
        if not is_open:
            subject = f"re: {subject}"

        context = {
            "lab": lab,
            "entered": is_open,
            "review_link": get_url_from_view_path(
                reverse('view_classification', kwargs={'record_id': vc.id})),
        }

        html = render_to_string('classification/emails/discordance_email.html', context)
        text = render_to_string('classification/emails/discordance_email.txt', context)

        primary_users_to_notify: List[str] = list()

        # find users for the lab and super users
        lab_users = Group.objects.get(name=lab.group_name).user_set.filter(is_active=True)
        user: User
        for user in lab_users:
            if UserSettings.get_for_user(user).email_discordance_updates:
                primary_users_to_notify.append(user.email)

        sent = EmailLog.send_mail(subject=subject,
                                  html=html,
                                  text=text,
                                  from_email=settings.DISCORDANCE_EMAIL,
                                  recipient_list=primary_users_to_notify,
                                  allow_users_to_see_others=True)

        report_event('Discordance email sent', extra_data={"classification_id": vc.id, "lab": lab.name, "open": is_open, "recipient_count": len(primary_users_to_notify)})
        return sent

    return False
"""