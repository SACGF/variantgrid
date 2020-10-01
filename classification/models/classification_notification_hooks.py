from typing import List

from django.conf import settings
from django.contrib.auth.models import User, Group
from django.dispatch import receiver
from django.template.loader import render_to_string
from django.urls import reverse

from email_manager.models import EmailLog
from flags.models import flag_comment_action, Flag, FlagResolution, FlagComment, FlagStatus
from library.django_utils import get_url_from_view_path
from library.log_utils import report_event
from snpdb.models import UserSettings
from classification.models import flag_types, Classification, DiscordanceReport


@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):
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
