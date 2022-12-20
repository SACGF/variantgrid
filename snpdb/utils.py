from enum import Enum
from typing import Optional, List

from django.conf import settings
from django.contrib.auth.models import Group, User

from email_manager.models import EmailLog
from library.log_utils import NotificationBuilder, send_notification
from library.utils import empty_to_none
from snpdb.models import Lab, UserSettings, Tag, TagColorsCollection


class LabNotificationBuilder(NotificationBuilder):
    """
    For notifying members of the lab about an important event.
    Will email individual users with email enabled, or just the lab_email if one has been provided
    """

    class NotificationType(Enum):
        DISCORDANCE = "Discordance"
        GENERAL_UPDATE = "General Update"
        SLACK_ONLY = "Slack Only"

    def __init__(self, lab: Lab, message: str, notification_type: NotificationType = NotificationType.DISCORDANCE):
        if not isinstance(lab, Lab):
            raise ValueError(f"Expected lab, got {lab}")
        self.lab = lab
        self.notification_type = notification_type
        super().__init__(message=message)

    def as_html(self):
        base_html = super().as_html()
        return f"<p>Hello {settings.SITE_NAME} User,</p>{base_html}<p>Thanks,<br/>The {settings.SITE_NAME} Team.</p>"

    def as_text(self):
        base_text = super().as_text()
        return f"Hello {settings.SITE_NAME} User,\n\n{base_text}\n\nThanks,\nThe {settings.SITE_NAME} Team."

    def send(self):
        self.sent = True
        if slack_web_hook := empty_to_none(self.lab.slack_webhook):
            # send to lab's external slack
            send_notification(message=self.message, blocks=self.as_slack(), slack_webhook_url=slack_web_hook)

        if self.notification_type == LabNotificationBuilder.NotificationType.SLACK_ONLY:
            return

        recipient_list: List[str] = []
        if lab_email := self.lab.email:
            recipient_list.append(lab_email)
        else:
            lab_users = Group.objects.get(name=self.lab.group_name).user_set.filter(is_active=True, email__isnull=False)
            user: User
            for user in lab_users:
                passes_check = False
                if self.notification_type == LabNotificationBuilder.NotificationType.DISCORDANCE:
                    passes_check = UserSettings.get_for_user(user).email_discordance_updates
                elif self.notification_type == LabNotificationBuilder.NotificationType.GENERAL_UPDATE:
                    passes_check = UserSettings.get_for_user(user).email_weekly_updates
                if passes_check:
                    recipient_list.append(user.email)

        if recipient_list:
            EmailLog.send_mail(subject=self.message,
                               html=self.as_html(),
                               text=self.as_text(),
                               from_email=settings.DISCORDANCE_EMAIL,
                               recipient_list=recipient_list,
                               allow_users_to_see_others=True)

    @property
    def webhook_url(self) -> Optional[str]:
        return self.lab.slack_webhook


def get_all_tags_and_user_colors(user, tag_colors_collection=None):
    """ Returns Hash of { tag_name : color }
        with color being None if not set for user """

    if tag_colors_collection is None:
        user_settings = UserSettings.get_for_user(user)
        tag_colors_collection = user_settings.tag_colors

    user_colors_by_tag = {}
    if tag_colors_collection:
        user_colors_by_tag = tag_colors_collection.get_user_colors_by_tag()

    user_tag_colors = {}
    for tag in Tag.objects.all().order_by("pk"):
        user_tag_colors[tag] = user_colors_by_tag.get(tag.id)
    return user_tag_colors


def get_tag_styles_and_colors(user, tag_colors_collection: TagColorsCollection = None):
    user_tag_styles = []
    user_tag_colors = {}

    for tag, style in get_all_tags_and_user_colors(user, tag_colors_collection=tag_colors_collection).items():
        user_tag_styles.append((tag.id, style))
        rgb = None
        if style:
            rgb = style.get("background-color")
        user_tag_colors[tag.id] = rgb
    return user_tag_styles, user_tag_colors
