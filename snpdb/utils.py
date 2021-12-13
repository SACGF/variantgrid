from enum import Enum
from typing import Optional, List

from django.conf import settings
from django.contrib.auth.models import Group, User

from email_manager.models import EmailLog
from library.log_utils import NotificationBuilder, send_notification
from snpdb.models import Lab, UserSettings


class LabNotificationBuilder(NotificationBuilder):

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
        if slack_web_hook := self.lab.slack_webhook:
            # send to lab's external slack
            send_notification(message=self.message, blocks=self.as_slack(), slack_webhook_url=slack_web_hook)

        if self.notification_type == LabNotificationBuilder.NotificationType.SLACK_ONLY:
            return

        recipient_list: List[str] = list()
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
