from typing import Optional

from django.conf import settings
from django.core.mail import send_mail, get_connection, EmailMultiAlternatives
from django.db import models
from django.urls import reverse
from model_utils.models import TimeStampedModel

from library.preview_request import PreviewData, PreviewKeyValue, PreviewModelMixin


class EmailLog(PreviewModelMixin, TimeStampedModel):
    subject = models.TextField()
    html = models.TextField()
    text = models.TextField()
    from_email = models.TextField(null=True, blank=True)
    recipient_list = models.TextField()
    probably_sent = models.BooleanField()
    single_email = models.BooleanField()

    def get_absolute_url(self):
        return reverse('email_detail', kwargs={"email_id": self.pk})

    @staticmethod
    def send_mail(subject: str,
                  html: str,
                  text: str,
                  from_email: Optional[str],
                  recipient_list: list[str],
                  allow_users_to_see_others: bool = False) -> bool:

        if html:
            html = html.strip()
        if text:
            text = text.strip()
        maybe_sent = False
        if from_email and recipient_list and settings.SEND_EMAILS:
            if allow_users_to_see_others:
                maybe_sent = send_mail(
                    subject=subject,
                    html_message=html,
                    message=text,
                    from_email=from_email,
                    recipient_list=recipient_list,
                    fail_silently=True
                )
            else:
                connection = get_connection(
                    fail_silently=True
                )
                email_messages = []
                for recipient in recipient_list:
                    mail = EmailMultiAlternatives(subject=subject, body=text, from_email=from_email,
                                                  to=[recipient], connection=connection)
                    if html:
                        mail.attach_alternative(html, 'text/html')
                    email_messages.append(mail)
                maybe_sent = connection.send_messages(email_messages)

        EmailLog.objects.create(
            subject=subject,
            html=html,
            text=text,
            from_email=from_email,
            recipient_list='; '.join(recipient_list) if recipient_list else '',
            probably_sent=maybe_sent,
            single_email=allow_users_to_see_others
        )
        return maybe_sent

    @classmethod
    def preview_category(cls) -> str:
        return "Emails"

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-envelope"

    @property
    def preview(self) -> 'PreviewData':
        title = self.subject
        extras = [
            PreviewKeyValue(key="To", value=self.recipient_list)
        ]
        return PreviewData.for_object(
            obj=self,
            category="Email",
            title=title,
            icon=self.preview_icon(),
            internal_url=self.get_absolute_url(),
            summary_extra=extras

        )
