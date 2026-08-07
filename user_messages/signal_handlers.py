from django.conf import settings

from email_manager.models.email_model import EmailLog


def email_new_message_handler(sender, instance, created, **kwargs):
    """
    On a newly created message, email the recipient via the logged EmailLog infrastructure
    (honours ``settings.SEND_EMAILS``). Disable with ``settings.USER_MESSAGES_EMAIL_NOTIFY = False``.
    """
    if not created:
        return
    if not getattr(settings, 'USER_MESSAGES_EMAIL_NOTIFY', True):
        return

    recipient = instance.recipient
    if not (recipient and recipient.email):
        return

    subject = f"New Message: {instance.subject}"
    html = f"<p>You have received a new message from {instance.sender}:</p>{instance.body_html}"
    EmailLog.send_mail(subject=subject,
                       html=html,
                       text=instance.body,
                       from_email=settings.DEFAULT_FROM_EMAIL,
                       recipient_list=[recipient.email])
