from django.conf import settings
from django.db import models
from django.urls import reverse
from django.utils import timezone
from django.utils.safestring import SafeString
from django.utils.translation import gettext_lazy as _
from markdown import markdown

from library.utils import sanitize_html


class MessageManager(models.Manager):

    def inbox_for(self, user):
        """ Messages received by the given user that are not marked as deleted. """
        return self.filter(recipient=user, recipient_deleted_at__isnull=True)

    def outbox_for(self, user):
        """ Messages sent by the given user that are not marked as deleted. """
        return self.filter(sender=user, sender_deleted_at__isnull=True)

    def trash_for(self, user):
        """ Messages sent or received by the given user that are marked as deleted. """
        return self.filter(
            recipient=user, recipient_deleted_at__isnull=False,
        ) | self.filter(
            sender=user, sender_deleted_at__isnull=False,
        )


class Message(models.Model):
    """
    A private message from user to user. The ``body`` is authored/stored as Markdown and rendered
    safely for display via :attr:`body_html`.
    """
    subject = models.CharField(_("Subject"), max_length=140)
    body = models.TextField(_("Body"))
    sender = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='sent_messages', verbose_name=_("Sender"), on_delete=models.PROTECT)
    recipient = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='received_messages', null=True, blank=True, verbose_name=_("Recipient"), on_delete=models.SET_NULL)
    parent_msg = models.ForeignKey('self', related_name='next_messages', null=True, blank=True, verbose_name=_("Parent message"), on_delete=models.SET_NULL)
    sent_at = models.DateTimeField(_("sent at"), null=True, blank=True)
    read_at = models.DateTimeField(_("read at"), null=True, blank=True)
    replied_at = models.DateTimeField(_("replied at"), null=True, blank=True)
    sender_deleted_at = models.DateTimeField(_("Sender deleted at"), null=True, blank=True)
    recipient_deleted_at = models.DateTimeField(_("Recipient deleted at"), null=True, blank=True)

    objects = MessageManager()

    class Meta:
        ordering = ['-sent_at']
        verbose_name = _("Message")
        verbose_name_plural = _("Messages")

    def __str__(self):
        return self.subject

    @property
    def body_html(self) -> SafeString:
        """ Render the Markdown body to sanitized HTML that is safe to display unescaped. """
        return sanitize_html(markdown(self.body or ""))

    def new(self):
        """ Whether the recipient has read the message or not. """
        return self.read_at is None

    def replied(self):
        """ Whether the recipient has written a reply to this message. """
        return self.replied_at is not None

    def get_absolute_url(self):
        return reverse('messages_detail', args=[self.id])

    def save(self, **kwargs):
        if not self.id:
            self.sent_at = timezone.now()
        super().save(**kwargs)


def inbox_count_for(user):
    """ Number of unread messages for the given user (does not mark them seen). """
    return Message.objects.filter(recipient=user, read_at__isnull=True, recipient_deleted_at__isnull=True).count()
