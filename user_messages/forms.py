from django import forms
from django.utils.translation import gettext_lazy as _

from user_messages.fields import CommaSeparatedUserField
from user_messages.models import Message


class ComposeForm(forms.Form):
    """ A simple default form for private messages. """
    recipient = CommaSeparatedUserField(label=_("Recipient"))
    subject = forms.CharField(label=_("Subject"), max_length=140)
    body = forms.CharField(
        label=_("Body"),
        help_text=_("Markdown is supported."),
        widget=forms.Textarea(attrs={'rows': '12', 'cols': '55'}))

    def save(self, sender):
        recipients = self.cleaned_data['recipient']
        subject = self.cleaned_data['subject']
        body = self.cleaned_data['body']
        message_list = []
        for r in recipients:
            msg = Message.objects.create(sender=sender, recipient=r, subject=subject, body=body)
            message_list.append(msg)
        return message_list
