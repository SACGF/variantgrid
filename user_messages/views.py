from django.contrib import messages
from django.contrib.auth import get_user_model
from django.http import Http404, HttpResponseRedirect
from django.shortcuts import render, get_object_or_404
from django.urls import reverse
from django.utils import timezone
from django.utils.translation import gettext as _

from user_messages.forms import ComposeForm
from user_messages.models import Message

User = get_user_model()


def inbox(request):
    """ Displays a list of received messages for the current user. """
    message_list = Message.objects.inbox_for(request.user)
    return render(request, 'user_messages/inbox.html', {'message_list': message_list})


def outbox(request):
    """ Displays a list of sent messages by the current user. """
    message_list = Message.objects.outbox_for(request.user)
    return render(request, 'user_messages/outbox.html', {'message_list': message_list})


def trash(request):
    """ Displays a list of deleted messages. """
    message_list = Message.objects.trash_for(request.user)
    return render(request, 'user_messages/trash.html', {'message_list': message_list})


def compose(request, recipient=None):
    """ Displays and handles the form to compose a new message. """
    if request.method == "POST":
        form = ComposeForm(request.POST)
        if form.is_valid():
            form.save(sender=request.user)
            messages.info(request, _("Message successfully sent."))
            return HttpResponseRedirect(request.GET.get('next') or reverse('messages_inbox'))
    else:
        form = ComposeForm(initial={"subject": request.GET.get("subject", "")})
        if recipient is not None:
            form.fields['recipient'].initial = list(User.objects.filter(**{
                f'{User.USERNAME_FIELD}__in': [r.strip() for r in recipient.split('+')]
            }))
    return render(request, 'user_messages/compose.html', {'form': form})


def delete(request, message_id):
    """
    Marks a message as deleted by sender or recipient. The message is not really removed from the
    database, because two users must delete a message before it's safe to remove it completely.
    """
    user = request.user
    now = timezone.now()
    message = get_object_or_404(Message, id=message_id)
    deleted = False
    if message.sender == user:
        message.sender_deleted_at = now
        deleted = True
    if message.recipient == user:
        message.recipient_deleted_at = now
        deleted = True
    if deleted:
        message.save()
        messages.info(request, _("Message successfully deleted."))
        return HttpResponseRedirect(request.GET.get('next') or reverse('messages_inbox'))
    raise Http404


def undelete(request, message_id):
    """ Recovers a message from trash. """
    user = request.user
    message = get_object_or_404(Message, id=message_id)
    undeleted = False
    if message.sender == user:
        message.sender_deleted_at = None
        undeleted = True
    if message.recipient == user:
        message.recipient_deleted_at = None
        undeleted = True
    if undeleted:
        message.save()
        messages.info(request, _("Message successfully recovered."))
        return HttpResponseRedirect(request.GET.get('next') or reverse('messages_inbox'))
    raise Http404


def view(request, message_id):
    """
    Shows a single message. The user is only allowed to see the message if they are either the
    sender or the recipient, otherwise a 404 is raised. If the user is the recipient and the
    message is unread, ``read_at`` is set to now.
    """
    user = request.user
    message = get_object_or_404(Message, id=message_id)
    if (message.sender != user) and (message.recipient != user):
        raise Http404
    if message.read_at is None and message.recipient == user:
        message.read_at = timezone.now()
        message.save()
    return render(request, 'user_messages/view.html', {'message': message})
