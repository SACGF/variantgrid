from random import randint

from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver

from eventlog.models import Event, ViewEvent
from library.guardian_utils import bot_group
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckRecentActivity


@receiver(signal=health_check_signal)
def active_users_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    event_log_activity = Event.objects.filter(date__gte=health_request.since, date__lt=health_request.now).values_list('user', flat=True).distinct()
    view_event_activity = ViewEvent.objects.filter(created__gte=health_request.since, created__lt=health_request.now).values_list('user', flat=True).distinct()

    active_users = User.objects.filter(Q(pk__in=event_log_activity) | Q(pk__in=view_event_activity)).exclude(groups=bot_group()).order_by('username') \
        .values_list('username', flat=True)

    count = active_users.count()
    emoji: str
    if count == 0:
        emoji = ":ghost:"
    else:
        emojis = [":nerd_face:", ":thinking_face:", ":face_with_monocle:", ":face_with_cowboy_hat:"]
        emoji = emojis[randint(0, len(emojis) - 1)]

    return HealthCheckRecentActivity(
        emoji=emoji,
        name="Active Users",
        amount=count,
        extra=", ".join(list(active_users)),
        stand_alone=True  # always give active users its own line
    )
