from collections import Counter
from random import randint

from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver

from classification.models import DiscordanceReportTriage
from eventlog.models import Event, ViewEvent
from library.guardian_utils import bot_group
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckRecentActivity
from snpdb.models import UserPreview


@receiver(signal=health_check_signal)
def active_users_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    event_log_activity = Event.objects.filter(date__gte=health_request.since, date__lt=health_request.now).values_list('user', flat=True).distinct()
    view_event_activity = ViewEvent.objects.filter(created__gte=health_request.since, created__lt=health_request.now).values_list('user', flat=True).distinct()

    active_users_queryset = User.objects.filter(Q(pk__in=event_log_activity) | Q(pk__in=view_event_activity)).exclude(groups=bot_group()).order_by('username')
    active_users = active_users_queryset.values_list('username', flat=True)

    count = active_users.count()
    emoji: str
    if count == 0:
        emoji = ":ghost:"
    else:
        emojis = [":nerd_face:", ":thinking_face:", ":face_with_monocle:", ":face_with_cowboy_hat:"]
        emoji = emojis[randint(0, len(emojis) - 1)]

    user_previews = [UserPreview(user).preview for user in active_users_queryset]

    return HealthCheckRecentActivity(
        emoji=emoji,
        name="Active Users",
        amount=count,
        extra=", ".join(list(active_users)),
        preview=user_previews,
        stand_alone=True  # always give active users its own line
    )


@receiver(signal=health_check_signal)
def email_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    from email_manager.models import EmailLog
    recent_emails = EmailLog.objects.filter(created__gte=health_request.since, created__lt=health_request.now)

    if count := recent_emails.count():
        subject_counts = Counter(email.subject for email in recent_emails)

        email_subj = ", ".join(
            f"{subject} x *{count}*" if count > 1 else subject for subject, count in subject_counts.items()
        )
        
        return HealthCheckRecentActivity(
            emoji=":email:",
            name="Emails Sent",
            amount=count,
            extra=email_subj,
            stand_alone=True,
            preview=[email.preview for email in recent_emails]
        )


@receiver(signal=health_check_signal)
def discordance_triage_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    recent_reports = DiscordanceReportTriage.objects.filter(
        modified__gte=health_request.since,
        modified__lt=health_request.now, user__isnull=False)

    if count := recent_reports.count():
        return HealthCheckRecentActivity(
            emoji=":triangular_ruler:",
            name="Discordance Triage Reports",
            amount=count,
            extra=", ".join([f'DR_{report.discordance_report.pk}' for report in recent_reports]),
            stand_alone=True,
            preview=[report.preview for report in recent_reports]
        )
