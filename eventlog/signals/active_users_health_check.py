from collections import Counter

from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver

from classification.models import DiscordanceReportTriage, ClinVarExportBatch, ClinVarExportBatchStatus
from eventlog.models import Event, ViewEvent
from library.guardian_utils import bot_group
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckRecentActivity
from snpdb.models import UserPreview


@receiver(signal=health_check_signal)
def active_users_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    event_log_activity = Event.objects.filter(date__gte=health_request.since, date__lt=health_request.now).values_list('user', flat=True).distinct()
    view_event_activity = ViewEvent.objects.filter(created__gte=health_request.since, created__lt=health_request.now).values_list('user', flat=True).distinct()

    superusers_queryset = User.objects.filter(
        Q(is_superuser=True) | Q(groups__name='bot'),
        Q(pk__in=event_log_activity) | Q(pk__in=view_event_activity)
    ).distinct().order_by('username')

    non_admin_queryset = User.objects.filter(Q(pk__in=event_log_activity) | Q(pk__in=view_event_activity)).exclude(is_superuser=True).exclude(groups=bot_group()).order_by('username')

    superusers = superusers_queryset.values_list('username', flat=True)
    non_admin_users = non_admin_queryset.values_list('username', flat=True)

    superusers_count = superusers.count()
    non_admin_count = non_admin_users.count()

    superusers_emoji = ":crown:" if superusers_count > 0 else ":ghost:"
    non_admin_emoji = ":nerd_face:" if non_admin_count > 0 else ":ghost:"

    superusers_previews = [UserPreview(user).preview for user in superusers_queryset]
    non_admin_previews = [UserPreview(user).preview for user in non_admin_queryset]

    health_checks = []

    if superusers_count > 0:
        health_checks.append(
            HealthCheckRecentActivity(
                emoji=superusers_emoji,
                name="Active Admin Users",
                amount=superusers_count,
                extra=", ".join(list(superusers)),
                preview=superusers_previews,
                stand_alone=True  # always give superusers their own line
            )
        )

    if non_admin_count > 0:
        health_checks.append(
            HealthCheckRecentActivity(
                emoji=non_admin_emoji,
                name="Active Users",
                amount=non_admin_count,
                extra=", ".join(list(non_admin_users)),
                preview=non_admin_previews,
                stand_alone=True
            )
        )
    return health_checks


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


@receiver(signal=health_check_signal)
def clinvar_export_batch_healthcheck(sender, health_request: HealthCheckRequest, **kwargs):

    accepted_statuses = [ClinVarExportBatchStatus.AWAITING_UPLOAD, ClinVarExportBatchStatus.UPLOADING]
    recent_batches = ClinVarExportBatch.objects.filter(status__in=accepted_statuses,
                                                       created__gte=health_request.since,
                                                       created__lt=health_request.now)

    if count := recent_batches.count():
        return HealthCheckRecentActivity(
            emoji=":package:",
            name="ClinVar Export Batches",
            amount=count,
            extra=", ".join([f'{batch.clinvar_key.name} - *{batch.get_status_display()}*' for batch in recent_batches]),
            stand_alone=True,
        )
