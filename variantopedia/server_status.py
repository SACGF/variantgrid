from datetime import timedelta
from typing import Optional

from django.contrib.auth.models import User
from django.utils import timezone
from django.utils.timesince import timesince

from analysis.models import Analysis
from classification.models import Classification
from eventlog.models import Event
from library.enums.log_level import LogLevel
from library.guardian_utils import bot_group
from snpdb.models import UserPageAck, VCF


def get_dashboard_notices(user: User, days_ago: int) -> dict:
    """ returns {} if nothing to show """
    MAX_PAST_DAYS = 30

    start_time = 0
    notice_header = ""

    if days_ago:
        # if days ago was passed in, don't update upa
        days_ago = min(days_ago, MAX_PAST_DAYS)
        start_time = timezone.now() - timedelta(days=days_ago)
        notice_header = f"Since the last {days_ago} day{'s' if days_ago > 1 else ''}"
    # else:
    #     upa, created = UserPageAck.objects.get_or_create(user=user, page_id="server_status")
    #     if created:
    #         if user.last_login:
    #             start_time = user.last_login
    #             notice_header = f"Since last login ({timesince(start_time)} ago)"
    #     else:
    #         start_time = upa.modified
    #         notice_header = f"Since last visit to this page ({timesince(start_time)} ago)"
    #         upa.save()  # Update last modified timestamp
    #
    #     max_days_ago = timezone.now() - timedelta(days=MAX_PAST_DAYS)
    #     if max_days_ago > start_time:
    #         start_time = max_days_ago
    #         notice_header = f"Since the last {MAX_PAST_DAYS} days"

    if user.is_superuser:
        events = Event.objects.filter(date__gte=start_time, severity=LogLevel.ERROR)
    else:
        events = Event.objects.none()

    active_users = list(User.objects.filter(pk__in=Event.objects.filter(date__gte=start_time).exclude(user__groups=bot_group()).values_list('user', flat=True).distinct()).order_by('username').values_list('username', flat=True))

    vcfs = VCF.filter_for_user(user, True).filter(date__gte=start_time)
    analyses = Analysis.filter_for_user(user)
    analyses_created = analyses.filter(created__gte=start_time)
    analyses_modified = analyses.filter(created__lt=start_time, modified__gte=start_time)
    classifications_of_interest = Classification.dashboard_report_classifications_of_interest(since=start_time)
    new_classification_count = Classification.dashboard_report_new_classifications(since=start_time)

    return {
        "active_users": active_users,
        "notice_header": notice_header,
        "events": events,
        "classifications_of_interest": classifications_of_interest,
        "classifications_created": new_classification_count,
        "vcfs": vcfs,
        "analyses_created": analyses_created,
        "analyses_modified": analyses_modified,
    }
