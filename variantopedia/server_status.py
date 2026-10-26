import itertools
from datetime import datetime, timedelta
from typing import Optional

from django.contrib.auth.models import User
from django.utils.timezone import localtime

from eventlog.notifications import NotificationBuilder
from library.health_check import (
    HealthCheckAge,
    HealthCheckRecentActivity,
    HealthCheckRequest,
    HealthCheckStat,
    health_check_overall_stats_signal,
    health_check_signal,
)
from library.utils import flatten_nested_lists


def get_dashboard_notices(user: User, days_ago: int) -> dict:
    MAX_PAST_DAYS = 30
    days_ago = min(days_ago, MAX_PAST_DAYS)
    notice_header = f"Since the last {days_ago} day{'s' if days_ago > 1 else ''}"

    now = localtime()
    since = now - timedelta(days=days_ago)
    health_request = HealthCheckRequest(since=since, now=now)
    results: list[HealthCheckRecentActivity] = []

    for _, result in health_check_signal.send_robust(sender=None, health_request=health_request):
        if isinstance(result, list):
            for single in result:
                if isinstance(single, HealthCheckRecentActivity):
                    results.append(single)
        elif isinstance(result, HealthCheckRecentActivity):
            results.append(result)

    results = sorted(results, key=lambda hc: hc.sort_order())

    recent_lines_with_previews: list[HealthCheckRecentActivity] = []
    overall_lines: list[HealthCheckRecentActivity] = []

    for result in results:
        if result.is_recent_activity():
            if result.preview:
                recent_lines_with_previews.append(result)
        else:
            overall_lines.append(result)

    return {
        "notice_header": notice_header,
        'recent_lines_with_previews': recent_lines_with_previews,
        'overall_lines': overall_lines,
        'days_ago': days_ago
    }


def populate_health_check(notification: NotificationBuilder, since: Optional[datetime] = None):
    # calling send_robust will respond with a list of tuples of receiver, the result of calling the receiver
    # the results could be an Exception, a HealthCheckStat, or a list of HealthCheckStats
    # flatten the results, report all the exceptions, then present all the HealthChecks

    now = localtime()
    if since is None:
        since = now - timedelta(days=1)
    health_request = HealthCheckRequest(since=since, now=now)

    results = []
    for caller, result in health_check_signal.send_robust(sender=None, health_request=health_request):
        if isinstance(result, Exception):
            notification.add_markdown(f"Exception generating health check by {caller}: {result}")
        else:
            results.extend(result if isinstance(result, list) else [result])

    for caller, result in health_check_overall_stats_signal.send_robust(sender=None, health_request=health_request):
        if isinstance(result, Exception):
            notification.add_markdown(f"Exception generating health check by {caller}: {result}")
        else:
            results.extend(result if isinstance(result, list) else [result])

    checks: list[HealthCheckStat] = flatten_nested_lists(results)
    checks.sort(key=lambda hc: hc.sort_order())

    recent_lines = []
    overall_lines = []
    for check_type, checks_typed in itertools.groupby(checks, key=type):
        checks_typed_list = list(checks_typed)

        if issubclass(check_type, HealthCheckAge):
            age_lines = HealthCheckAge.to_lines(checks_typed_list, health_request)
            if check_type.is_recent_activity():
                recent_lines.extend(age_lines)
            else:
                overall_lines.extend(age_lines)
        else:
            section_lines = [check.as_markdown() for check in checks_typed_list]
            if check_type.is_recent_activity():
                recent_lines.extend(section_lines)
            else:
                overall_lines.extend(section_lines)

    if recent_lines:
        notification.add_markdown("\n".join(recent_lines), indented=True)

    if overall_lines:
        notification.add_markdown("*Overall*")
        notification.add_markdown("\n".join(overall_lines), indented=True)
