from datetime import timedelta
from django.contrib.auth.models import User
from django.utils.timezone import localtime
from library.health_check import HealthCheckRecentActivity, HealthCheckRequest, health_check_signal


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
