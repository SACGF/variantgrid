from datetime import timedelta

from django.contrib.auth.models import User
from django.utils.timezone import localtime

from library.health_check import HealthCheckRequest, HealthCheckRecentActivity, health_check_signal
from library.utils import flatten_nested_lists


def get_dashboard_notices(user: User, days_ago: int) -> dict:
    """ returns {} if nothing to show """
    MAX_PAST_DAYS = 30
    notice_header = ""

    # if days ago was passed in, don't update upa
    days_ago = min(days_ago, MAX_PAST_DAYS)
    notice_header = f"Since the last {days_ago} day{'s' if days_ago > 1 else ''}"

    now = localtime()
    since = now - timedelta(days=days_ago)
    health_request = HealthCheckRequest(since=since, now=now)
    results = []
    previews = {}
    for _, result in health_check_signal.send_robust(sender=None, health_request=health_request):
        if not isinstance(result, Exception):
            if isinstance(result, list):
                results.extend(result)
                for activity in result:
                    if isinstance(activity, HealthCheckRecentActivity):
                        if activity.preview is not None:
                            previews.setdefault(activity.name, []).append(activity.preview)
            elif isinstance(result, HealthCheckRecentActivity):
                results.append(result)
                if result.preview is not None:
                    previews.setdefault(result.name, []).append(result.preview)

    checks = flatten_nested_lists(results)
    checks = sorted(checks, key=lambda hc: hc.sort_order())

    recent_lines_with_previews = []
    overall_lines = []

    for check in checks:
        line_content = check.as_html()
        line_previews = previews.get(check.name, [])
        if check.is_recent_activity():
            recent_lines_with_previews.append((line_content, line_previews))
        else:
            overall_lines.append(line_content)

    return {
        "notice_header": notice_header,
        'recent_lines_with_previews': checks,
        'overall_lines': overall_lines,
        'days_ago': days_ago
    }
