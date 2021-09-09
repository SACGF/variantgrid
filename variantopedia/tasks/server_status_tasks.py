from random import randint

import celery
from django.conf import settings
from django.urls import reverse

from classification.models import Classification
from library.django_utils import get_url_from_view_path
from library.guardian_utils import admin_bot
from library.log_utils import NotificationBuilder
from library.utils import pretty_label, count
from variantgrid.perm_path import get_visible_url_names
from variantgrid.tasks.server_monitoring_tasks import get_disk_messages
from variantopedia.server_status import get_dashboard_notices


@celery.shared_task
def notify_server_status():
    if not settings.HEALTH_CHECK_ENABLED:
        return

    dashboard_notices = get_dashboard_notices(admin_bot(), days_ago=1)
    url = get_url_from_view_path(reverse('server_status')) + '?days=1'

    emoji = ":male-doctor:" if randint(0, 1) else ":female-doctor:"
    nb = NotificationBuilder(message="Health Check", emoji=emoji)
    nb.add_header("Health Check")
    nb.add_markdown(f"URL : <{url}|{url}>")
    nb.add_markdown("*Disk usage*")
    disk_usage = list()
    for _, message in get_disk_messages(info_messages=True):
        disk_usage.append(f":floppy_disk: {message}")
    nb.add_markdown("\n".join(disk_usage), indented=True)
    nb.add_markdown("*In the last 24 hours*")

    keys = set(dashboard_notices.keys())
    keys.discard('events')
    keys.discard('notice_header')
    visible_urls = get_visible_url_names()
    if not visible_urls.get('analyses'):
        for exclude_key in ['vcfs', 'analyses_created', 'analyses_modified']:
            keys.discard(exclude_key)
    sorted_keys = sorted(list(keys))

    lines = list()

    for key in sorted_keys:
        values = dashboard_notices.get(key)
        count_display = count(values)

        display_individuals = False
        emoji = ":blue_book:"
        if 'analyses' in key:
            emoji = ":orange_book:"
        elif 'vcf' in key:
            emoji = ":green_book:"
        elif 'active_users' in key:
            if count_display == 0:
                emoji = ":ghost:"
            else:
                emojis = [":nerd_face:", ":thinking_face:", ":face_with_monocle:", ":face_with_cowboy_hat:"]
                emoji = emojis[randint(0, len(emojis) - 1)]
                display_individuals = True

        if count_display:
            count_display = f"*{count_display}*"

        line = f"{emoji} {count_display} : {pretty_label(key)}"
        if display_individuals:
            line = f"{line} : {', '.join(str(value) for value in values)}"

        lines.append(line)
    nb.add_markdown("\n".join(lines), indented=True)
    nb.add_markdown("*In Total*")
    nb.add_markdown(
        f":people_holding_hands: {Classification.dashboard_total_shared_classifications():,} : Classifications Shared" +
        f"\n:cry: {Classification.dashboard_total_unshared_classifications():,} : Classifications Un-Shared",
        indented=True
    )
    nb.send()
