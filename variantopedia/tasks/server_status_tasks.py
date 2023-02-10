
from random import randint

import celery
from django.conf import settings
from django.urls import reverse

from library.django_utils import get_url_from_view_path
from library.health_check import populate_health_check
from library.log_utils import NotificationBuilder


@celery.shared_task
def notify_server_status():
    if not settings.HEALTH_CHECK_ENABLED:
        return
    notify_server_status_now()


def notify_server_status_now():
    url = get_url_from_view_path(reverse('server_status')) + '?activeTab=server_status_activity_detail_1'
    nb = NotificationBuilder(message="Health Check")
    heading_emoji = ":male-doctor:" if randint(0, 1) else ":female-doctor:"
    nb.add_header(f"{heading_emoji} Health Check")
    nb.add_markdown(f"*In the <{url}|last 24 hours>*")
    populate_health_check(nb)
    nb.send()
