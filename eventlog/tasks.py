"""
Scheduled housekeeping for the event log: prune_api_request_log (daily, from variantgrid/celery.py).
"""
from datetime import timedelta

import celery
from django.conf import settings
from django.utils import timezone

from eventlog.models import ApiRequestLog


@celery.shared_task(queue='db_workers')
def prune_api_request_log():
    cutoff = timezone.now() - timedelta(days=settings.API_REQUEST_LOG_DAYS)
    ApiRequestLog.objects.filter(created__lt=cutoff).delete()
