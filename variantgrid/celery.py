# noinspection HttpUrlsUsage
"""
    Set default Django settings module for the 'celery' program. Setup as per:
    http://celery.readthedocs.org/en/latest/django/first-steps-with-django.html
"""

import logging
import os

import celery
import rollbar
from celery import Celery
from celery.schedules import crontab
from django.conf import settings

from library.constants import HOUR_SECS, MINUTE_SECS
from library.django_utils.rollbar_middleware import RollbarIgnoreException

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'variantgrid.settings')


app = Celery('variantgrid')

# Using a string here means the worker will not have to
# pickle the object when using Windows.
app.config_from_object('django.conf:settings', namespace='CELERY')
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

app.conf.timezone = settings.TIME_ZONE

# Note: I'm having timezone issues, crontab doesn't work but schedule seconds does

app.conf.beat_schedule = {}

# Runs Sync as declared in the Sync table
# As of 2019-July it's only uploading to shariant
SYNC_DETAILS = getattr(settings, "SYNC_DETAILS", {})
if SYNC_DETAILS and any(sd["enabled"] for sd in SYNC_DETAILS.values()):
    app.conf.beat_schedule['sync_data'] = {
        'task': 'sync.tasks.sync_tasks.sync_all',
        'schedule': HOUR_SECS,
    }

# TODO - move this into settings???
if all([settings.SEQAUTO_ENABLED, settings.SEQAUTO_SCAN_DISKS, settings.UPLOAD_ENABLED]):
    app.conf.beat_schedule['seqauto-nightly-scan'] = {
        'task': 'seqauto.tasks.scan_run_jobs.scan_run_jobs',
        'schedule': HOUR_SECS * 2,
    }

SAPATH_ENABLED = any((a.startswith("sapath") for a in settings.INSTALLED_APPS))
if SAPATH_ENABLED:
    helix_user = getattr(settings, "SAPATH_HELIX_USER", None)
    if helix_user:
        app.conf.beat_schedule['sapath-helix-load-if-changed'] = {
            'task': 'sapath.tasks.import_helix_task.sapath_helix_load_if_changed',
            'schedule': HOUR_SECS,  # Check every hour, only update if hash changed
        }

app.conf.beat_schedule['notify-server-status'] = {
    'task': 'variantopedia.tasks.server_status_tasks.notify_server_status',
    'schedule': crontab(hour=19, minute=0),
}


# send update emails once a day (if there has been activity)
if settings.DISCORDANCE_EMAIL:
    app.conf.beat_schedule['discordance-emails-weekly'] = {
        'task': 'classification.views.classification_email_view.send_summary_emails',
        'schedule': crontab(hour=10, minute=0, day_of_week='mon')
    }

# Server monitoring tasks - send RollBar warnings
if settings.SERVER_MIN_DISK_WARNING_GIGS:
    app.conf.beat_schedule['warn-low-disk-space'] = {
        'task': 'variantgrid.tasks.server_monitoring_tasks.warn_low_disk_space',
        'schedule': HOUR_SECS,
    }

# Ping an uptime service monitor every 30 minutes to say the schedular is still running
if settings.HEARTBEAT_URL:
    app.conf.beat_schedule['heartbeat'] = {
        'task': 'variantgrid.tasks.server_monitoring_tasks.heartbeat',
        'schedule': MINUTE_SECS * 30,
    }


@app.task(bind=True)
def debug_task(self):
    logging.info('Request: %r', self.request)


@app.task(bind=True)
def fail_task(self):
    raise RollbarIgnoreException("Born to lose, I've lived my life in vain")


# Rollbar https://github.com/rollbar/rollbar-celery-example
rollbar.init(**settings.ROLLBAR)


def celery_base_data_hook(request, data):
    data['framework'] = 'celery'


rollbar.BASE_DATA_HOOK = celery_base_data_hook


@celery.signals.task_failure.connect
def on_task_failure(**kwargs):
    if exception := kwargs.get("exception"):
        if isinstance(exception, RollbarIgnoreException):
            return

    rollbar.report_exc_info(extra_data=kwargs)
