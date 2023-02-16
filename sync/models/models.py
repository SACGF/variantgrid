from typing import List, Optional

from django.conf import settings
from django.db import models
from django.db.models.deletion import CASCADE
from django.utils.timesince import timesince
from django_extensions.db.models import TimeStampedModel

from sync.models.enums import SyncStatus


class SyncDestination(models.Model):
    """
    Configuration for a sync.
    Examples
    - syncing records up to shariant
    - syncing records down from shariant
    - syncing records to clinvar (in future)
    """

    name = models.TextField(null=False, unique=True)
    config = models.JSONField(null=False, blank=True, default=dict)
    enabled = models.BooleanField(null=False, blank=True, default=True)

    def run(self, full_sync: bool = False, max_rows: Optional[int] = None):
        from sync.sync_runner import run_sync
        run_sync(self, full_sync=full_sync, max_rows=max_rows)

    @property
    def sync_details(self) -> dict:
        key = self.config["sync_details"]
        sd = settings.SYNC_DETAILS.get(key)
        if not sd:
            actual_keys = ", ".join(settings.SYNC_DETAILS.keys()) if settings.SYNC_DETAILS else "No options configured"
            raise ValueError(f"Could not find sync details '{key}' options are: {actual_keys}")
        return sd

    def report(self) -> Optional[str]:
        """ Some environments may temporarily disable them (in which case we want to know)
            So report if any have ever tried. If never tried return None """
        sd_info = None
        last_attempt = self.syncrun_set.order_by('-created').first()
        if self.enabled or last_attempt:
            if last_attempt:
                time_since_last_attempt = timesince(last_attempt.modified)
            else:
                time_since_last_attempt = "never"

            if last_success := self.syncrun_set.filter(status=SyncStatus.SUCCESS).order_by('-created').first():
                time_since_last_success = timesince(last_success.modified)
            else:
                time_since_last_success = "never"

            sd_info = f"{self} - last success: {time_since_last_success}"
            if last_attempt != last_success:
                sd_info += f", last attempt: {time_since_last_attempt}"
            if not self.enabled:
                sd_info += " (*currently disabled*)"

        return sd_info

    @staticmethod
    def get_reports() -> List:
        reports = []
        for sd in SyncDestination.objects.order_by("pk").all():
            if report := sd.report():
                reports.append(report)
        return reports

    def __str__(self):
        return self.name


class SyncRun(TimeStampedModel):
    """
    An instance of
    """
    destination = models.ForeignKey(SyncDestination, on_delete=CASCADE)
    status = models.CharField(max_length=1, choices=SyncStatus.choices, null=False)
    meta = models.JSONField(null=True, blank=True, default=None)

    def __str__(self):
        return f'{str(self.destination)} {str(self.created)}'
