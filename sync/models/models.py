from typing import Optional

from django.conf import settings
from django.db import models
from django.db.models.deletion import CASCADE
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
