from django.db import models
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel
from library.utils import empty_dict
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
    config = models.JSONField(null=False, blank=True, default=empty_dict)

    def run(self, full_sync: bool = False):
        from sync.sync_runner import run_sync
        run_sync(self, full_sync=full_sync)

    def __str__(self):
        return self.name


class SyncRun(TimeStampedModel):
    """
    An instance of
    """
    destination = models.ForeignKey(SyncDestination, on_delete=CASCADE)
    status = models.CharField(max_length=1, choices=SyncStatus.CHOICES, null=False)
    meta = models.JSONField(null=True, blank=True, default=None)

    def __str__(self):
        return f'{str(self.destination)} {str(self.created)}'
