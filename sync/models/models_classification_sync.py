from typing import Optional
from urllib.parse import urljoin

from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
from django_extensions.db.models import TimeStampedModel

from classification.models.classification import Classification, ClassificationModification
from sync.models.models import SyncDestination, SyncRun


class ClassificationModificationSyncRecord(TimeStampedModel):
    """
    For tracking uploading/downloading of ClassificationModifications.
    Only create a record if the overall connection is working, e.g. if the username and password are incorrect
    do not create a failure for every one of these
    """
    run = models.ForeignKey(SyncRun, on_delete=CASCADE)
    classification_modification = models.ForeignKey(ClassificationModification, on_delete=CASCADE)
    success = models.BooleanField(default=True)
    meta = models.JSONField(null=True, blank=True, default=None)

    @staticmethod
    def filter_out_synced(qs: QuerySet, destination: SyncDestination):
        cmsr_qs = ClassificationModificationSyncRecord.objects.filter(run__destination=destination, success=True)
        up_to_date = cmsr_qs.values_list('classification_modification', flat=True)
        return qs.exclude(pk__in=up_to_date)

    @property
    def remote_url(self) -> Optional[str]:
        """ URL of the record on the destination server, None if we don't know where it landed """
        if self.success:
            # the v2 API response nests the remote record's pk under "meta"
            if remote_pk := ((self.meta or {}).get("meta") or {}).get("id"):
                host = self.run.destination.sync_details["host"]
                return urljoin(host, Classification.get_url_for_pk(remote_pk))
        return None
