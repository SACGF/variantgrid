from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
from django_extensions.db.models import TimeStampedModel

from classification.models.classification import ClassificationModification
from sync.models.models import SyncRun, SyncDestination


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
        up_to_date = ClassificationModificationSyncRecord.objects.filter(run__destination=destination, success=True).values_list('classification_modification', flat=True)
        return qs.exclude(pk__in=up_to_date)
