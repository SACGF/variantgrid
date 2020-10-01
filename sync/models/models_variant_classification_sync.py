from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
from django_extensions.db.models import TimeStampedModel

from sync.models.models import SyncRun, SyncDestination
from classification.models.variant_classification import VariantClassificationModification


class VariantClassificationModificationSyncRecord(TimeStampedModel):
    """
    For tracking uploading/downloading of VariantClassificationModifications.
    Only create a record if the overall connection is working, e.g. if the username and password are incorrect
    do not create a failure for every one of these
    """
    run = models.ForeignKey(SyncRun, on_delete=CASCADE)
    variant_classification_modification = models.ForeignKey(VariantClassificationModification, on_delete=CASCADE)
    success = models.BooleanField(default=True)
    meta = models.JSONField(null=True, blank=True, default=None)

    @staticmethod
    def filter_out_synced(qs: QuerySet, destination: SyncDestination):
        up_to_date = VariantClassificationModificationSyncRecord.objects.filter(run__destination=destination, success=True).values_list('variant_classification_modification', flat=True)
        return qs.exclude(pk__in=up_to_date)
