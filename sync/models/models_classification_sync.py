import re
from typing import Optional
from urllib.parse import urljoin

from django.db import models
from django.db.models import Q
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
from django_extensions.db.models import TimeStampedModel

from classification.models.classification import Classification, ClassificationModification
from sync.models.models import SyncDestination, SyncRun

_REMOTE_PK_PATTERN = re.compile(r'[\w\-]+')


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
    def filter_out_synced(qs: QuerySet[ClassificationModification], destination: SyncDestination) -> QuerySet[ClassificationModification]:
        cmsr_qs = ClassificationModificationSyncRecord.objects.filter(run__destination=destination, success=True)
        synced = cmsr_qs.values_list('classification_modification', flat=True)
        # withdrawing makes no new modification, so a withdrawn record's synced version stays due until the
        # remote answers for it as withdrawn (or deleted, for a record it held below all users)
        withdrawal_synced = cmsr_qs.filter(Q(meta__withdrawn=True) | Q(meta__deleted=True)) \
            .values_list('classification_modification', flat=True)
        up_to_date = Q(pk__in=synced) & (Q(classification__withdrawn=False) | Q(pk__in=withdrawal_synced))
        return qs.exclude(up_to_date)

    @property
    def remote_url(self) -> Optional[str]:
        """ URL of the record on the destination server, None if we don't know where it landed """
        if self.success:
            # the v2 API response nests the remote record's pk under "meta"
            remote_pk = ((self.meta or {}).get("meta") or {}).get("id")
            # the pk comes from the remote server, so don't let it add anything but a path segment to the link
            if remote_pk and _REMOTE_PK_PATTERN.fullmatch(str(remote_pk)):
                host = self.run.destination.sync_details["host"]
                base = host if host.endswith('/') else host + '/'
                return urljoin(base, Classification.get_url_for_pk(remote_pk).lstrip('/'))
        return None
