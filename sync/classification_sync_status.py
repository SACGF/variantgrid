"""
Reports where a single classification stands with regards to each upload SyncDestination,
so the classification page can say what actually reached a destination rather than inferring
it from the share level.
"""
import logging
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional

from classification.enums.classification_enums import ShareLevel
from classification.models.classification import Classification, ClassificationModification
from sync.models.models import SyncDestination
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.sync_runner import ClassificationUploadSyncRunner, sync_runner_for_destination


class ClassificationSyncState(str, Enum):
    SYNCED = "synced"
    CHANGES_PENDING = "changes_pending"
    PENDING = "pending"
    EXCLUDED = "excluded"

    @property
    def label(self) -> str:
        return {
            ClassificationSyncState.SYNCED: "Uploaded",
            ClassificationSyncState.CHANGES_PENDING: "Uploaded, later changes upload on the next sync",
            ClassificationSyncState.PENDING: "Uploads on the next sync",
            ClassificationSyncState.EXCLUDED: "Held locally",
        }[self]

    @property
    def icon(self) -> str:
        return {
            ClassificationSyncState.SYNCED: "fas fa-check-circle text-success",
            ClassificationSyncState.CHANGES_PENDING: "fas fa-check-circle text-success",
            ClassificationSyncState.PENDING: "fas fa-clock text-secondary",
            ClassificationSyncState.EXCLUDED: "fas fa-info-circle text-secondary",
        }[self]


@dataclass(frozen=True)
class ClassificationSyncDestinationStatus:
    destination_name: str
    state: ClassificationSyncState
    remote_url: Optional[str] = None
    last_synced: Optional[datetime] = None
    reasons: list[str] = field(default_factory=list)

    @property
    def icon(self) -> str:
        return self.state.icon

    @property
    def label(self) -> str:
        return self.state.label

    @property
    def note(self) -> Optional[str]:
        """ Warns when the remote link points at a copy the record would no longer be sent as """
        if self.state == ClassificationSyncState.EXCLUDED and self.remote_url:
            return "The remote record is from an earlier upload - the record as it stands now would not be sent."
        return None


def classification_sync_status(classification: Classification) -> list[ClassificationSyncDestinationStatus]:
    """ One status per enabled upload destination that selects local classifications """
    if classification.lab.external:
        return []

    last_published = classification.last_published_version
    if not last_published or last_published.share_level not in ShareLevel.DISCORDANT_LEVEL_KEYS:
        # the share level message already says the record stays with the lab, no need to repeat it
        return []

    statuses = []
    for sync_destination in SyncDestination.objects.filter(config__direction='upload', enabled=True).order_by('name'):
        try:
            runner = sync_runner_for_destination(sync_destination)
            if not isinstance(runner, ClassificationUploadSyncRunner):
                continue
            runner.configure(sync_destination)
            statuses.append(_status_for_destination(last_published, sync_destination, runner))
        except ValueError:
            # a destination we can't resolve a runner or host for is one we can't report on
            logging.warning("Could not report sync status for %s", sync_destination, exc_info=True)
    return statuses


def _last_successful_sync_record(sync_destination: SyncDestination, **kwargs) -> Optional[ClassificationModificationSyncRecord]:
    return ClassificationModificationSyncRecord.objects \
        .filter(run__destination=sync_destination, success=True, **kwargs) \
        .order_by('-created').first()


def _status_for_destination(
        last_published: ClassificationModification,
        sync_destination: SyncDestination,
        runner: ClassificationUploadSyncRunner) -> ClassificationSyncDestinationStatus:

    if sync_record := _last_successful_sync_record(sync_destination, classification_modification=last_published):
        return ClassificationSyncDestinationStatus(
            destination_name=sync_destination.name,
            state=ClassificationSyncState.SYNCED,
            remote_url=sync_record.remote_url,
            last_synced=sync_record.created
        )

    would_sync = runner.records_to_sync(full_sync=True).filter(pk=last_published.pk).exists()
    earlier_record = _last_successful_sync_record(
        sync_destination, classification_modification__classification=last_published.classification)

    if would_sync:
        if earlier_record:
            return ClassificationSyncDestinationStatus(
                destination_name=sync_destination.name,
                state=ClassificationSyncState.CHANGES_PENDING,
                remote_url=earlier_record.remote_url,
                last_synced=earlier_record.created
            )
        return ClassificationSyncDestinationStatus(
            destination_name=sync_destination.name,
            state=ClassificationSyncState.PENDING,
            remote_url=runner.remote_url_for(last_published)
        )

    return ClassificationSyncDestinationStatus(
        destination_name=sync_destination.name,
        state=ClassificationSyncState.EXCLUDED,
        remote_url=earlier_record.remote_url if earlier_record else None,
        last_synced=earlier_record.created if earlier_record else None,
        reasons=runner.exclusion_reasons(last_published)
    )
