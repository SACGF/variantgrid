"""
Generic pre-drop archival pipeline for any RelatedModelsPartitionModel subclass.

The helper kicks off the long-running Celery task and returns immediately;
the task does the actual pg_dump -> verify -> drop sequence.

@see claude/issue_1537_archive_plan.md
"""

import os
from datetime import datetime

from django.contrib.auth.models import User
from django.db import transaction
from django.utils import timezone

from eventlog.models import create_event
from library.django_utils.data_archive_mixin import DataArchiveMixin
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.django_utils.partition_archive_paths import get_partition_archive_path
from snpdb.models.models_partition_archive import PartitionArchive


class PartitionArchivePreconditionError(Exception):
    """ Raised when a partition archive can't be scheduled (already archived, wrong shape). """


def archive_partitioned_model(source: RelatedModelsPartitionModel,
                              user: User,
                              reason: str = "") -> PartitionArchive:
    """ Schedule a partition-dump+drop for `source`.

        Creates a PENDING PartitionArchive row and queues the Celery task.
        Caller (admin action) returns immediately; the task does the work.
    """
    if not isinstance(source, DataArchiveMixin):
        raise PartitionArchivePreconditionError(
            f"{type(source).__name__} does not use DataArchiveMixin -- cannot stamp archive metadata"
        )
    if source.data_archived:
        raise PartitionArchivePreconditionError(f"{source} is already archived")

    table_names = list(source.RECORDS_BASE_TABLE_NAMES)
    if not table_names:
        raise PartitionArchivePreconditionError(f"{source} declares no RECORDS_BASE_TABLE_NAMES")

    meta = source._meta
    in_progress_qs = PartitionArchive.objects.filter(
        source_app_label=meta.app_label,
        source_model=meta.object_name,
        source_pk=str(source.pk),
        status__in=[PartitionArchive.Status.PENDING, PartitionArchive.Status.IN_PROGRESS],
    )
    if in_progress_qs.exists():
        raise PartitionArchivePreconditionError(
            f"{source} already has a PENDING/IN_PROGRESS PartitionArchive row"
        )

    iso = datetime.utcnow().strftime("%Y-%m-%dT%H-%M-%S")
    archive_name = f"{meta.app_label}_{meta.model_name}_{source.pk}_{iso}"
    dump_path = get_partition_archive_path(archive_name + ".dump")

    with transaction.atomic():
        archive = PartitionArchive.objects.create(
            archive_name=archive_name,
            dump_path=dump_path,
            source_app_label=meta.app_label,
            source_model=meta.object_name,
            source_pk=str(source.pk),
            source_table_names=table_names,
            status=PartitionArchive.Status.PENDING,
            dumped_by=user,
        )

    transaction.on_commit(lambda: _schedule_archive_task(archive.pk, reason))
    return archive


def _schedule_archive_task(archive_pk: int, reason: str):
    from snpdb.tasks.partition_archive_tasks import perform_partition_archive
    perform_partition_archive.delay(archive_pk, reason)


def mark_partition_archive_restored(archive: PartitionArchive, user: User) -> None:
    if archive.status != PartitionArchive.Status.COMPLETE:
        raise ValueError(
            f"Archive status is {archive.status} -- only COMPLETE archives can be marked restored"
        )
    if archive.restored_date is not None:
        raise ValueError("Archive is already marked restored")

    archive.restored_date = timezone.now()
    archive.restored_by = user
    archive.save(update_fields=["restored_date", "restored_by"])

    source = archive.resolve_source()
    if source is not None and getattr(source, "data_archived", False):
        source.data_archived_date = None
        source.data_archived_by = None
        source.data_archive_reason = None
        source.data_restorable_from = None
        source.save(update_fields=["data_archived_date", "data_archived_by",
                                   "data_archive_reason", "data_restorable_from"])

    create_event(user, "partition_archive_restored",
                 details=f"archive_pk={archive.pk} archive_name={archive.archive_name!r}",
                 app_name="snpdb")


def clear_partition_archive_dump(archive: PartitionArchive, user: User) -> None:
    if archive.cleared_date is not None:
        raise ValueError("Dump is already cleared")
    if archive.dump_path and os.path.exists(archive.dump_path):
        os.unlink(archive.dump_path)
    archive.cleared_date = timezone.now()
    archive.cleared_by = user
    archive.save(update_fields=["cleared_date", "cleared_by"])
    create_event(user, "partition_archive_cleared",
                 details=f"archive_pk={archive.pk} archive_name={archive.archive_name!r}",
                 app_name="snpdb")
