"""
Archive/restore for GeneCoverageCollection.

@see claude/issue_1536_data_archive_plan.md
"""

import logging
import os

from django.contrib.auth.models import User
from django.utils import timezone

from eventlog.models import create_event
from genes.models import GeneCoverageCollection
from snpdb.archive import ArchivePreconditionError
from snpdb.models import DataState


def archive_gene_coverage_collection(gcc: GeneCoverageCollection, user: User, reason: str = "") -> None:
    """ Drop GeneCoverage* partition data, stamp the mixin. Keeps the row. """
    if gcc.data_archived:
        return
    upload_path = gcc.path
    if not upload_path or not os.path.exists(upload_path):
        raise ArchivePreconditionError(
            "Source coverage file doesn't exist — cannot archive. "
            "If you want to free up space you can permanently delete this collection."
        )

    logging.info("archive_gene_coverage_collection: dropping partition for GCC %s", gcc.pk)
    gcc.delete_related_objects()  # drops partition tables (no recreate)

    gcc.data_archived_date = timezone.now()
    gcc.data_archived_by = user
    gcc.data_archive_reason = reason
    gcc.data_restorable_from = upload_path
    gcc.save()

    create_event(user, "gene_coverage_collection_archived",
                 details=f"gcc_id={gcc.pk} restorable_from={upload_path!r} reason={reason!r}")


def restore_gene_coverage_collection(gcc: GeneCoverageCollection, user: User):
    """ Re-import an archived GeneCoverageCollection from its recorded source path. """
    from genes.tasks.gene_coverage_tasks import reload_gene_coverage_collection

    if not gcc.data_archived:
        raise ValueError("GeneCoverageCollection is not archived")
    if not gcc.data_restorable_from:
        raise ValueError("No restore path recorded")
    if not os.path.exists(gcc.data_restorable_from):
        raise ValueError(f"Restore source missing: {gcc.data_restorable_from}")

    # Recreate the partition tables before reload writes into them.
    gcc.path = gcc.data_restorable_from
    gcc.create_partition()
    gcc.data_archived_date = None
    gcc.data_archived_by = None
    gcc.data_archive_reason = None
    gcc.data_restorable_from = None
    gcc.data_state = DataState.RUNNING
    gcc.save()

    create_event(user, "gene_coverage_collection_restored",
                 details=f"gcc_id={gcc.pk} restored_from={gcc.path!r}")
    return reload_gene_coverage_collection.delay(gcc.pk)
