"""
Archive/restore for VCF (and exception types reused across apps).

@see claude/issue_1536_data_archive_plan.md
"""

import logging
import os
from typing import Optional

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.db import transaction
from django.utils import timezone

from eventlog.models import create_event
from snpdb.models import VCF
from snpdb.models.models_enums import ImportStatus
from snpdb.variant_zygosity_count import update_all_variant_zygosity_counts_for_vcf


class DataArchivedError(Exception):
    """ Raised when code tries to read underlying data that has been archived.

        Callers either propagate to a UI surface (banner / badge) or catch in
        a source-node helper that maps to NodeStatus.ERROR_CONFIGURATION.
    """

    def __init__(self, obj):
        self.obj = obj
        date = getattr(obj, "data_archived_date", None)
        by = getattr(obj, "data_archived_by", None)
        reason = getattr(obj, "data_archive_reason", None) or ""
        date_str = date.strftime("%Y-%m-%d") if date else "unknown date"
        by_str = str(by) if by else "unknown user"
        msg = f"Source data archived on {date_str} by {by_str}"
        if reason:
            msg += f": {reason}"
        super().__init__(msg)


class ArchivePreconditionError(Exception):
    """ Raised when archive can't proceed (e.g. uploaded source file missing). """


def _resolve_vcf_upload_path(vcf: VCF) -> Optional[str]:
    """ Returns the filesystem path of the uploaded source VCF, or None if missing. """
    try:
        uf = vcf.uploadedvcf.uploaded_file
        return uf.get_filename()
    except (ObjectDoesNotExist, AttributeError):
        return None


def vcf_can_be_archived(vcf: VCF) -> bool:
    """ True if the uploaded source file still exists on disk. """
    if vcf.data_archived or vcf.data_archive_in_progress:
        return False
    path = _resolve_vcf_upload_path(vcf)
    return bool(path) and os.path.exists(path)


def mark_vcf_archive_started(vcf: VCF) -> None:
    """ Stamp the VCF as archive-in-progress before queueing the Celery task, so an
        immediate page reload reflects it and won't offer Archive/Restore again. """
    vcf.data_archive_started_date = timezone.now()
    vcf.save(update_fields=["data_archive_started_date"])


def check_vcf_archive_precondition(vcf: VCF, force: bool = False) -> Optional[str]:
    """ Validate that the VCF can be archived (cheap - safe to call synchronously
        from a web request before queueing the slow work).

        Returns the resolved upload path if the data will be restorable, else None
        (non-recoverable force archive). Raises ArchivePreconditionError if the
        uploaded source file is missing and force was not requested.
    """
    upload_path = _resolve_vcf_upload_path(vcf)
    restorable = bool(upload_path) and os.path.exists(upload_path)
    if not restorable and not force:
        raise ArchivePreconditionError(
            "Uploaded file doesn't exist — cannot archive. "
            "If you want to free up space you can permanently delete this VCF."
        )
    return upload_path if restorable else None


def archive_vcf(vcf: VCF, user: User, reason: str = "", force: bool = False) -> None:
    """ Drop CohortGenotype data, decrement zygosity counts, stamp the mixin.

        Keeps the VCF row, samples, and cohort intact. Sub-cohorts of this VCF
        will see DataArchivedError next time they hit cohort_genotype_collection.

        force=True archives even when the uploaded source file is missing. The
        data is dropped with no recorded restore source, so it cannot be
        restored later (non-recoverable). Used to free space for VCFs whose
        underlying files have already been deleted.

        This drops partition data and walks zygosity counts, which can take a long
        time - run it from a Celery task (see snpdb.tasks.vcf_archive_tasks) rather
        than inside a web request.
    """
    if vcf.data_archived:
        return
    upload_path = check_vcf_archive_precondition(vcf, force=force)
    restorable = upload_path is not None

    # Single transaction: if anything fails mid-flight, the whole thing rolls back to a
    # clean un-archived state rather than leaving data partially dropped but unstamped.
    # The "archived" stamp commits atomically with the data drop. The zygosity subtraction
    # is intentionally inside the transaction (and not swallowed): if a later step fails,
    # the subtraction rolls back too, so a retry can't decrement the global counts twice.
    with transaction.atomic():
        update_all_variant_zygosity_counts_for_vcf(vcf, '-')

        logging.info("archive_vcf: dropping internal data for VCF %s (restorable=%s)", vcf.pk, restorable)
        vcf.delete_internal_data(recreate_partitions=False)

        vcf.data_archived_date = timezone.now()
        vcf.data_archived_by = user
        vcf.data_archive_reason = reason
        vcf.data_restorable_from = upload_path if restorable else None
        vcf.save()

    # Only runs if the archive committed.
    # Refresh existing analyses so cached q-objects/node statuses reflect the archive.
    # reload_analysis_nodes (called per-analysis inside) bumps node.version → invalidates
    # the cache key and re-evaluates errors (which now include archive checks).
    from analysis.tasks.auto_analysis_tasks import reload_auto_analyses_for_vcf
    reload_auto_analyses_for_vcf.delay(vcf.pk)

    create_event(user, "vcf_archived",
                 details=f"vcf_id={vcf.pk} name={vcf.name!r} restorable_from={vcf.data_restorable_from!r} reason={reason!r}")


def restore_vcf(vcf: VCF, user: User):
    """ Re-import an archived VCF from its recorded source path.

        The success path (ImportGenotypeVCFSuccessTask) clears the mixin
        fields and emits the vcf_restored event.
    """
    from upload.uploaded_file_type import retry_upload_pipeline

    if not vcf.data_archived:
        raise ValueError("VCF is not archived")
    if not vcf.data_restorable_from:
        raise ValueError("No restore path recorded")
    if not os.path.exists(vcf.data_restorable_from):
        raise ValueError(f"Restore source missing: {vcf.data_restorable_from}")

    upload_pipeline = vcf.uploadedvcf.uploaded_file.uploadpipeline
    vcf.import_status = ImportStatus.IMPORTING
    vcf.save()
    for sample in vcf.sample_set.all():
        sample.import_status = ImportStatus.IMPORTING
        sample.save()
    return retry_upload_pipeline(upload_pipeline)
