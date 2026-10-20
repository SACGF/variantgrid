"""
AnnotationRuns nothing will ever dispatch (#1940): ones with no range lock (#1654), and ones on a HISTORICAL
VariantAnnotationVersion, which the dispatcher never leases. abandoned_annotation_runs_by_reason splits a run
queryset by why it is stuck (`vg status`, the deployment check); clean_up_abandoned_annotation_runs, run by the
clean_up_abandoned_annotation_runs command, deletes the unfinished ones that own no annotation rows and ERRORs the
rest. Promotion cancels a demoted version's runs itself (AnnotationRun.cancel_incomplete_runs).
"""
import logging
from collections import Counter

from django.db.models import QuerySet
from django.utils import timezone

from annotation.enums.annotation_enums import AnnotationStatus
from annotation.models.models import (
    VariantAnnotation,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion

ABANDONED_ERROR_EXCEPTION = "Abandoned: no range lock or HISTORICAL version, so it will never be dispatched"


def incomplete_annotation_runs_qs() -> QuerySet[AnnotationRun]:
    return AnnotationRun.objects.exclude(status__in=AnnotationStatus.get_completed_states())


def abandoned_annotation_runs_by_reason(qs: QuerySet[AnnotationRun]) -> dict[str, QuerySet[AnnotationRun]]:
    return {
        "no range lock": qs.filter(annotation_range_lock__isnull=True),
        "historical version": qs.filter(
            annotation_range_lock__version__status=VariantAnnotationVersion.Status.HISTORICAL),
    }


def annotation_run_owns_annotation_rows(annotation_run: AnnotationRun) -> bool:
    """ One run at a time on purpose: on a single annotation_run_id each partition uses its index, while over a
        list of ids the planner seq-scans 100M+ rows. """
    return any(klass.objects.filter(annotation_run=annotation_run).exists()
               for klass in (VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap))


def clean_up_abandoned_annotation_runs(dry_run: bool = False) -> Counter:
    """ Delete each unfinished abandoned run that owns no annotation rows, along with any range lock that
        leaves empty. A run that does own rows keeps them: an incomplete one is ERRORed, an ERROR one is left
        as it is. Runs under a live lease are skipped. Returns counts by outcome. """
    outcomes = Counter()
    unfinished_qs = AnnotationRun.objects.exclude(status=AnnotationStatus.FINISHED)
    emptied_range_lock_ids = set()
    for reason, qs in abandoned_annotation_runs_by_reason(unfinished_qs).items():
        for annotation_run in qs.order_by("pk"):
            if annotation_run.lease_expires and annotation_run.lease_expires >= timezone.now():
                outcome = "skipped (live lease)"
            elif not annotation_run_owns_annotation_rows(annotation_run):
                outcome = "deleted"
                if not dry_run:
                    if annotation_run.annotation_range_lock_id:
                        emptied_range_lock_ids.add(annotation_run.annotation_range_lock_id)
                    annotation_run.delete()
            elif annotation_run.status != AnnotationStatus.ERROR:
                outcome = "marked error"
                if not dry_run:
                    annotation_run.error_exception = ABANDONED_ERROR_EXCEPTION
                    annotation_run.save()
            else:
                outcome = "left (error with annotation rows)"
            outcomes[f"{reason}: {outcome}"] += 1

    if emptied_range_lock_ids:
        empty_locks_qs = AnnotationRangeLock.objects.filter(pk__in=emptied_range_lock_ids, annotationrun__isnull=True)
        num_deleted, _ = empty_locks_qs.delete()
        outcomes["range locks deleted"] = num_deleted
    logging.info("Abandoned annotation runs%s: %s", " (dry run)" if dry_run else "", dict(outcomes))
    return outcomes
