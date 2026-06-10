import logging
import sys

from django.db import transaction
from django.db.models.aggregates import Count, Max, Min
from django.utils import timezone

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationStatus,
    Variant,
    VariantAnnotationVersion,
)
from annotation.vep_annotation import get_vep_variant_annotation_version_kwargs
from library.django_utils import highest_pk
from snpdb.models.models_genome import GenomeBuild


def get_or_create_variant_annotation_version_from_current_vep(genome_build: GenomeBuild) -> tuple[VariantAnnotationVersion, bool]:
    kwargs = get_vep_variant_annotation_version_kwargs(genome_build)
    # New rows start as NEW; promotion to ACTIVE happens once tables are populated
    variant_annotation_version, created = VariantAnnotationVersion.objects.get_or_create(
        **kwargs, defaults={"status": VariantAnnotationVersion.Status.NEW})
    now = timezone.now()
    if created:
        logging.info("New Variant Annotation version created!")
        variant_annotation_version.created = now
    else:
        logging.info("Existing annotation ok...")
        variant_annotation_version.last_checked_date = now

    variant_annotation_version.save()
    return variant_annotation_version, created


def vav_diff_vs_kwargs(vav: VariantAnnotationVersion, vep_vav_kwargs) -> str:
    diff = {}
    for k, v in vep_vav_kwargs.items():
        existing_v = getattr(vav, k)
        if existing_v != v:
            diff[k] = (existing_v, v)

    diff_items = []
    for field, (db_val, vep_val) in diff.items():
        diff_items.append(f"{field} - {db_val=} != {vep_val=}")
    return ", ".join(diff_items)

def _get_unannotated_count_min_max(annotation_version, search_min: int,
                                   annotation_batch_min=None, annotation_batch_max=None):
    """ Get at most annotation_batch_size unannotated variants above search_min """

    if annotation_batch_max is None:
        annotation_batch_max = sys.maxsize

    # Original code used a sort then limit to get a certain number of variants above a minimum but this was slow
    # when the variant table was huge. So we attempt to do this quickly by grabbing blocks of the primary key.
    variant_max = highest_pk(Variant)

    unannotated_count = 0
    min_variant_id = None
    max_variant_id = None
    while True:
        remaining_max = annotation_batch_max - unannotated_count
        search_max = min(variant_max, search_min + remaining_max)
        if search_min > variant_max:
            break

        logging.debug("Searching for unannotated variants in range: %d-%d", search_min, search_max)
        qs = get_variants_qs_for_annotation(annotation_version, min_variant_id=search_min, max_variant_id=search_max)
        results = qs.aggregate(count=Count("id"), min_variant_id=Min("id"), max_variant_id=Max("id"))
        if results["count"]:
            unannotated_count += results["count"]
            if min_variant_id is None:
                min_variant_id = results["min_variant_id"]
            max_variant_id = results["max_variant_id"]
            if annotation_batch_min and unannotated_count >= annotation_batch_min:
                break

        search_min = search_max + 1

    return unannotated_count, min_variant_id, max_variant_id


def get_annotation_range_lock_and_unannotated_count(variant_annotation_version: VariantAnnotationVersion,
                                                    annotation_batch_min=None, annotation_batch_max=None):
    logging.info("AnnotationRangeLock: variant annotation version = %s", variant_annotation_version)

    data = AnnotationRangeLock.objects.filter(version=variant_annotation_version).aggregate(Max("max_variant"))
    max_locked_variant_id = data.get('max_variant__max') or 0
    qs = Variant.objects.all()
    if max_locked_variant_id:
        qs = qs.filter(pk__gt=max_locked_variant_id)

    # This query is fast (it's only when you join to tables looking for unannotated that it's slow)
    min_data = qs.aggregate(Min("pk"))
    # If max_locked_variant_id is the highest Variant on the system this will be None so use max + 1
    search_min = min_data["pk__min"] or max_locked_variant_id + 1

    annotation_version = variant_annotation_version.get_any_annotation_version()
    unannotated_variants_count, min_variant_id, max_variant_id = _get_unannotated_count_min_max(annotation_version,
                                                                                                search_min,
                                                                                                annotation_batch_min,
                                                                                                annotation_batch_max)

    annotation_range_lock = None
    if unannotated_variants_count:
        annotation_range_lock = AnnotationRangeLock(version=variant_annotation_version,
                                                    min_variant_id=min_variant_id,
                                                    max_variant_id=max_variant_id,
                                                    count=unannotated_variants_count)

    logging.info("AnnotationRangeLock: range: %s, count: %d", annotation_range_lock, unannotated_variants_count)
    return annotation_range_lock, unannotated_variants_count


def _range_lock_is_dispatchable(range_lock: AnnotationRangeLock, now=None) -> bool:
    """ A range lock is mergeable only while all of its runs are still pending (#2667):
        CREATED, un-leased, no task_id and not external. A lock with no runs yet (orphaned mid-create)
        is pure metadata and safe to merge too. """
    runs = list(range_lock.annotationrun_set.all())
    return all(run.is_dispatchable(now) for run in runs)


def _absorb_range_lock(survivor: AnnotationRangeLock, absorbed: AnnotationRangeLock):
    """ Extend `survivor` to cover `absorbed`'s range then delete `absorbed`. Wrapped in an atomic
        transaction with select_for_update so a crash can't leave `survivor` overlapping an
        un-deleted neighbour (#2667 'Transactionality of merge'). Mutates `survivor` in place.

        Synchronous + cheap: nothing has been dumped pre-launch, so a lock is pure metadata and the
        cascade-deleted CREATED runs own no annotation rows. """
    new_max_variant_id = absorbed.max_variant_id
    new_count = (survivor.count or 0) + (absorbed.count or 0)
    with transaction.atomic():
        locked_survivor = AnnotationRangeLock.objects.select_for_update().get(pk=survivor.pk)
        locked_absorbed = AnnotationRangeLock.objects.select_for_update().get(pk=absorbed.pk)
        # Belt-and-braces lock of the runs being cascade-deleted (against manual/admin actions)
        list(AnnotationRun.objects.select_for_update().filter(
            annotation_range_lock_id__in=[survivor.pk, absorbed.pk]))
        locked_survivor.max_variant_id = new_max_variant_id
        locked_survivor.count = new_count
        locked_survivor.save()
        locked_absorbed.delete()  # cascades its CREATED runs
    survivor.max_variant_id = new_max_variant_id
    survivor.count = new_count


def merge_pending_range_locks(variant_annotation_version: VariantAnnotationVersion, batch_max=None) -> int:
    """ #2667: Greedily combine consecutive pending range locks into larger ones (capped at batch_max)
        so a freed worker picks up one efficient merged batch instead of many tiny drip-runs.

        Inverse of subdivide_annotation_range_lock. Only operates on locks whose runs are all
        dispatchable; "adjacent" = next lock by min_variant (locks tile the space in increasing-pk
        order). Returns the number of locks absorbed (for logging/tests). """
    if batch_max is None:
        batch_max = sys.maxsize
    now = timezone.now()
    locks = list(AnnotationRangeLock.objects.filter(version=variant_annotation_version)
                 .order_by("min_variant_id"))

    total_absorbed = 0
    i = 0
    while i < len(locks):
        survivor = locks[i]
        if not _range_lock_is_dispatchable(survivor, now):
            i += 1
            continue
        j = i + 1
        while j < len(locks):
            candidate = locks[j]
            if not _range_lock_is_dispatchable(candidate, now):
                break
            if (survivor.count or 0) + (candidate.count or 0) > batch_max:
                break
            _absorb_range_lock(survivor, candidate)
            total_absorbed += 1
            j += 1
        i = j

    if total_absorbed:
        logging.info("merge_pending_range_locks(%s): absorbed %d range lock(s)",
                     variant_annotation_version, total_absorbed)
    return total_absorbed


def get_lowest_unannotated_variant_id(variant_annotation_version):
    # Get min_variant_id from annotation lock that hasn't completed
    # There can be multiple AnnotationRuns (Statndard/SV) for a range lock
    qs = AnnotationRun.objects.filter(annotation_range_lock__version=variant_annotation_version)
    unannotated_qs = qs.exclude(status=AnnotationStatus.FINISHED)
    data = unannotated_qs.aggregate(first_unannotated_variant_id=Min("annotation_range_lock__min_variant_id"))
    first_unannotated_variant_id = data["first_unannotated_variant_id"]
    if first_unannotated_variant_id is None:
        # All annotation locks completed - get 1 past the highest max
        annotated_qs = AnnotationRangeLock.objects.filter(version=variant_annotation_version,
                                                          annotationrun__status=AnnotationStatus.FINISHED)
        data = annotated_qs.aggregate(max_annotated_variant_id=Max("max_variant_id"))
        max_annotated_variant_id = data["max_annotated_variant_id"] or 0
        first_unannotated_variant_id = max_annotated_variant_id + 1

    return first_unannotated_variant_id
