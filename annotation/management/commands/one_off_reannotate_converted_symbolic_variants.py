"""
Re-annotates, through the STRUCTURAL_VARIANT pipeline, variants whose stored annotation came from the STANDARD
pipeline although their version now routes them to STRUCTURAL_VARIANT (#2786). These are the explicit del/dup/inv
that one_off_fix_symbolic_variants converted in place to a symbolic <DEL>/<DUP>/<INV> at or over the version's
structural_variant_min_size: the conversion kept the pk, and so the annotation VEP wrote for the explicit
sequence, and the scheduler only ever fills new pk ranges.

Only the versions the dispatcher annotates (the latest NEW and ACTIVE per build) are repaired. The variants are
grouped by the range lock whose STANDARD run annotated them; each group's STANDARD rows (VariantAnnotation,
VariantTranscriptAnnotation, VariantGeneOverlap) are deleted, and a new AnnotationRangeLock over the group's pk
bounds - inside the original lock - gets the runs the scheduler would give it. The VEP pipelines take only the
unannotated variants of their type in range (AnnotSV takes every SV in range, and only UPDATEs), so the new lock's
STRUCTURAL_VARIANT run picks up exactly the deleted variants, and the dispatcher annotates them like any other
backlog; a failed run is retried from the annotation runs page. A group whose range still has an unfinished VEP
run is left alone, for a later run of this command, which finds only what is still left.

Other versions are counted but not touched: the dispatcher never runs a HISTORICAL version's annotation, and its
VEP no longer matches the installed one. --dry-run prints the counts and changes nothing.
Registered as a deploy-time ManualOperation (annotation/migrations/0017_reannotate_converted_symbolic_variants.py).
"""
from dataclasses import dataclass, field

from django.core.management.base import BaseCommand
from django.db import transaction

from annotation.annotation_version_querysets import (
    get_queryset_for_annotation_version,
    get_variants_qs_for_annotation,
)
from annotation.annotation_versions import _reset_run_counts_after_extend
from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.models.models import (
    VariantAnnotation,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from annotation.tasks.annotation_scheduler_task import (
    _dispatchable_variant_annotation_versions,
    _handle_range_lock,
    _scheduled_pipeline_versions,
    dispatch_annotation_runs,
)
from library.django_utils.data_archive_mixin import DataArchivedError

STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL_VARIANT = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


@dataclass
class LockGroup:
    """ Misrouted variants annotated by one range lock's STANDARD run """
    range_lock_id: int
    range_lock: AnnotationRangeLock = None
    variant_ids: list[int] = field(default_factory=list)
    standard_run_ids: set[int] = field(default_factory=set)


def find_misrouted_variants(annotation_version) -> list[LockGroup]:
    """ Variants the version routes to STRUCTURAL_VARIANT that hold a STANDARD run's annotation """
    qs = get_variants_qs_for_annotation(annotation_version, pipeline_type=STRUCTURAL_VARIANT, annotated=True)
    qs = qs.filter(variantannotation__annotation_run__pipeline_type=STANDARD)
    rows = qs.order_by("pk").values_list("pk", "variantannotation__annotation_run_id",
                                         "variantannotation__annotation_run__annotation_range_lock_id")
    groups: dict[int, LockGroup] = {}
    for variant_id, run_id, range_lock_id in rows:
        group = groups.setdefault(range_lock_id, LockGroup(range_lock_id))
        group.variant_ids.append(variant_id)
        group.standard_run_ids.add(run_id)
    range_locks = AnnotationRangeLock.objects.in_bulk(groups)
    for range_lock_id, group in groups.items():
        group.range_lock = range_locks[range_lock_id]
    return list(groups.values())


def has_unfinished_vep_runs(vav: VariantAnnotationVersion, min_variant_id: int, max_variant_id: int) -> bool:
    """ A VEP run that may still write or retry rows over this range would race the new lock's run """
    qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav,
                                      annotation_range_lock__min_variant_id__lte=max_variant_id,
                                      annotation_range_lock__max_variant_id__gte=min_variant_id,
                                      pipeline_type__in=[STANDARD, STRUCTURAL_VARIANT])
    return qs.exclude(status=AnnotationStatus.FINISHED).exists()


def delete_standard_annotation(annotation_version, group: LockGroup):
    for klass in [VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap]:
        qs = get_queryset_for_annotation_version(klass, annotation_version)
        qs.filter(variant_id__in=group.variant_ids, annotation_run_id__in=group.standard_run_ids).delete()


def create_range_locks(vav: VariantAnnotationVersion, group: LockGroup, pipeline_versions) -> int:
    """ Lock(s) over the group's bounds with the scheduler's runs. Returns the number of locks created """
    variant_ids = group.variant_ids
    original = group.range_lock
    min_variant_id, max_variant_id = variant_ids[0], variant_ids[-1]
    if (min_variant_id, max_variant_id) != (original.min_variant_id, original.max_variant_id):
        bounds = [(min_variant_id, max_variant_id)]
    elif min_variant_id != max_variant_id:
        # A second lock over the original's bounds breaks unique_together, so the first variant gets its own
        bounds = [(min_variant_id, min_variant_id), (variant_ids[1], max_variant_id)]
    else:
        # The original lock holds only this variant, so its STRUCTURAL_VARIANT run finished empty: reopen it
        _reset_run_counts_after_extend(original)
        return 0

    for lock_min, lock_max in bounds:
        count = sum(lock_min <= variant_id <= lock_max for variant_id in variant_ids)
        range_lock = AnnotationRangeLock.objects.create(version=vav, min_variant_id=lock_min,
                                                        max_variant_id=lock_max, count=count)
        _handle_range_lock(range_lock, pipeline_versions)
    return len(bounds)


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true", help="Print the counts without changing anything")
        parser.add_argument("--variant-annotation-version", type=int,
                            help="Only this version (default: every version whose data isn't archived)")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        dispatchable_ids = {vav.pk for vav in _dispatchable_variant_annotation_versions()}
        qs = VariantAnnotationVersion.objects.filter(data_archived_date__isnull=True).select_related("genome_build")
        if vav_id := options["variant_annotation_version"]:
            qs = qs.filter(pk=vav_id)
        for vav in qs.order_by("pk"):
            prefix = "DRY RUN " if dry_run else ""
            self.stdout.write(f"{prefix}{vav}: {self._handle_version(vav, vav.pk in dispatchable_ids, dry_run)}")

    @staticmethod
    def _handle_version(vav: VariantAnnotationVersion, dispatchable: bool, dry_run: bool) -> str:
        annotation_version = vav.get_any_annotation_version()
        if annotation_version is None:
            return "skipped: no AnnotationVersion"
        try:
            groups = find_misrouted_variants(annotation_version)
        except DataArchivedError as e:
            return f"skipped: {e}"
        num_variants = sum(len(g.variant_ids) for g in groups)
        summary = f"{num_variants} variants with STANDARD annotation now routed to STRUCTURAL_VARIANT"
        if not num_variants:
            return summary
        if not dispatchable:
            return f"{summary} - left as they are: the dispatcher doesn't annotate this {vav.get_status_display()} version"
        if blocker := vav.get_annotation_run_blocker():
            return f"{summary} - left as they are until annotation is unblocked: {blocker}"

        ready = []
        waiting = []
        for group in groups:
            if has_unfinished_vep_runs(vav, group.variant_ids[0], group.variant_ids[-1]):
                waiting.append(group)
            else:
                ready.append(group)
        num_ready = sum(len(g.variant_ids) for g in ready)
        if waiting:
            num_waiting = sum(len(g.variant_ids) for g in waiting)
            summary += f" ({num_waiting} wait for unfinished annotation runs over their range - run this again later)"
        if dry_run or not ready:
            return f"{summary}, {num_ready} to re-annotate from {len(ready)} range locks"

        pipeline_versions = _scheduled_pipeline_versions(vav.genome_build)
        num_locks = 0
        with transaction.atomic():
            for group in ready:
                delete_standard_annotation(annotation_version, group)
                num_locks += create_range_locks(vav, group, pipeline_versions)
        dispatch_annotation_runs.si(vav.pk).apply_async()
        return (f"{summary}, {num_ready} STANDARD annotations deleted and queued for STRUCTURAL_VARIANT annotation "
                f"({num_locks} new range locks)")
