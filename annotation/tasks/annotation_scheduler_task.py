import logging
import os
from datetime import timedelta

import celery
from django.conf import settings
from django.core.cache import cache
from django.db.models import Case, F, IntegerField, Q, Value, When
from django.utils import timezone

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count, \
    merge_pending_range_locks
from annotation.celery_utils import annotation_worker_slots
from annotation.models import AnnotationRun, AnnotationStatus, VariantAnnotationPipelineType, \
    VariantAnnotationVersion
from annotation.models.models import AnnotationVersion, AnnotationRangeLock
from annotation.tasks.annotate_variants import annotate_variants
from library.log_utils import log_traceback
from snpdb.models import GenomeBuild, ImportStatus, Sample, VCF, Variant, JobsControl

# #1646: max count_annotation_run tasks the dispatcher kicks per cycle, so creating a new annotation
# version (thousands of runs at once) doesn't burst the count queue - the rest drain over later cycles.
COUNT_KICK_BATCH = 500


@celery.shared_task(queue='scheduling_single_worker')
def annotation_scheduler(status: str = None):
    """ Run on scheduling_single_worker queue to avoid race conditions.
        `status` is a VariantAnnotationVersion.Status value (default ACTIVE).
        Never operates on HISTORICAL VAVs. """
    if status is None:
        status = VariantAnnotationVersion.Status.ACTIVE
    if status == VariantAnnotationVersion.Status.HISTORICAL:
        raise ValueError("annotation_scheduler must not be run against HISTORICAL VariantAnnotationVersions")

    LOCK_EXPIRE = 60 * 5  # 5 minutes
    lock_id = "annotation-scheduler-lock"

    # cache.add fails if the key already exists
    acquire_lock = lambda: cache.add(lock_id, "true", LOCK_EXPIRE)
    release_lock = lambda: cache.delete(lock_id)

    try:
        if acquire_lock():
            try:
                logging.info("Got the lock for annotation scheduler (status=%s)", status)
                for genome_build in GenomeBuild.builds_with_annotation():
                    annotation_version = AnnotationVersion.latest(genome_build, status=status, validate=False)
                    if annotation_version is None:
                        continue
                    variant_annotation_version = annotation_version.variant_annotation_version
                    while True:
                        range_lock = _handle_variant_annotation_version(variant_annotation_version)
                        if range_lock is None:
                            break
                    # #2667: scheduler only creates pending state - the dispatcher decides what launches
                    # (capacity-limited, merging the backlog first). Kick it for this VAV. The dispatcher
                    # also fires the #1646 count tasks, so it is the single place that keeps counts moving.
                    dispatch_annotation_runs.si(variant_annotation_version.pk).apply_async()
            finally:
                logging.info("Releasing lock")
                release_lock()
        else:
            logging.info("Someone else has %s", lock_id)
    except:
        log_traceback()


def _handle_range_lock(range_lock, pipeline_type=None):
    pipeline_types = []
    if pipeline_type is not None:
        pipeline_types.append(pipeline_type)
    else:
        pipeline_types = list(VariantAnnotationPipelineType)

    for pipeline_type in pipeline_types:
        annotation_run, _ = AnnotationRun.objects.get_or_create(annotation_range_lock=range_lock,
                                                                pipeline_type=pipeline_type)
        if annotation_run.external:
            # External annotation (#1568): VEP is managed off-VM via the annotation_external command.
            # Belt-and-braces guard - the scheduler only operates on ACTIVE versions and external runs
            # live on NEW versions, so this should not normally be reached.
            logging.info("Skipping external AnnotationRun %s (awaiting external annotation)", annotation_run.pk)
            continue
        # #2667: the run is left pending (CREATED, un-leased). dispatch_annotation_runs launches it
        # when there is worker capacity, merging the backlog into bigger batches first.


@celery.shared_task(queue='db_workers')
def count_annotation_run(annotation_run_id):
    """ count variants a pending run will process and store in count
        If empty, finish without requiring a full dump/dispatch

        We don't do this on scheduling queue as that has to be light/fast """
    annotation_run = AnnotationRun.objects.filter(pk=annotation_run_id).first()
    if annotation_run is None or annotation_run.external:
        return
    if annotation_run.status != AnnotationStatus.CREATED or annotation_run.count is not None:
        return  # already counted, or already leased/launched/finished
    range_lock = annotation_run.annotation_range_lock
    if range_lock is None:
        # Something else cleared this - this job is just a way to have a quick exit skipping is fine
        return
    annotation_version = range_lock.version.get_any_annotation_version()
    count = get_variants_qs_for_annotation(annotation_version, pipeline_type=annotation_run.pipeline_type,
                                           min_variant_id=range_lock.min_variant_id,
                                           max_variant_id=range_lock.max_variant_id).count()
    update = {"count": count}
    if count == 0:
        # Empty - finish without a dump. dump_count=0 keeps get_status() consistent (-> FINISHED).
        now = timezone.now()
        update.update(status=AnnotationStatus.FINISHED, dump_count=0, dump_start=now, dump_end=now)
    updated = AnnotationRun.objects.filter(
        pk=annotation_run_id,
        status=AnnotationStatus.CREATED,
        task_id__isnull=True,
        lease_expires__isnull=True,
        count__isnull=True,
    ).update(**update)
    if updated and count == 0:
        logging.info("Finished empty AnnotationRun %s (%s) via count task", annotation_run_id,
                     annotation_run.pipeline_type)


def _trigger_counts_for_uncounted_runs(vav: VariantAnnotationVersion):
    """ #1646: kick count_annotation_run (db_workers) for pending (CREATED, un-leased) runs of this VAV
        with no `count` yet - freshly created runs, and runs a merge reopened. Lowest-pk-first (counted
        before they'd be dispatched) and capped per call so a fresh backlog of thousands doesn't burst
        the queue in one go - the remainder is picked up on following dispatch cycles. Idempotent: the
        count task's guarded update no-ops on a run that got counted or picked up in the meantime, and
        the count__isnull filter drops runs once counted, so re-kicking does not build up. """
    run_ids = AnnotationRun.objects.filter(
        annotation_range_lock__version=vav,
        external=False,
        status=AnnotationStatus.CREATED,
        task_id__isnull=True,
        count__isnull=True,
    ).order_by("annotation_range_lock__min_variant_id").values_list("pk", flat=True)[:COUNT_KICK_BATCH]
    for run_id in run_ids:
        count_annotation_run.apply_async((run_id,))  # @UndefinedVariable


def _handle_variant_annotation_version(variant_annotation_version):
    # If we crash in the wrong place, we may end up with unassigned range locks (need 1 of each type)
    arl_qs = AnnotationRangeLock.objects.filter(version=variant_annotation_version)
    for pipeline_type in VariantAnnotationPipelineType:
        # Look for missing any AnnotationRun for this lock
        for range_lock in arl_qs.exclude(annotationrun__pipeline_type=pipeline_type):
            logging.warning("Assigned orphaned annotation range lock: %s, run type: %s",
                            range_lock, VariantAnnotationPipelineType(pipeline_type).label)
            _handle_range_lock(range_lock, pipeline_type)

    range_lock, unannotated_count = get_annotation_range_lock_and_unannotated_count(variant_annotation_version,
                                                                                    settings.ANNOTATION_VEP_BATCH_MIN,
                                                                                    settings.ANNOTATION_VEP_BATCH_MAX)
    if range_lock is not None:
        range_lock.save()
        _handle_range_lock(range_lock)
    else:
        if unannotated_count:
            logging.warning("Unannotated variants but couldn't get a lock.")
            logging.warning("If you are sure no annotation requests are running you can reset things with 'reset_annotation_states'")
        else:
            logging.info("No unannotated variants left!")

            waiting_vcfs = VCF.objects.filter(import_status=ImportStatus.IMPORTING)
            waiting_samples = Sample.objects.filter(import_status=ImportStatus.IMPORTING)
            if waiting_vcfs or waiting_samples:
                logging.warning("Waiting vcfs or samples - if you are sure no annotation requests are running you can set vcfs to be available with 'annotation_set_all_complete' ")
            else:
                logging.info("No waiting vcfs or samples (this is good!)")

    return range_lock


@celery.shared_task(queue='scheduling_single_worker')
def dispatch_annotation_runs(variant_annotation_version_id=None):
    """ #2667: Single-authority dispatcher (mirrors analysis create_and_launch_analysis_tasks). The
        scheduler only creates pending AnnotationRuns; this decides what actually launches.

        Latency-first: while workers are free, pending work launches immediately as-is. Merging only
        kicks in once all workers are busy and a backlog has formed - then a freed worker picks up one
        efficient merged batch instead of many tiny runs.

        With no id, sweeps every dispatchable VAV - the latest NEW and latest ACTIVE per build (see
        _dispatchable_variant_annotation_versions). All lease writes happen here on
        scheduling_single_worker so they serialise - no row-lock gymnastics. """
    if JobsControl.is_paused():
        return  # operational brake (e.g. crash safety auto-pause) - don't launch annotation runs
    try:
        if variant_annotation_version_id is not None:
            vav = VariantAnnotationVersion.objects.get(pk=variant_annotation_version_id)
            _dispatch_for_vav(vav)
        else:
            for vav in _dispatchable_variant_annotation_versions():
                _dispatch_for_vav(vav)
    except:
        log_traceback()


def _dispatchable_variant_annotation_versions() -> list[VariantAnnotationVersion]:
    """ Versions the no-arg beat sweep dispatches: the latest NEW and latest ACTIVE per build (never
        HISTORICAL). NEW is included because that is where a new annotation version is actively being
        built - it carries the bulk of the pending runs (including the external upload-only runs imported
        onto it), so without a heartbeat here it only ever advances on discrete kicks (VCF import /
        manual / per-completion re-kick) and a stalled dispatch never self-recovers (stale leases go
        unreclaimed, freed capacity unused). Safe to include: the dispatcher already filters external=False
        and annotate_variants no-ops external runs awaiting an off-VM VCF, and get_annotation_run_blocker()
        still guards a not-yet-ready version. """
    statuses = [VariantAnnotationVersion.Status.NEW, VariantAnnotationVersion.Status.ACTIVE]
    vavs = []
    seen = set()
    for genome_build in GenomeBuild.builds_with_annotation():
        for status in statuses:
            annotation_version = AnnotationVersion.latest(genome_build, status=status, validate=False)
            if annotation_version is not None:
                vav = annotation_version.variant_annotation_version
                if vav.pk not in seen:
                    seen.add(vav.pk)
                    vavs.append(vav)
    return vavs


def _dispatch_for_vav(vav: VariantAnnotationVersion):
    now = timezone.now()
    reclaim_stalled_annotation_runs(vav, now)
    # #1646: fill `count` (and finish empties) off-thread before we consider dispatching. Fired every
    # cycle - this is the reliable driver, since the scheduler only runs on VCF import / manual trigger
    # while the dispatcher runs on every completion. Runs even when capacity is 0 so a saturated
    # system still closes empty runs (on db_workers, not competing with VEP).
    _trigger_counts_for_uncounted_runs(vav)

    slots = annotation_worker_slots()
    in_flight = _in_flight_runs_qs(vav, now).count()
    capacity = slots - in_flight
    logging.info("dispatch_annotation_runs(%s): slots=%d in_flight=%d capacity=%d",
                 vav, slots, in_flight, capacity)
    if capacity <= 0:
        return

    if _dispatchable_runs_qs(vav, now).count() > capacity:
        # A real backlog exists (more pending work than free capacity) - merge it into bigger batches
        # before leasing. With free capacity for everything, we skip this and launch as-is (low latency).
        merge_pending_range_locks(vav, settings.ANNOTATION_VEP_BATCH_MAX)

    # #1646: launch the cheap jobs first so slots free up faster. Two tiers:
    #   0 - resume upload-only runs (past VEP, just the quick DB upload)
    #   1 - the heavy CREATED dump+VEP+upload runs
    # Empty runs never reach here - the count task (stage 2) finishes them without a dispatch. Within
    # each tier, lowest-pk-first so the 'lowest unannotated id' watermark keeps moving and waiting VCF
    # pipelines unblock in order.
    worker_id = f"dispatch:{dispatch_annotation_runs.request.id or 'sync'}"
    dispatchable = _dispatchable_runs_qs(vav, now).annotate(
        dispatch_priority=Case(When(status=AnnotationStatus.ANNOTATION_COMPLETED, then=Value(0)),
                               default=Value(1), output_field=IntegerField())
    ).order_by("dispatch_priority", "annotation_range_lock__min_variant_id")
    launched = 0
    for annotation_run in dispatchable:
        if launched >= capacity:
            break
        _lease_and_launch_run(annotation_run, worker_id, now)
        launched += 1
    logging.info("dispatch_annotation_runs(%s): launched %d run(s)", vav, launched)


def _dispatchable_runs_qs(vav: VariantAnnotationVersion, now):
    """ Runs ready to lease + launch, un-leased, no task_id, not external. Two lanes (#1646):
        pending (CREATED) runs get the full dump+VEP+upload; resume-upload runs (past VEP - annotation
        complete with an annotated VCF present) get re-launched upload-only. """
    created = Q(status=AnnotationStatus.CREATED)
    resume = Q(status=AnnotationStatus.ANNOTATION_COMPLETED, vcf_annotated_filename__isnull=False)
    return AnnotationRun.objects.filter(
        annotation_range_lock__version=vav,
        external=False,
        task_id__isnull=True,
    ).filter(created | resume).filter(Q(lease_expires__isnull=True) | Q(lease_expires__lt=now))


def _in_flight_runs_qs(vav: VariantAnnotationVersion, now):
    """ Runs occupying a worker slot: live lease, holding task_id, or a running status. A completed
        run (FINISHED/ERROR) never occupies a slot, even if its lease hasn't been cleared yet - so it
        is excluded regardless, otherwise stale leases on finished runs would shrink capacity.

        #1646: a reclaimed resume-upload run (ANNOTATION_COMPLETED, annotated VCF present, un-leased,
        no task_id) is waiting for the dispatcher, not running - exclude it so it doesn't shrink
        capacity while also appearing in _dispatchable_runs_qs. """
    completed = AnnotationStatus.get_completed_states()
    settled = [AnnotationStatus.CREATED] + list(completed)
    live_lease = Q(lease_expires__gte=now)
    has_task = Q(task_id__isnull=False)
    running = ~Q(status__in=settled)
    resume_ready = Q(status=AnnotationStatus.ANNOTATION_COMPLETED, vcf_annotated_filename__isnull=False,
                     task_id__isnull=True) & (Q(lease_expires__isnull=True) | Q(lease_expires__lt=now))
    return AnnotationRun.objects.filter(annotation_range_lock__version=vav, external=False) \
        .filter(live_lease | has_task | running) \
        .exclude(status__in=completed) \
        .exclude(resume_ready)


def _lease_and_launch_run(annotation_run: AnnotationRun, worker_id: str, now):
    lease_expires = now + timedelta(seconds=settings.ANNOTATION_RUN_LEASE_SECONDS)
    AnnotationRun.objects.filter(pk=annotation_run.pk).update(
        leased_by=worker_id,
        lease_expires=lease_expires,
        attempt_count=F("attempt_count") + 1,
    )
    annotate_variants.apply_async((annotation_run.pk,))  # @UndefinedVariable


def reclaim_stalled_annotation_runs(vav: VariantAnnotationVersion, now=None):
    """ #2667: Reclaim runs whose lease expired (worker died) back to dispatchable CREATED state.
        A run that has burned ANNOTATION_MAX_RUN_ATTEMPTS leases is failed to ERROR rather than
        retried forever. attempt_count is bumped at lease time (_lease_and_launch_run). """
    if now is None:
        now = timezone.now()
    expired_qs = AnnotationRun.objects.filter(
        annotation_range_lock__version=vav,
        external=False,
        lease_expires__lt=now,
    ).exclude(status__in=AnnotationStatus.get_completed_states())

    for annotation_run in expired_qs:
        if annotation_run.attempt_count >= settings.ANNOTATION_MAX_RUN_ATTEMPTS:
            logging.warning("Failing AnnotationRun %s - lease expired after %d attempts",
                            annotation_run.pk, annotation_run.attempt_count)
            _fail_stalled_run(annotation_run)
        else:
            logging.info("Reclaiming stalled AnnotationRun %s (attempt %d)",
                         annotation_run.pk, annotation_run.attempt_count)
            _reset_run_for_redispatch(annotation_run)


def _fail_stalled_run(annotation_run: AnnotationRun):
    annotation_run.error_exception = "AnnotationRun lease expired (worker lost); exceeded max attempts."
    annotation_run.leased_by = None
    annotation_run.lease_expires = None
    annotation_run.task_id = None
    annotation_run.save()  # get_status() -> ERROR (error_exception set)


def _reset_run_for_redispatch(annotation_run: AnnotationRun):
    """ Return a stalled run to a dispatchable state so the dispatcher can re-launch it. Always clears
        the lease/task lock. Two flavours of pipeline progress:

        #1646 resume-upload: the run is already past VEP (annotated VCF present on disk) and only the
        quick DB upload remains. Keep the expensive VEP output - just scrub any partially-imported rows
        and reset the upload markers, leaving status ANNOTATION_COMPLETED so the dispatcher re-launches
        it upload-only (annotate_variants skips dump+VEP straight to import_vcf_annotations). This avoids
        throwing away minutes of VEP per stalled run.

        Otherwise (pre-VEP progress): scrub all partial output (on-disk dumps + partial rows) back to a
        clean CREATED state so a re-dump won't collide. """
    annotated_filename = annotation_run.vcf_annotated_filename
    if annotated_filename and os.path.exists(annotated_filename):
        # Past VEP - resume upload-only. Scrub partially-imported rows; keep dump/VEP/AnnotSV artifacts.
        annotation_run.delete_related_objects()
        annotation_run.upload_start = None
        annotation_run.upload_end = None
        annotation_run.annotsv_imported = False  # re-import re-updates the recreated rows (idempotent)
    elif annotation_run.status != AnnotationStatus.CREATED:
        # Worker got part-way before dying (pre-VEP) - scrub partial output (a CREATED run has none).
        annotation_run.delete_related_objects()
        for filename in (annotation_run.vcf_dump_filename, annotation_run.vcf_annotated_filename):
            if filename and os.path.exists(filename):
                os.remove(filename)
        annotation_run.dump_start = None
        annotation_run.dump_end = None
        annotation_run.dump_count = None
        annotation_run.annotation_start = None
        annotation_run.annotation_end = None
        annotation_run.upload_start = None
        annotation_run.upload_end = None
        annotation_run.vcf_dump_filename = None
        annotation_run.vcf_annotated_filename = None
    annotation_run.task_id = None
    annotation_run.leased_by = None
    annotation_run.lease_expires = None
    annotation_run.error_exception = None
    annotation_run.save()  # get_status() -> CREATED, or ANNOTATION_COMPLETED for a resume-upload run


def subdivide_annotation_range_lock(arl: AnnotationRangeLock,
                                    minimum_size=AnnotationRangeLock.MIN_SIZE_FOR_SUBDIVISION) -> AnnotationRangeLock:
    """ Sometimes we have an annotation run crash or timeout etc, this splits it in half so we can try again
        With a smaller run

        Range lock passed in becomes bottom half of original range
        Returns new range lock for top half of original range
    """

    size = int(arl.max_variant_id) - int(arl.min_variant_id)
    if size < minimum_size:
        raise ValueError(f"Cannot subdivide {arl} below minimum size of {minimum_size}")

    logging.info("Subdividing %s", arl)

    # Delete all existing annotation runs
    res = arl.annotationrun_set.all().delete()
    logging.info("Deleted attached annotation runs: %s", res)

    half_size = size / 2
    halfway_point = int(arl.min_variant_id) + half_size
    first_at_or_above_halfway = Variant.objects.filter(pk__gte=halfway_point).order_by("pk").first()
    first_below_halfway = Variant.objects.filter(pk__gte=arl.min_variant_id, pk__lt=halfway_point).order_by("pk").last()

    av = arl.version.get_any_annotation_version()
    unannotated_qs = get_variants_qs_for_annotation(av, min_variant_id=first_at_or_above_halfway.pk,
                                                    max_variant_id=arl.max_variant.pk)

    top_half_count = unannotated_qs.count()
    new_lock = AnnotationRangeLock.objects.create(version=arl.version,
                                                  min_variant=first_at_or_above_halfway,
                                                  max_variant=arl.max_variant,
                                                  count=top_half_count)

    arl.max_variant = first_below_halfway
    arl.count -= top_half_count
    arl.save()

    logging.info("Shifted old lock down: %s", arl)
    logging.info("New lock created: %s", new_lock)
    return new_lock
