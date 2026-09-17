import logging
import os
import threading
from datetime import timedelta

import celery
from celery import chain
from celery.canvas import Signature
from django.conf import settings
from django.db import connection, transaction
from django.db.models import F
from django.utils import timezone

from annotation.models import AnnotationStatus
from annotation.models.models import AnnotationRun, InvalidAnnotationVersionError
from annotation.pipelines import get_runner
from annotation.signals.manual_signals import (
    annotation_run_complete_signal,
    annotation_run_discarded_signal,
)
from eventlog.models import create_event
from library.enums.log_level import LogLevel
from library.log_utils import get_traceback, log_traceback, report_message

# #2667: kick the single-authority dispatcher by name to avoid importing annotation_scheduler_task
# (which imports this module). Mirror of analysis _trigger_rescheduling (#346).
DISPATCH_ANNOTATION_RUNS_TASK = "annotation.tasks.annotation_scheduler_task.dispatch_annotation_runs"


def _record_lock_failure(annotation_run_id, task_id: str, msg: str):
    """ #1658: record a failed task_id lock grab without touching any other column.

        Reaching this means another task holds the lock - we have just proved we do NOT own this run. A
        full save() of our instance here writes it whole over the winner's row, and the instance was loaded
        *before* the grab, so its task_id is the None that made us try in the first place: we would blank
        the winner's execution lock. The run then looks dispatchable (null task_id + dead lease) while the
        winner is still mid-VEP, so a third attempt can launch, and the winner's own lease heartbeat sees
        the task_id mismatch and aborts it. Same cascade as a reclaim, but no reclaim needed - a
        double-dispatch is enough, and the winner may be only seconds in.

        Re-read under a row lock rather than reusing our stale copy, so a set_task_log the winner commits
        concurrently isn't dropped on the way through this read-modify-write. """
    with transaction.atomic():
        annotation_run = AnnotationRun.objects.select_for_update().filter(pk=annotation_run_id).first()
        if annotation_run is None:
            return  # deleted under us - nothing to record against
        celery_task_logs = annotation_run.celery_task_logs or {}
        celery_task_logs.setdefault(task_id, {})["lock_failed"] = msg
        AnnotationRun.objects.filter(pk=annotation_run_id).update(celery_task_logs=celery_task_logs)


class AnnotationRunReclaimedError(Exception):
    """ #1658: this attempt lost its run - the lease expired, the dispatcher reclaimed the run and handed
        it to a fresh attempt, and this one was aborted (or finished) too late to own the result. Celery
        should still record the attempt as failed, but a distinct type keeps "did VEP break?" separable
        from "did we lose the race?" in logs and in any future retry logic. """


def _trigger_dispatch(variant_annotation_version_id):
    """ A run just completed (success OR failure) - a worker has freed up. Kick the dispatcher to
        launch the next (merged) batch. Both kicks serialise through scheduling_single_worker and
        fast-exit if nothing is dispatchable, so firing twice is safe (the dispatcher is single
        authority and leases atomically).

        Why two kicks: race condition releasing the lock. This worker has just cleared the run's
        task_id/lease and saved, but the dispatcher (a separate worker) computes free capacity by
        reading in-flight runs from the DB. If its read transaction starts before our completion
        commit is visible to it, it still counts this run as in-flight, sees no free slot, and exits -
        leaving the freed slot idle until the next event. The immediate kick handles the common case
        (commit already visible, lowest latency); the 3s-delayed kick re-runs the dispatcher once the
        lock release is unambiguously committed, so the freed capacity is actually picked up. """
    sig = Signature(DISPATCH_ANNOTATION_RUNS_TASK, args=(variant_annotation_version_id,))
    sig.apply_async()
    sig.apply_async(countdown=3)


def _dispatch_trigger_sig(variant_annotation_version_id) -> Signature:
    """ Immutable by-name dispatcher kick, for chaining after a retry's cleanup task so the dispatcher
        only re-picks the run once the cleanup has committed (a retry must not race its own scrub). The
        run then routes by status - upload-only stays ANNOTATION_COMPLETED (import lane), a full retry is
        a fresh CREATED run (VEP lane). #1649 """
    return Signature(DISPATCH_ANNOTATION_RUNS_TASK, args=(variant_annotation_version_id,), immutable=True)


class AnnotationRunLeaseHeartbeat:
    """ Renews an AnnotationRun's dispatcher lease on a background thread while the task actively runs,
        so a long-but-healthy run (e.g. a multi-hour structural-variant VEP dump, which has no SV-count
        cap) is never mistaken for a dead worker and reclaimed into a duplicate. Because a live run keeps
        its own lease fresh, ANNOTATION_RUN_LEASE_SECONDS can stay short - a genuinely dead worker stops
        heartbeating and is reclaimed within one window.

        The renew is guarded on our own task_id: if the run was reclaimed/reassigned out from under us the
        update matches no rows, so we never stomp a lease another worker now holds - and we stop. Each
        successful renew also syncs annotation_run.lease_expires in memory so the pipeline's own full
        save() calls can't write back the stale load-time value over the fresh lease. """

    def __init__(self, annotation_run: 'AnnotationRun', task_id: str):
        self.annotation_run = annotation_run
        self.task_id = task_id
        self.lease_seconds = settings.ANNOTATION_RUN_LEASE_SECONDS
        self.interval = settings.ANNOTATION_RUN_LEASE_HEARTBEAT_SECONDS
        self._stop = threading.Event()
        self._thread = None
        # #1658: handle to the currently-running subprocess (VEP), registered by execute_cmd on the task's
        # main thread and read by this heartbeat thread to abort it if the run is reclaimed under us. Guarded
        # by a lock because the two threads touch it concurrently.
        self._process = None
        self._process_lock = threading.Lock()

    def set_process(self, process):
        """ execute_cmd process_callback (#1658): register the live subprocess so the heartbeat can kill it
            if the lease is lost. Pass None once the subprocess has returned to drop the stale handle. """
        with self._process_lock:
            self._process = process

    def _abort_process(self):
        """ #1658: kill the registered subprocess so its blocking communicate() on the main thread returns a
            non-zero code, tripping the pipeline's `return_code != 0` guard and failing this losing attempt
            cleanly - without writing results over the new owner's run. Best-effort; already-exited is a no-op. """
        with self._process_lock:
            process = self._process
        if process is None or process.poll() is not None:
            return
        logging.warning("Lease heartbeat: killing subprocess for reclaimed AnnotationRun %s (task %s)",
                        self.annotation_run.pk, self.task_id)
        try:
            process.kill()
        except Exception:
            log_traceback()

    def _renew(self):
        now = timezone.now()
        lease_expires = now + timedelta(seconds=self.lease_seconds)
        updated = AnnotationRun.objects.filter(
            pk=self.annotation_run.pk, task_id=self.task_id).update(lease_expires=lease_expires)
        if updated:
            # Keep the in-memory instance in step so a full save() in the pipeline doesn't clobber the
            # renewed lease with the stale value loaded before the heartbeat started.
            self.annotation_run.lease_expires = lease_expires
        else:
            logging.warning("Lease heartbeat: AnnotationRun %s no longer held by task %s - stopping",
                            self.annotation_run.pk, self.task_id)
            self._stop.set()
            # #1658: escalate loss into aborting the run - terminate VEP so it doesn't run to completion and
            # double-write DB state for a run this worker no longer owns.
            self._abort_process()

    def _run(self):
        try:
            while not self._stop.wait(self.interval):
                try:
                    self._renew()
                except Exception:
                    # A transient DB blip must not kill the heartbeat (or the run) - log and retry next tick.
                    log_traceback()
        finally:
            connection.close()  # release this thread's own DB connection

    def __enter__(self):
        self._thread = threading.Thread(target=self._run, daemon=True,
                                        name=f"lease-heartbeat-{self.annotation_run.pk}")
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval)
        return False


@celery.shared_task
def delete_annotation_run(annotation_run_id):
    try:
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        annotation_run.status = AnnotationStatus.DELETING
        annotation_run.save()
        annotation_run.delete()
    except:
        log_traceback()
        raise


@celery.shared_task
def delete_annotation_run_uploaded_data(annotation_run_id):
    """ Deletes related objects but not actual run (used for retry annotation upload) """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    annotation_run.delete_related_objects()


@celery.shared_task
def assign_range_lock_to_annotation_run(annotation_run_id, annotation_range_lock_id):
    # has 1-to-1 so any previous AnnotationRuns linking to AnnotationRangeLock must have been deleted
    AnnotationRun.objects.filter(pk=annotation_run_id).update(annotation_range_lock_id=annotation_range_lock_id)


@celery.shared_task
def reset_annotation_run_for_retry(annotation_run_id):
    """ #1654: clear a run's partially-imported annotation rows and reset it in place to a clean CREATED
        state (keeping its range lock). Runs on a worker because the row-clear can take a minute or two for
        a large range; the retry request just queues this ahead of the dispatch trigger. """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    annotation_run.reset_for_retry()


@celery.shared_task
def annotate_variants(annotation_run_id):
    """ VEP lane (annotation_workers): dump unannotated variants + run the VEP subprocess, then stop at
        ANNOTATION_COMPLETED and release the lease. The DB upload is deliberately NOT done here - the
        dispatcher re-picks the completed run in its resume lane and runs import_annotation_run on
        db_workers, so a throttled VEP slot is never held through the bulk insert. External runs
        (#1568) never reach this - they're dumped externally and rejoin post-VEP as ANNOTATION_COMPLETED.
        See #1649. """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    my_task_id = annotate_variants.request.id
    logging.info("annotate_variants: %s", annotation_run)

    # External annotation (#1568): VEP for these runs is managed externally via the annotation_external
    # command. Never auto-run VEP here while waiting for the operator to import an annotated VCF.
    if annotation_run.external and annotation_run.vcf_annotated_filename is None:
        logging.info("Skipping external AnnotationRun %s (awaiting external annotation)", annotation_run.pk)
        return

    # task_id used as Celery lock. attempt_count is bumped here - at execution, not dispatch - so a run
    # whose task merely sat in an overloaded queue past its lease never burns its retry budget.
    num_modified = AnnotationRun.objects.filter(pk=annotation_run.pk, task_id__isnull=True).update(
        task_id=my_task_id, attempt_count=F("attempt_count") + 1)
    if num_modified != 1:
        msg = f"Celery couldn't get task_id lock on AnnotationRun: {annotation_run.pk}"
        _record_lock_failure(annotation_run.pk, my_task_id, msg)
        raise ValueError(msg)

    try:
        # Reload to get updated task_id
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        if blocker := annotation_run.variant_annotation_version.get_annotation_run_blocker():
            # We need this so that transcript/versions are in DB so FKs link
            msg = f"{annotation_run.variant_annotation_version} {blocker}"
            raise InvalidAnnotationVersionError(msg)
        annotation_run.task_id = my_task_id
        annotation_run.set_task_log("start", timezone.now())
        annotation_run.save()

        # Renew the lease while the (potentially many-hour, e.g. structural-variant) VEP work runs, so a
        # live worker is never reclaimed under us into a duplicate run.
        with AnnotationRunLeaseHeartbeat(annotation_run, my_task_id) as lease_heartbeat:
            if annotation_run.vcf_annotated_filename is None:
                get_runner(annotation_run.pipeline_type).annotate(annotation_run,
                                                                  lease_heartbeat=lease_heartbeat)
        # DB upload now runs as a separate db_workers task (import_annotation_run), launched by the
        # dispatcher when this run reaches ANNOTATION_COMPLETED. #1649
    except Exception as e:
        tb = get_traceback()
        annotation_run.error_exception = tb
        annotation_run.set_task_log("error_exception", tb)
        # The `finally` below writes the whole instance again, so this save is really here for its return
        # value - the ownership answer. Deliberately a conditional UPDATE rather than a cheaper .exists()
        # check: the row count is authoritative at the instant it runs, where a read-then-decide has a
        # window. The extra UPDATE is not worth optimising away - annotation wall-clock is VEP, not locks.
        if not annotation_run.save_if_owner(my_task_id):
            # #1658: we lost the race - the lease expired, the run was reclaimed and a fresh attempt owns
            # it now (the exception above is typically our own heartbeat killing VEP). Report at warning
            # level: this is the lease design working as intended under worker churn, not a pipeline
            # failure, and raising a Rollbar error here describes a VEP break that did not happen - which
            # both cries wolf and buries the reclaims that genuinely are pathological.
            msg = f"AnnotationRun {annotation_run.pk} was reclaimed - task {my_task_id} no longer owns it"
            report_message(message=msg, level='warning', extra_data={'output': str(e), 'error': tb})
            create_event(None, "AnnotationRun reclaimed", tb, severity=LogLevel.WARNING)
            raise AnnotationRunReclaimedError(msg) from e

        error_message = f"{e}: {annotation_run.pipeline_stderr}"
        name = 'Annotation pipeline run ' + str(annotation_run.id)
        report_message(message=name,
                       level='error',
                       extra_data={'output': error_message,
                                   'error': tb})

        create_event(None, "AnnotationRun failed", tb, severity=LogLevel.ERROR)
        raise
    finally:
        annotation_run.set_task_log("end", timezone.now())
        # #2667: release the dispatcher lease so the now-completed run no longer looks in-flight,
        # then kick the dispatcher to refill this freed worker slot (covers success and failure).
        # #1658: conditional on still holding the run. A stalled attempt whose run was reclaimed mid-VEP
        # would otherwise blank the *new* owner's task_id/lease - re-opening the run to the dispatcher
        # while the new owner is still running VEP, and tripping the new owner's own heartbeat into
        # aborting itself. One reclaim cascades.
        annotation_run.task_id = None
        annotation_run.leased_by = None
        annotation_run.lease_expires = None
        if not annotation_run.save_if_owner(my_task_id):
            logging.warning("AnnotationRun %s no longer held by task %s - leaving the new owner's lease "
                            "alone and discarding this attempt's state", annotation_run.pk, my_task_id)
            _cleanup_reclaimed_run_files(annotation_run)
        # The dispatcher is single authority and fast-exits when nothing is dispatchable, so kicking it
        # from a losing attempt too is harmless.
        _trigger_dispatch(annotation_run.annotation_range_lock.version_id)


@celery.shared_task
def import_annotation_run(annotation_run_id):
    """ Import lane (db_workers): bulk-load an already-annotated VCF into the DB. Reached for every run
        once it is ANNOTATION_COMPLETED with an annotated VCF present - whether VEP just finished locally,
        an external run (#1568) was imported, or an upload-only retry reset it. Kept off annotation_workers
        so quick DB inserts never consume a throttled VEP slot. See #1649. """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    my_task_id = import_annotation_run.request.id

    # task_id used as Celery lock (as annotate_variants) so a double-dispatch can't double-import.
    # attempt_count bumped at execution, not dispatch (as annotate_variants).
    num_modified = AnnotationRun.objects.filter(pk=annotation_run.pk, task_id__isnull=True).update(
        task_id=my_task_id, attempt_count=F("attempt_count") + 1)
    if num_modified != 1:
        msg = f"Celery couldn't get task_id lock on AnnotationRun: {annotation_run.pk}"
        _record_lock_failure(annotation_run.pk, my_task_id, msg)
        raise ValueError(msg)

    try:
        # Reload to get updated task_id
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        # Re-check state at execution time (the is_upload_resumable predicate): a queued import task can
        # be stale - the run already imported via an earlier dispatch (and the cleanup receiver removed
        # its annotated file), or was reclaimed. A stale duplicate is a silent no-op, never an ERROR;
        # the finally below releases the lock.
        if annotation_run.status != AnnotationStatus.ANNOTATION_COMPLETED \
                or not (annotation_run.vcf_annotated_filename
                        and os.path.exists(annotation_run.vcf_annotated_filename)):
            logging.info("AnnotationRun %s no longer needs importing (status=%s) - stale queued task %s",
                         annotation_run.pk, annotation_run.status, my_task_id)
            return
        annotation_run.set_task_log("import_start", timezone.now())
        annotation_run.save()

        # Renew the lease while the bulk DB insert runs, so a live worker is never reclaimed under us.
        with AnnotationRunLeaseHeartbeat(annotation_run, my_task_id):
            get_runner(annotation_run.pipeline_type).import_results(annotation_run)
        # A completed upload outranks any stale error - reclaim may have failed the run out while this
        # task sat queued, and get_status() tests error_exception first.
        annotation_run.error_exception = None
        # The run is only truly complete once imported (moved here from annotate_variants). #1649
        # #1670: annotation_run is carried so the cleanup receiver can reclaim this run's output now that
        # its rows are committed. Existing receivers take **kwargs, so the extra key is transparent.
        annotation_run_complete_signal.send(sender=os.path.basename(__file__),
                                            variant_annotation_version=annotation_run.annotation_range_lock.version,
                                            pipeline_type=annotation_run.pipeline_type,
                                            annotation_run=annotation_run)
    except Exception as e:
        tb = get_traceback()
        annotation_run.error_exception = tb
        annotation_run.set_task_log("error_exception", tb)
        # As annotate_variants: a conditional UPDATE rather than a cheaper ownership read, so the answer
        # has no window. Cost is irrelevant next to the import it follows.
        if not annotation_run.save_if_owner(my_task_id):
            # #1658: reclaimed mid-import - the new attempt owns the run, so this attempt's traceback must
            # not flip it to ERROR. A lost race is worker churn, not a pipeline failure (see annotate_variants).
            msg = f"AnnotationRun {annotation_run.pk} was reclaimed - task {my_task_id} no longer owns it"
            create_event(None, "AnnotationRun import reclaimed", tb, severity=LogLevel.WARNING)
            raise AnnotationRunReclaimedError(msg) from e

        create_event(None, "AnnotationRun import failed", tb, severity=LogLevel.ERROR)
        raise
    finally:
        annotation_run.set_task_log("import_end", timezone.now())
        # #2667: release the dispatcher lease so the now-completed run no longer looks in-flight, then
        # kick the dispatcher to refill this freed import slot (covers success and failure).
        # #1658: conditional on still holding the run, so a reclaimed attempt can't blank the new owner's
        # lease and re-open the run to the dispatcher while that owner is still importing.
        annotation_run.task_id = None
        annotation_run.leased_by = None
        annotation_run.lease_expires = None
        if not annotation_run.save_if_owner(my_task_id):
            logging.warning("AnnotationRun %s no longer held by task %s - leaving the new owner's lease "
                            "alone and discarding this attempt's state", annotation_run.pk, my_task_id)
        # The import lane owns no files of its own - the annotated VCF belongs to the run, and the new
        # owner is importing from it - so there is nothing to clean up here.
        _trigger_dispatch(annotation_run.annotation_range_lock.version_id)


def _cleanup_reclaimed_run_files(annotation_run):
    """ #1658: discard the on-disk output of an attempt that lost its run.

        Reclaim (_reset_run_for_redispatch) deletes what it can derive from the dump stem and then NULLs
        the filename fields, so by the time the losing attempt unwinds the DB no longer references its
        per-task files. There is no sweep of ANNOTATION_VCF_DUMP_DIR, so anything reclaim could not name
        has to be collected here. Before per-task paths (#1658) the annotated name was a fixed function of
        the run, so a later attempt overwrote it or a later reclaim deleted it by name; now this attempt
        is the last party still holding its own paths.

        Passes the in-memory instance deliberately: it still carries the filenames reclaim has since
        NULLed on the row. #1670: the removal itself lives in annotation.signals.annotation_run_cleanup. """
    annotation_run_discarded_signal.send(sender=os.path.basename(__file__), annotation_run=annotation_run)


def annotation_run_retry(annotation_run: AnnotationRun, upload_only=False) -> AnnotationRun:
    if upload_only:
        if annotation_run.vcf_annotated_filename is None:
            msg = "Retry annotation run upload only requires annotation VCF to be written"
            raise ValueError(msg)
        if not os.path.exists(annotation_run.vcf_annotated_filename):
            # #1670: a run that imported successfully has had its VEP output reclaimed, so the row still
            # names a file that is gone. Say so here rather than failing deep inside cyvcf2.
            msg = (f"Annotation VCF '{annotation_run.vcf_annotated_filename}' no longer exists "
                   f"(removed after successful annotation) - retry the full run instead of upload only")
            raise ValueError(msg)

    annotation_range_lock = annotation_run.annotation_range_lock
    if annotation_range_lock is None:
        msg = "Can't retry annotation run with no annotation lock!"
        raise ValueError(msg)

    if upload_only:
        # A manual retry is a fresh start - reset the lease bookkeeping so the dispatcher's attempt-cap
        # (reclaim_stalled_annotation_runs) counts from zero again. Otherwise an upload-only retry reuses
        # the same run, whose attempt_count is already at ANNOTATION_MAX_RUN_ATTEMPTS, so the next stall
        # fails it immediately with "exceeded max attempts".
        annotation_run.error_exception = None  # Clear so status won't be error
        annotation_run.task_id = None  # Allow celery jobs to get lock on it
        annotation_run.leased_by = None
        annotation_run.lease_expires = None
        annotation_run.attempt_count = 0
        annotation_run.save()
        # Delete uploaded data, then hand back to the dispatcher: the run stays ANNOTATION_COMPLETED with
        # its annotated VCF present -> import lane -> db_workers (import_annotation_run). #1649
        tasks = [
            delete_annotation_run_uploaded_data.si(annotation_run.pk),
        ]
    else:
        # #1654: full retry resets the run in place to a clean CREATED state (clearing any partially-
        # imported annotation rows) rather than delete + recreate. Reusing the same run - keeping its
        # range lock - means no rangeless AnnotationRun is ever committed, so a mid-retry crash can no
        # longer strand it invisibly in Created. The row-clear can take a minute or two for a large
        # range, so it runs on a worker (ahead of the dispatch trigger) instead of blocking this request;
        # the run stays ERROR (visible, retryable) until the reset flips it to CREATED. #1649
        tasks = [
            reset_annotation_run_for_retry.si(annotation_run.pk),
        ]

    # #1649: retry no longer launches annotate_variants inline - it hands back to the single-authority
    # dispatcher (preserving lease/reclaim/attempt-cap), which launches the right lane by status.
    tasks.append(_dispatch_trigger_sig(annotation_range_lock.version_id))
    task = chain(tasks)
    task.apply_async()

    return annotation_run


