import logging
from datetime import timedelta

import celery
from auditlog.context import disable_auditlog
from celery.canvas import Signature
from celery.contrib.abortable import AbortableTask
from celery.result import AsyncResult
from django.db.models import OuterRef, Subquery, F, Q
from django.db.utils import OperationalError, IntegrityError
from django.utils import timezone

from analysis.exceptions import NodeConfigurationException, NodeParentErrorsException, CeleryTasksObsoleteException, \
    NodeOutOfDateException
from analysis.models.nodes.analysis_node import AnalysisNode, NodeStatus, NodeVersion, NodeCache, NodeTask, NodeColors
from eventlog.models import create_event
from library.constants import MINUTE_SECS
from library.enums.log_level import LogLevel
from library.log_utils import log_traceback, get_traceback
from snpdb.models import ProcessingStatus

CREATE_AND_LAUNCH_TASK = "analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks"

# Bounds both transient-error retries (_backoff_node) and dead-worker reclaims
# (analysis_update_tasks.lease_ready_nodes) with the single NodeTask.attempt_count.
MAX_NODE_ATTEMPTS = 3


def next_backoff(attempts):
    """ mocha next_backoff(), shorter ceiling. VG node loads are local DB work, not network
        calls, so transient failures are rare and clear fast: 10s, 20s, 40s ... cap 5 min. """
    return min(10 * (2 ** attempts), MINUTE_SECS * 5)


def _trigger_rescheduling(analysis_id):
    """ Event-driven kick (mocha: poll calls dispatch_pending_actions.delay()).
        The immediate trigger handles the common path; the short-delay trigger covers the race
        where this task commits just after a concurrent dispatcher read statuses. Both are
        no-ops if nothing is ready (lease + fast-exit), and serialise through
        scheduling_single_worker. """
    sig = Signature(CREATE_AND_LAUNCH_TASK, args=(analysis_id,))
    sig.apply_async()
    sig.apply_async(countdown=3)


def _clear_lease(node_id, version):
    NodeTask.objects.filter(node_version__node_id=node_id, node_version__version=version).update(
        lease_expires=None, leased_by=None, celery_task=None)


def _backoff_node(node_id, version, analysis_id) -> bool:
    """ The worker writes the outcome of its OWN node (not assignment): on a transient error,
        set it back to DIRTY with a future run_after so the single-worker dispatcher re-leases
        it after the delay. Returns False once attempts are exhausted. """
    now = timezone.now()
    node_task = NodeTask.objects.filter(node_version__node_id=node_id, node_version__version=version).first()
    attempts = node_task.attempt_count if node_task else MAX_NODE_ATTEMPTS
    if attempts >= MAX_NODE_ATTEMPTS:
        return False
    delay = next_backoff(attempts)
    NodeTask.objects.filter(node_version__node_id=node_id, node_version__version=version).update(
        run_after=now + timedelta(seconds=delay), lease_expires=None, leased_by=None, celery_task=None)
    with disable_auditlog():
        AnalysisNode.objects.filter(pk=node_id, version=version).update(status=NodeStatus.DIRTY)
    # Prompt re-lease after the backoff window (the periodic sweep is the catch-all safety net).
    Signature(CREATE_AND_LAUNCH_TASK, args=(analysis_id,)).apply_async(countdown=delay)
    return True


def wait_for_task(celery_task):
    # Normally celery would die with "Never call result.get() within a task!"
    # See http://docs.celeryq.org/en/latest/userguide/tasks.html#task-synchronous-subtasks
    with celery.result.allow_join_result():
        #logging.info("*** Wait_for_task: %s", celery_task)
        result = AsyncResult(celery_task)
        result.get(timeout=MINUTE_SECS * 5)


@celery.shared_task(base=AbortableTask)
def update_node_task(node_id, version):
    """ Runs exactly ONE node (all its dependencies are already satisfied - the dispatcher only
        leases ready nodes), writes the node's own outcome, clears its lease, then re-triggers
        the dispatcher so newly-unblocked children get leased. """
    analysis_id = None
    with disable_auditlog():
        try:
            node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
            analysis_id = node.analysis_id  # capture before load()/deletion
        except AnalysisNode.DoesNotExist:  # @UndefinedVariable
            # Node was deleted or version bumped - this task is obsolete, nothing to re-trigger
            return

        try:
            errors = None
            node_errors = node.get_errors()
            if not node_errors:
                try:
                    # Even if no errors now, parent nodes can be removed on us during load causing failure
                    # Also will throw NodeOutOfDateException if node already bumped (before calling expensive load())
                    node.set_node_task_and_status(update_node_task.request.id, NodeStatus.LOADING)
                    node.load()
                    # Check if we need to clear shadow color
                    if node.shadow_color == NodeColors.ERROR and node.is_valid:
                        node.update(shadow_color=None)
                except NodeOutOfDateException:
                    logging.warning("Node %d/%d out of date - exiting", node.pk, node.version)
                    return  # version bumped - reload_analysis_nodes already re-triggered; do NOT re-trigger here
                except OperationalError:
                    # Transient (DB blip / lock timeout). VG node loads are local DB work, not
                    # network calls, so this is rare and usually clears fast - back off and retry
                    # rather than give up, until MAX_NODE_ATTEMPTS.
                    if _backoff_node(node_id, version, analysis_id):
                        return  # node set DIRTY + run_after; delayed re-trigger scheduled
                    status = NodeStatus.CANCELLED  # out of retries
                except NodeConfigurationException:
                    status = NodeStatus.ERROR_CONFIGURATION
                except NodeParentErrorsException:
                    status = NodeStatus.ERROR_WITH_PARENT
                except Exception:
                    errors = get_traceback()
                    status = NodeStatus.ERROR
                else:
                    status = None  # load() set READY itself
            else:
                status = AnalysisNode.get_status_from_errors(node_errors)

            if status is not None:
                try:
                    logging.info("Node %d/%d status: %s errored: %s", node.pk, node.version, status, errors)
                    shadow_color = NodeColors.ERROR if NodeStatus.is_error(status) else None
                    node.update(status=status, errors=errors, shadow_color=shadow_color)
                except (IntegrityError, NodeOutOfDateException) as e:
                    logging.warning("Node %d/%d out of date: {%s} - exiting", node.pk, node.version, e)
                    pass  # out of date or deleted - just ignore
        finally:
            _clear_lease(node_id, version)

    if analysis_id is not None:
        _trigger_rescheduling(analysis_id)


@celery.shared_task(bind=True)
def wait_for_cache_task(self, node_cache_id):
    """ Kept for backwards-compatibility (in-flight messages at deploy time); no longer added to
        scheduling chains - an in-flight cache is now a dependency gate (issue #346). """
    MAX_CHECKS = 60
    try:
        node_cache = NodeCache.objects.get(pk=node_cache_id)
    except NodeCache.DoesNotExist:
        raise CeleryTasksObsoleteException()  # Kills dependent tasks w/o reporting in Rollbar

    status = node_cache.variant_collection.status
    if status in ProcessingStatus.FINISHED_STATES:
        return
    if status not in (ProcessingStatus.CREATED, ProcessingStatus.PROCESSING):
        raise ValueError(f"{node_cache} collection status={status}")
    if self.request.retries >= MAX_CHECKS:
        raise ValueError(f"Timed out after {self.request.retries} checks waiting for {node_cache}")
    logging.debug(f"Waiting on {node_cache}")
    raise self.retry(countdown=1, max_retries=MAX_CHECKS)


@celery.shared_task
def node_cache_task(node_id, version):
    """ Builds a node's NodeCache, then re-triggers the dispatcher so nodes blocked only on this
        cache get leased (issue #346). """
    logging.info("node_cache_task: %s/%d", node_id, version)

    analysis_id = None
    try:
        try:
            node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
            analysis_id = node.analysis_id
        except AnalysisNode.DoesNotExist:  # @UndefinedVariable
            # Node was deleted - this task is obsolete
            return

        node_cache = NodeCache.objects.get(node_version=node.node_version)
        variant_collection = node_cache.variant_collection
        if variant_collection.status != ProcessingStatus.CREATED:
            return

        if not (node.is_valid and node.modifies_parents()):
            logging.debug("Not doing anything for node %s", node.pk)
            variant_collection.status = ProcessingStatus.SKIPPED
            variant_collection.save()
            return

        variant_collection.status = ProcessingStatus.PROCESSING
        variant_collection.save()

        # Node could be READY, we're just generating cache afterwards (node doesn't need but children do)
        if node.status != NodeStatus.READY:
            node.set_node_task_and_status(node_cache_task.request.id, NodeStatus.LOADING_CACHE)

        try:
            node.write_cache(variant_collection)
            processing_status = ProcessingStatus.SUCCESS
            status = NodeStatus.LOADING
        except:
            log_traceback()
            processing_status = ProcessingStatus.ERROR
            status = NodeStatus.ERROR

        variant_collection.status = processing_status
        variant_collection.save()

        if node.status != NodeStatus.READY:
            with disable_auditlog():
                node.status = status
                node.save()
    finally:
        if analysis_id is not None:
            _trigger_rescheduling(analysis_id)


@celery.shared_task
def reschedule_stalled_analyses():
    """ mocha's periodic dispatch loop (Celery beat). Discovery only - finds analyses with work
        that isn't being actively worked: a DIRTY node (waiting to be dispatched) or a node whose
        lease has expired (abandoned by a dead worker), and kicks the single-worker dispatcher.
        This is what makes lease-expiry self-heal after a worker is killed.

        It does NO assignment and makes no timing decisions: backoff (run_after) and the
        reclaim / re-lease / terminal-fail decisions are all enforced authoritatively in
        lease_ready_nodes. This query is deliberately over-inclusive - a kick with nothing ready
        fast-exits in the single worker - so it never misses stalled work. """
    now = timezone.now()
    stalled_lease = Q(nodeversion__nodetask__lease_expires__lt=now)  # abandoned by a dead worker
    dirty = Q(status=NodeStatus.DIRTY)  # waiting to be dispatched (run_after honoured at lease time)
    analysis_ids = (AnalysisNode.objects
                    .filter(status__in=NodeStatus.LOADING_STATUSES)
                    .filter(stalled_lease | dirty)
                    .values_list("analysis_id", flat=True).distinct())

    for analysis_id in analysis_ids:
        # -> scheduling_single_worker (fast-exit if nothing ready)
        Signature(CREATE_AND_LAUNCH_TASK, args=(analysis_id,)).apply_async()


@celery.shared_task
def delete_analysis_old_node_versions(analysis_id):
    """ This is called after creating new versions, so delete anything not latest """
    logging.debug(f"Deleting stale cache for analysis_id={analysis_id}")
    node_versions_qs = NodeVersion.objects.filter(node__analysis_id=analysis_id)
    latest = node_versions_qs.filter(node_id=OuterRef('node_id')).order_by('-version')
    sub_query = Subquery(latest.values('pk')[:1])
    node_versions_qs.annotate(latest_version=sub_query).exclude(pk=F('latest_version')).delete()


@celery.shared_task(bind=True)
def wait_for_node(self, node_id):
    """ Used to build a dependency on a node that's already loading.

        No longer added to scheduling chains (issue #346 removed the worker-starvation pattern);
        kept because analysis_grid_export_tasks._wait_for_output_node() still calls it
        synchronously as a blocking export gate.

        Uses Celery retry (not sleep) to free the worker between checks, preventing deadlocks
        when all workers are occupied waiting for parent nodes.
     """
    TIME_BETWEEN_CHECKS = [5, 5, 10, 10, 30, 30, 60, MINUTE_SECS * 2]
    EVENT_NAME = "wait_for_node"

    try:
        # Any exception will exit this task which is ok
        node = AnalysisNode.objects.get(pk=node_id)

        if NodeStatus.is_ready(node.status):
            logging.info(f"Node {node} status={node.get_status_display()} was ready")
            return

        try:
            node_task = NodeTask.objects.get(node_version__node=node, node_version__version=node.version)
            if node_task.celery_task:
                wait_for_task(node_task.celery_task)
                return
        except NodeTask.DoesNotExist:
            pass

        # loading (eg QUEUED) - there won't be a celery task yet
        if node.status in NodeStatus.LOADING_STATUSES:
            retry_index = self.request.retries
            if retry_index < len(TIME_BETWEEN_CHECKS):
                sleep_time = TIME_BETWEEN_CHECKS[retry_index]
                total_time = sum(TIME_BETWEEN_CHECKS[:retry_index])
                details = f"Waiting on parent node {node_id} which is QUEUED"
                details += f" - waiting for {sleep_time} secs, {total_time} so far!"
                create_event(None, EVENT_NAME, details, severity=LogLevel.WARNING)
                raise self.retry(countdown=sleep_time, max_retries=len(TIME_BETWEEN_CHECKS))
            else:
                total_time = sum(TIME_BETWEEN_CHECKS)
                details = f"Waited on parent node {node_id} for {total_time} seconds. "
                details += "Didn't become available, dying so we don't cause a deadlock!"
                create_event(None, EVENT_NAME, details, severity=LogLevel.ERROR)
        else:
            details = f"Waiting on parent node {node_id} NON LOADING status {node.get_status_display()} no celery task!"
            create_event(None, EVENT_NAME, details, severity=LogLevel.ERROR)
    except celery.exceptions.Retry:
        raise
    except:
        log_traceback()
