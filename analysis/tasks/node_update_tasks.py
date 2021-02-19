import celery
from celery.canvas import chord
from celery.contrib.abortable import AbortableTask
from celery.result import AsyncResult, allow_join_result
from django.db.utils import OperationalError
import logging
from time import sleep

from analysis.exceptions import NodeConfigurationException, NodeParentErrorsException
from analysis.models.nodes.analysis_node import AnalysisNode, NodeStatus, NodeVersion, NodeCache
from eventlog.models import create_event
from library.constants import MINUTE_SECS
from library.enums.log_level import LogLevel
from library.log_utils import log_traceback, get_traceback
from snpdb.models import ProcessingStatus


@celery.task
def dummy_task():
    pass


def add_jobs_to_groups(jobs, groups):
    if jobs:
        if len(jobs) > 1:
            # Use a chord, as groups don't synchronize inside a chain
            g = chord(jobs, dummy_task.si())  # @UndefinedVariable
        else:
            g = jobs[0]  # @UndefinedVariable

        groups.append(g)


def wait_for_task(celery_task):
    # Normally celery would die with "Never call result.get() within a task!"
    # See http://docs.celeryq.org/en/latest/userguide/tasks.html#task-synchronous-subtasks
    with allow_join_result():
        #logging.info("*** Wait_for_task: %s", celery_task)
        result = AsyncResult(celery_task)
        result.get(timeout=MINUTE_SECS * 5)


@celery.task(base=AbortableTask)
def update_node_task(node_id, version):
    try:
        node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
    except AnalysisNode.DoesNotExist:  # @UndefinedVariable
        # Node was deleted - this task is obsolete
        return

    errors = None
    node_errors = node.get_errors()
    if not node_errors:
        try:
            # Even if no errors now, parent nodes can be removed on us during load causing failure
            node.set_node_task_and_status(update_node_task.request.id, NodeStatus.LOADING)
            node.load()
            return  # load saves everything, no need to save again below
        except OperationalError:
            status = NodeStatus.CANCELLED
        except NodeConfigurationException:
            status = NodeStatus.ERROR_CONFIGURATION
        except NodeParentErrorsException:
            status = NodeStatus.ERROR_WITH_PARENT
        except Exception:
            errors = get_traceback()
            status = NodeStatus.ERROR
    else:
        status = AnalysisNode.get_status_from_errors(node_errors)

    if status:
        node.status = status
        node.celery_task = None
        node.db_pid = None
        node.errors = errors
        node.save()  # Will set background color etc based on errors


@celery.task
def wait_for_cache_task(node_cache_id):
    MAX_CHECKS = 60
    num_checks = 0
    while True:
        node_cache = NodeCache.objects.get(pk=node_cache_id)
        status = node_cache.variant_collection.status
        if status in ProcessingStatus.FINISHED_STATES:
            print(f"{node_cache} DONE")
            break
        elif status not in (ProcessingStatus.CREATED, ProcessingStatus.PROCESSING):
            raise ValueError(f"{node_cache} collection status={status}")
        num_checks += 1
        if num_checks > MAX_CHECKS:
            raise ValueError(f"Timed out after {num_checks} waiting for {node_cache}")
        sleep(1)
        logging.debug(f"Waiting on {node_cache}")


@celery.task
def node_cache_task(node_id, version):
    logging.info("node_cache_task: %s/%d", node_id, version)

    try:
        node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
    except AnalysisNode.DoesNotExist:  # @UndefinedVariable
        # Node was deleted - this task is obsolete
        return

    node_cache = NodeCache.objects.get(node_version=node.node_version)
    variant_collection = node_cache.variant_collection
    if variant_collection.status != ProcessingStatus.CREATED:
        return

    if not (node.is_valid() and node.modifies_parents()):
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
        node.status = status
        node.celery_task = None
        node.db_pid = None
        node.save()


@celery.task
def delete_old_node_versions(node_id, version):
    """ This is called before bumping the version, so delete that version and earlier """
    logging.debug(f"Deleting stale cache for node_id={node_id}, version_id={version}")
    NodeVersion.objects.filter(node_id=node_id, version__lte=version).delete()


@celery.task
def wait_for_node(node_id):
    """ Used to build a dependency on a node that's already loading.

        The danger here is that we'll end up taking up a celery worker, waiting for
        the parent to become available... so we have to be careful
     """

    EVENT_NAME = "wait_for_node"

    #logging.info("wait_for_node: %s", node_id)
    try:
        TIME_BETWEEN_CHECKS = [1, 1, 2, 5, 60, MINUTE_SECS * 5]
        total_time = 0
        for sleep_time in TIME_BETWEEN_CHECKS:

            # Any exception will exit this task which is ok
            node = AnalysisNode.objects.get(pk=node_id)
            #logging.info("node status: %s", node.status)

            if NodeStatus.is_ready(node.status):
                #logging.info("Node was ready")
                return

            if node.celery_task:
                wait_for_task(node.celery_task)
                return
            else:
                # QUEUED - there won't be a celery task yet
                if node.status == NodeStatus.QUEUED:
                    details = f"Waiting on parent node {node_id} which is QUEUED"
                    details += f" - waiting for {sleep_time} secs, {total_time} so far!"
                    create_event(None, EVENT_NAME, details, severity=LogLevel.WARNING)
                else:
                    details = f"Waiting on parent node {node_id} status {node.status} no celery task!"
                    create_event(None, EVENT_NAME, details, severity=LogLevel.ERROR)
                    return

            logging.info("Sleeping for %d seconds", sleep_time)
            sleep(sleep_time)
            total_time += sleep_time

        # Timeout
        details = f"Waited on parent node {node_id} for {total_time} seconds. "
        details += "Didn't become available, dying so we don't cause a deadlock!"
        create_event(None, EVENT_NAME, details, severity=LogLevel.ERROR)
    except:
        log_traceback()
