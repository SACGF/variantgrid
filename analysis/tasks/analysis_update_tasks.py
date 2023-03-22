import itertools
import logging
import operator
import uuid
from collections import defaultdict
from functools import reduce
from typing import List, Set

import celery
import networkx as nx
from celery.canvas import Signature, chain
from django.db.models import Q

from analysis.models import Analysis, AnalysisEdge, NodeStatus, NodeTask
from analysis.models.nodes.node_utils import get_nodes_by_id, get_toposorted_nodes_from_parent_value_data
from analysis.tasks.node_update_tasks import wait_for_node
from library.log_utils import log_traceback


@celery.shared_task
def create_and_launch_analysis_tasks(analysis_id, run_async=True):
    """ This is run in a single worker queue so that we avoid race conditions"""
    try:
        tasks = _get_analysis_update_tasks(analysis_id)
    except:
        log_traceback()
        raise

    for t in tasks:
        if run_async:
            t.apply_async()
        else:
            result = t.apply()
            if not result.successful():
                raise Exception(result.result)


def _get_analysis_update_tasks(analysis_id) -> List:
    """ Runs update tasks on nodes that have status=DIRTY """

    tasks = []

    analysis = Analysis.objects.get(pk=analysis_id)
    nodes_by_id = get_nodes_by_id(analysis.analysisnode_set.all().select_subclasses())
    edges = AnalysisEdge.objects.filter(parent__analysis=analysis).values_list("parent", "child")

    all_nodes_graph = nx.DiGraph()
    all_nodes_graph.add_nodes_from(nodes_by_id)
    all_nodes_graph.add_edges_from(edges)

    logging.info("-" * 60)
    for connected_components in nx.weakly_connected_components(all_nodes_graph):
        sub_graph = all_nodes_graph.subgraph(connected_components)
        sub_graph_node_ids = list(sub_graph)
        # nx.topological_sort returns a flattened list, ie doesn’t break into groups which can run in parallel
        # so use other toposort library

        # We need a way to lock/claim the nodes - so someone else calling get_analysis_update_task()
        # doesn't also launch update tasks for them.
        analysis_update_uuid = uuid.uuid4()
        node_task_records = []
        logging.info("Dirty nodes:")
        for node_id in sub_graph_node_ids:
            node = nodes_by_id[node_id]
            if node.status == NodeStatus.DIRTY:
                node_task = NodeTask(node_id=node_id, version=node.version, analysis_update_uuid=analysis_update_uuid)
                logging.info(node_task)
                node_task_records.append(node_task)

        if not node_task_records:
            continue

        NodeTask.objects.bulk_create(node_task_records, ignore_conflicts=True)

        # Return the ones we got the lock for
        node_tasks = NodeTask.objects.filter(analysis_update_uuid=analysis_update_uuid)
        node_versions_to_update = dict(node_tasks.values_list("node_id", "version"))
        logging.info(f"Got lock for: {node_versions_to_update}")

        groups = []
        if node_versions_to_update:
            parent_value_data = defaultdict(set)
            for parent, child_list in nx.to_dict_of_lists(sub_graph).items():
                for child_node_id in child_list:
                    parent_value_data[child_node_id].add(parent)

            dependencies = _get_node_dependencies(nodes_by_id, parent_value_data)
            topo_sorted = get_toposorted_nodes_from_parent_value_data(nodes_by_id, parent_value_data)

            # Ensure cache loading tasks are only triggered once. Cache can come from different toposort level/groups
            # eg MergeNode asks parent Venn to cache (which is was already doing)
            all_cache_jobs = set()
            for grp in topo_sorted:
                group_cache_jobs = _add_jobs_for_group(node_versions_to_update, dependencies, grp, groups, all_cache_jobs)
                all_cache_jobs.update(group_cache_jobs)

            # Need to only set where version matches what we got lock for (as it may have updated)
            node_version_q_list = []
            for node_id, version in node_versions_to_update.items():
                node_version_q_list.append(Q(pk=node_id) & Q(version=version))
            q_node_version = reduce(operator.or_, node_version_q_list)
            analysis.analysisnode_set.filter(q_node_version).update(status=NodeStatus.QUEUED)

        if groups:
            t = _get_celery_workflow_task(groups)
            tasks.append(t)

    return tasks


def _get_celery_workflow_task(groups):
    # So I used to use a Celery chord to be able to execute things at the same level (group) in parallel
    # however, it ended up being 3x slower doing it the 'clever' way due to fighting for DB I/O resources
    # The celery messages could be enormous crashing out RabbitMQ and there would be celery chord jobs
    # left orphaned and executed every few seconds forever. So - just doing it the dumb way.
    task = chain(itertools.chain(*groups))
    return task


def _add_jobs_for_group(nodes_to_update, dependencies, grp, groups, existing_cache_jobs) -> Set:
    cache_jobs = set()
    after_jobs = []
    jobs = []

    for node in grp:
        if node.pk in nodes_to_update:
            update_job = node.get_update_task()

            # Cache jobs are separated, and put in a set to remove dupes such as when >=2 venn's have the same parents
            task_args_objs_set = node.get_cache_task_args_objs_set()
            new_cache_jobs = task_args_objs_set - existing_cache_jobs
            if new_cache_jobs:
                cache_jobs.update(new_cache_jobs)
                after_jobs.append(update_job)
            else:
                jobs.append(update_job)
        elif node.pk in dependencies:
            # Sometimes nodes may be already loading from another update - need to keep dependencies on existing task
            # and wait on loading parent tasks
            if NodeStatus.is_loading(node.status):
                jobs.append(wait_for_node.si(node.pk))  # @UndefinedVariable

    if cache_jobs:
        for task, args, _ in cache_jobs:
            if task:
                jobs.append(Signature(task, args=args, immutable=True))

    groups.append(jobs)
    groups.append(after_jobs)
    return cache_jobs


def _get_ancestor_set(node_id, parent_value_data):
    ancestor_set = set()
    for parent_id in parent_value_data[node_id]:
        ancestor_set.add(parent_id)
        ancestor_set.update(_get_ancestor_set(parent_id, parent_value_data))

    return ancestor_set


def _get_node_dependencies(nodes_by_id, parent_value_data):
    node_dependencies = set()
    for node_id in nodes_by_id:
        ancestors = _get_ancestor_set(node_id, parent_value_data)
        node_dependencies.update(ancestors)

    return node_dependencies
