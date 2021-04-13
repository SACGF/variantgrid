from typing import List, Set

from celery.canvas import chain, Signature
from collections import defaultdict
from django.db.models.query_utils import Q
from toposort import toposort
import random
import uuid
import networkx as nx

from analysis.models import Analysis
from analysis.models.nodes.analysis_node import AnalysisNode, AnalysisEdge, NodeStatus
from analysis.tasks.node_update_tasks import add_jobs_to_groups, wait_for_node
from library.celery_utils import execute_task


def get_nodes_by_id(nodes_qs):
    nodes_by_id = {}

    for node in nodes_qs:
        nodes_by_id[node.id] = node
    return nodes_by_id


def get_parent_value_dag_dictionary(nodes):
    """ Create dict that's used for toposort:

        Dependencies are expressed as a dictionary whose keys are items
        and whose values are a set of dependent items aka 'parent value DAG' """

    graph = defaultdict(set)

    for node_id in nodes:
        graph[node_id] = set()

    query = Q(parent__in=nodes) | Q(child__in=nodes)
    qs = AnalysisEdge.objects.filter(query).values_list('parent_id', 'child_id')  # @UndefinedVariable
    for parent_id, child_id in qs:
        graph[child_id].add(parent_id)
    return graph


def get_toposorted_nodes(nodes_qs):
    nodes = get_nodes_by_id(nodes_qs)
    parent_value_data = get_parent_value_dag_dictionary(nodes)
    return get_toposorted_nodes_from_parent_value_data(nodes, parent_value_data)


def get_toposorted_nodes_from_parent_value_data(nodes, parent_value_data):
    topo_sorted_nodes = []
    for grp in toposort(parent_value_data):
        nodes_group = []
        for node_id in grp:
            node = nodes.get(node_id)
            if node:
                nodes_group.append(node)
        if nodes_group:
            topo_sorted_nodes.append(nodes_group)

    return topo_sorted_nodes


def update_nodes(analysis_id, run_async=True):
    tasks = get_analysis_update_tasks(analysis_id)
    for t in tasks:
        execute_task(t, run_async=run_async)


def reload_analysis_nodes(analysis_id, run_async=True):
    for node in AnalysisNode.objects.filter(analysis_id=analysis_id).select_subclasses():
        node.update_children = False  # Will get them all in loop
        node.appearance_dirty = True
        node.queryset_dirty = True
        node.save()
    return update_nodes(analysis_id, run_async=run_async)


def get_ancestor_set(node_id, parent_value_data):
    ancestor_set = set()
    for parent_id in parent_value_data[node_id]:
        ancestor_set.add(parent_id)
        ancestor_set.update(get_ancestor_set(parent_id, parent_value_data))

    return ancestor_set


def get_node_dependencies(nodes_by_id, parent_value_data):
    node_dependencies = set()
    for node_id in nodes_by_id:
        ancestors = get_ancestor_set(node_id, parent_value_data)
        node_dependencies.update(ancestors)

    return node_dependencies


def _add_jobs_for_group(analysis_update_uuid, dependencies, grp, groups, existing_cache_jobs) -> Set:
    cache_jobs = set()
    after_jobs = []
    jobs = []

    for node in grp:
        if node.analysis_update_uuid == analysis_update_uuid:
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

    add_jobs_to_groups(jobs, groups)
    add_jobs_to_groups(after_jobs, groups)
    return cache_jobs


def get_analysis_update_tasks(analysis_id) -> List:
    """ Runs update tasks on nodes that have status=DIRTY """

    tasks = []

    analysis = Analysis.objects.get(pk=analysis_id)
    node_ids = analysis.analysisnode_set.all().values_list("pk", flat=True)
    edges = AnalysisEdge.objects.filter(parent__analysis=analysis).values_list("parent", "child")

    all_nodes_graph = nx.DiGraph()
    all_nodes_graph.add_nodes_from(node_ids)
    all_nodes_graph.add_edges_from(edges)

    for connected_components in nx.weakly_connected_components(all_nodes_graph):
        sub_graph = all_nodes_graph.subgraph(connected_components)
        sub_graph_node_ids = list(sub_graph)
        # nx.topological_sort returns a flattened list, ie doesn’t break into groups which can run in parallel
        # so use other toposort library

        # We need a way to lock/claim the nodes - so someone else calling get_analysis_update_task()
        # doesn't also launch update tasks for them.
        sub_graph_nodes_qs = analysis.analysisnode_set.filter(pk__in=sub_graph_node_ids).select_subclasses()
        dirty_nodes_qs = sub_graph_nodes_qs.filter(status=NodeStatus.DIRTY)
        analysis_update_uuid = uuid.uuid4()
        num_nodes_to_update = dirty_nodes_qs.update(analysis_update_uuid=analysis_update_uuid)

        groups = []
        if num_nodes_to_update:
            # Need to reload nodes to see if they have the analysis_update_uuid on them
            parent_value_data = defaultdict(set)
            for parent, child_list in nx.to_dict_of_lists(sub_graph).items():
                for child_node_id in child_list:
                    parent_value_data[child_node_id].add(parent)

            nodes_by_id = get_nodes_by_id(sub_graph_nodes_qs)
            dependencies = get_node_dependencies(nodes_by_id, parent_value_data)
            topo_sorted = get_toposorted_nodes_from_parent_value_data(nodes_by_id, parent_value_data)

            # Ensure cache loading tasks are only triggered once. Cache can come from different toposort level/groups
            # eg MergeNode asks parent Venn to cache (which is was already doing)
            all_cache_jobs = set()
            for grp in topo_sorted:
                group_cache_jobs = _add_jobs_for_group(analysis_update_uuid, dependencies, grp, groups, all_cache_jobs)
                all_cache_jobs.update(group_cache_jobs)

            update_nodes_qs = analysis.analysisnode_set.filter(analysis_update_uuid=analysis_update_uuid)
            update_nodes_qs.update(analysis_update_uuid=None, status=NodeStatus.QUEUED)

        if groups:
            t = chain(groups)
            tasks.append(t)

    return tasks


def get_rendering_dict(node):
    node_class = node.get_class_name()
    css_classes = node.get_css_classes()
    css_classes.append(node_class)
    css_classes.append("outputEndpoint")

    if node.max_inputs != 0:
        css_classes.append("inputEndpoint")

    node_args = node.get_rendering_args()

    if node.pk:
        node_id = int(node.pk)
    else:
        node_id = None

    style = f"left: {node.x}px; top: {node.y}px"
    attributes = {
        "node_id": node_id,
        "node_class": node.get_node_class_label(),
        "version_id": node.version,
        "appearance_version_id": node.appearance_version,
        "class": " ".join(css_classes),
        "id": node.get_css_id(),
        "style": style,
        "x": node.x,
        "y": node.y
    }
    return {
        "attributes": attributes,
        "node_class": node_class,
        "name": node.name,
        "args": node_args
    }


def create_node(node_class, analysis, **kwargs):
    save = kwargs.pop("save", True)
    node_kwargs = {"x": 10 + random.random() * 50,
                   "y": 10 + random.random() * 20}
    node_kwargs.update(kwargs)
    node_kwargs["analysis"] = analysis
    if save:
        node = node_class.objects.create(**node_kwargs)
    else:
        node = node_class(**node_kwargs)
    return node
