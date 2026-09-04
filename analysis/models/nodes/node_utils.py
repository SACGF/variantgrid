import json
import logging
import random
import time
from collections import defaultdict
from dataclasses import asdict

from auditlog.context import disable_auditlog
from celery.canvas import Signature
from django.db import connection
from django.db.models import F
from django.db.models.query_utils import Q
from django.utils import timezone
from toposort import toposort

from analysis.exceptions import NonFatalNodeError
from analysis.models import Analysis, NodeColors, NodeStatus
from analysis.models.nodes.analysis_node import AnalysisEdge, NodeVersion
from analysis.models.nodes.node_counts import get_tag_node_counts_dict, get_tagged_variant_ids_by_label
from analysis.tasks.node_update_tasks import delete_analysis_old_node_versions
from library.utils import add_exception_note
from snpdb.models.models_enums import TagFilter


def get_nodes_by_id(nodes_qs):
    nodes_by_id = {}

    for node in nodes_qs:
        nodes_by_id[node.id] = node
    return nodes_by_id


def get_child_position(parent):
    """ Children go along the flow from their parent - to the right in horizontal mode, below in vertical """
    if parent.analysis.analysis_horizontal_mode:
        child_x_offset, child_y_offset = 100, 50
    else:
        child_x_offset, child_y_offset = 50, 100
    x = parent.x + child_x_offset + random.randrange(-10, 10)
    y = parent.y + child_y_offset + random.randrange(-10, 10)
    return x, y


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


def update_analysis(analysis_id):
    """ Launches async job to update analysis """

    delete_analysis_old_node_versions.si(analysis_id).apply_async()

    task = Signature("analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks", args=(analysis_id,))
    task.apply_async()


def update_analysis_tag_node_counts(analysis: Analysis, tag_labels=None):
    """ Adding/removing a tag doesn't bump node versions, so the usual reload doesn't recount.
        Recount the tag node counts in place against the versions the nodes are already on. Building
        each node's queryset costs ~15 queries, so this runs in a task off the back of tagging
        @see analysis.tasks.variant_tag_tasks.
        tag_labels - restrict to these (default: every tag node count the analysis has configured) """
    configured_tag_labels = {label for label, _ in analysis.get_node_count_types() if TagFilter.get_tag_id(label)}
    if tag_labels is not None:
        configured_tag_labels &= set(tag_labels)
    if not configured_tag_labels:
        return

    # The tagged variants are the same for every node, so look them up once for the whole analysis
    tagged_variant_ids_by_label = get_tagged_variant_ids_by_label(analysis, configured_tag_labels)

    counts_by_node_version_id = {}
    for node in analysis.analysisnode_set.filter(status=NodeStatus.READY).select_subclasses():
        node_version = NodeVersion.objects.filter(node=node, version=node.version).first()
        if node_version is None:
            continue  # Node reloaded from under us - it'll count these itself
        try:
            counts_by_node_version_id[node_version.pk] = get_tag_node_counts_dict(node, tagged_variant_ids_by_label)
        except NonFatalNodeError:
            # Node is ready but an ancestor isn't (eg the analysis is mid-reload) so we can't build its
            # query - it counts these itself when it loads
            continue

    # Merge into the labels the load wrote rather than replacing them - this runs concurrently with
    # loads, and only computes the tag labels. A node that bumped its version while we were counting
    # has had this row deleted, so its UPDATE matches nothing and it counts these itself when it reloads.
    # "modified" is the client's signal that this recount landed @see nodes_status
    sql = """UPDATE analysis_nodeversion
             SET load_data = jsonb_set(load_data, '{counts}',
                                       COALESCE(load_data->'counts', '{}'::jsonb) || %s::jsonb),
                 modified = %s
             WHERE id = %s"""
    now = timezone.now()
    with connection.cursor() as cursor:
        for node_version_id, tag_node_counts in counts_by_node_version_id.items():
            cursor.execute(sql, [json.dumps(tag_node_counts), now, node_version_id])


def reload_analysis_nodes(analysis_id, only_errors=False):
    """ only_errors: only reload nodes with error statuses (children of error nodes
        have ERROR_WITH_PARENT status, so the error set is closed under descendants) """
    with disable_auditlog():
        start = time.time()
        analysis = Analysis.objects.get(pk=analysis_id)
        nodes_qs = analysis.analysisnode_set.all()
        nodes_by_id = get_nodes_by_id(nodes_qs.select_subclasses())
        parents = defaultdict(list)
        for parent, child in AnalysisEdge.objects.filter(parent__analysis=analysis).values_list("parent", "child"):
            parents[child].append(nodes_by_id[parent])

        analysis_errors = analysis.get_errors()
        valid_nodes = []
        invalid_nodes = []
        for node_id, node in nodes_by_id.items():
            node._cached_analysis_errors = analysis_errors
            node._cached_parents = parents.get(node_id, [])
            if only_errors and not NodeStatus.is_error(node.status):
                continue
            if node.get_errors():
                invalid_nodes.append(node_id)
            else:
                valid_nodes.append(node_id)
        num_nodes = len(valid_nodes) + len(invalid_nodes)

        update_kwargs = {
            "status": NodeStatus.DIRTY,
            "count": None,
            "errors": None,
            "cloned_from": None,
            "version": F("version") + 1,
            "appearance_version": F("appearance_version") + 1,
        }
        if valid_nodes:
            nodes_qs.filter(pk__in=valid_nodes).update(valid=True, shadow_color=NodeColors.VALID, **update_kwargs)
        if invalid_nodes:
            nodes_qs.filter(pk__in=invalid_nodes).update(valid=False, shadow_color=NodeColors.ERROR, **update_kwargs)

        node_versions = []
        reloaded_qs = nodes_qs.filter(pk__in=valid_nodes + invalid_nodes)
        for node_id, version in reloaded_qs.values_list("pk", "version"):
            node_versions.append(NodeVersion(node_id=node_id, version=version))
        if node_versions:
            NodeVersion.objects.bulk_create(node_versions, ignore_conflicts=True)

        Analysis.objects.filter(pk=analysis_id).update(modified=timezone.now())
        end = time.time()
        logging.info("%d saves took %.2f secs", num_nodes, end-start)
        return update_analysis(analysis_id)


def get_rendering_dict(node):
    node_class = node.get_class_name()
    # Need to add 'node-overlay' so we can find it again
    css_classes = ["node-overlay", node_class] + node.get_css_classes()
    node_args = node.get_rendering_args()

    if node.pk:
        node_id = int(node.pk)
    else:
        node_id = None

    try:
        node_class_label = node.get_node_class_label()
        # This has happened a few times, help debug problem
    except NotImplementedError as e:
        add_exception_note(e, f"Node {node_id=}: {e}")
        raise e

    style = f"left: {node.x}px; top: {node.y}px"
    attributes = {
        "node_id": node_id,
        "node_class": node_class_label,
        "node_classification": node.get_node_classification(),  # Card colour - see analysis_nodes.css
        "version_id": node.version,
        "appearance_version_id": node.appearance_version,
        "id": node.get_css_id(),
        "style": style,
        "input_endpoint": node.max_inputs != 0,
        "output_endpoint": True,  # Can always have output
        "x": node.x,
        "y": node.y
    }
    return {
        "attributes": attributes,
        "node_class": node_class,
        "overlay_css_classes": " ".join(css_classes),
        "name": node.name,
        "icon": asdict(node.get_node_icon()),
        "class_label_short": node.get_node_strip_label(),
        "chips": [asdict(chip) for chip in node.get_node_chips()],
        "args": node_args
    }
