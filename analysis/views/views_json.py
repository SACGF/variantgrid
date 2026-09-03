import json
import logging
import random
from collections import Counter, defaultdict

from celery.result import AsyncResult
from django.conf import settings
from django.db.models import F
from django.http.response import JsonResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.cache import never_cache
from django.views.decorators.http import require_POST

from analysis.models import (
    AnalysisTemplate,
    AnalysisVariable,
    Candidate,
    CandidateSearchRun,
    CandidateStatus,
    NodeCount,
    TagNode,
    VariantTag,
)
from analysis.models.enums import TagLocation, TagNodeMode
from analysis.models.nodes import node_utils
from analysis.models.nodes.analysis_node import AnalysisEdge, AnalysisNode, NodeStatus, NodeTask, NodeVersion
from analysis.models.nodes.filter_child import create_filter_child_node
from analysis.models.nodes.filters.built_in_filter_node import BuiltInFilterNode
from analysis.models.nodes.filters.selected_in_parent_node import NodeVariant, SelectedInParentNode
from analysis.models.nodes.filters.venn_node import VennNode
from analysis.models.nodes.node_types import get_menu_entries_by_key
from analysis.models.nodes.node_utils import (
    get_child_position,
    get_rendering_dict,
    get_toposorted_nodes,
    reload_analysis_nodes,
    update_analysis,
)
from analysis.serializers import CandidateSearchRunSerializer, VariantTagSerializer
from analysis.tasks.analysis_update_tasks import populate_clingen_alleles_from_analysis_node
from analysis.views.analysis_permissions import (
    get_analysis_or_404,
    get_node_subclass_or_404,
    get_node_subclass_or_non_fatal_exception,
)
from analysis.views.node_json_view import NodeJSONPostView
from library.django_utils import require_superuser
from ontology.models import OntologyTerm, OntologyVersion
from ontology.serializers import OntologyTermSerializer
from snpdb.models import BuiltInFilters, GenomeBuild, Sample, Tag
from snpdb.models.models_enums import TagFilter
from variantgrid.celery import app


@require_POST
def clone_analysis(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    new_analysis = analysis.clone(request.user)
    reload_analysis_nodes(new_analysis.pk)

    return JsonResponse({"analysis_id": new_analysis.pk})


@require_POST
def node_reveal_hidden(request, analysis_id, node_id):
    """ Running a template hides a node that errors while being configured, and everything below it -
        put that branch back so the user can see it and why it failed """
    node = get_node_subclass_or_404(request.user, node_id, write=True)
    branch_node_ids = {node.pk} | {n.pk for n in node.descendants_set()}
    hidden_qs = AnalysisNode.objects.filter(pk__in=branch_node_ids, visible=False)
    revealed_ids = list(hidden_qs.values_list("pk", flat=True))
    if revealed_ids:
        hidden_qs.update(visible=True, appearance_version=F("appearance_version") + 1)
        reload_analysis_nodes(node.analysis_id)

    # Same shape as nodes_copy, so the page can drop the branch onto the canvas without a reload
    nodes = []
    edges = []
    for revealed in AnalysisNode.objects.filter(pk__in=revealed_ids).select_subclasses():
        nodes.append(get_rendering_dict(revealed))
        for parent in revealed.analysisnode_ptr.parents():
            edges.append(revealed.get_connection_data(parent))
    return JsonResponse({"nodes": nodes, "edges": edges})


@never_cache
def analysis_node_versions(request, analysis_id):
    """ Returns a dict of {'node_versions' : [node.pk, node.version, node.appearance_version]} """
    analysis = get_analysis_or_404(request.user, analysis_id)
    nodes_qs = analysis.analysisnode_set.filter(visible=True)
    node_versions = nodes_qs.values_list("pk", "version", "appearance_version")
    return JsonResponse({"node_versions": list(node_versions)})


class NodeUpdate(NodeJSONPostView):
    def _get_node(self, request, **kwargs) -> AnalysisNode:
        node_id = kwargs["node_id"]
        return get_node_subclass_or_non_fatal_exception(request.user, node_id, write=True)

    def _get_data(self, request, node, *args, **kwargs):
        op = request.POST["op"]
        params = json.loads(request.POST['params'])
        # Load subclass, need to use OO to handle adding/deleting and save
        logging.debug("Loaded node %s version %d", node.pk, node.version)

        if op == "move":
            node.x = int(params['x'])
            node.y = int(params['y'])
            node.save()
        elif op == "update_connection":
            parent_id = params["parent_id"]
            parent = get_node_subclass_or_non_fatal_exception(request.user, parent_id, write=True)
            if params.get("remove"):
                node.remove_parent(parent)
            else:
                kwargs = {}
                if isinstance(node, VennNode):
                    side = params.get("side")  # Optional for venn
                    kwargs["side"] = side

                node.add_parent(parent, **kwargs)
            node.save()
            update_analysis(node.analysis_id)
        else:
            raise ValueError(f"Unknown operation '{op}'")

        return {}


MENU_ENTRIES_BY_KEY = None


@never_cache
def node_data(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    return JsonResponse(get_rendering_dict(node))


@require_POST
def node_create(request, analysis_id, node_type):
    global MENU_ENTRIES_BY_KEY
    if MENU_ENTRIES_BY_KEY is None:
        MENU_ENTRIES_BY_KEY = get_menu_entries_by_key()

    analysis = get_analysis_or_404(request.user, analysis_id, write=True)

    # A menu entry is a class plus the field values it stamps on the new node (eg a source level)
    node_class, menu_entry = MENU_ENTRIES_BY_KEY[node_type]
    # New nodes go at the start of the flow - the top in vertical mode, the left edge in horizontal
    if analysis.analysis_horizontal_mode:
        x = 50 + random.random() * 20
        y = 10 + random.random() * 50
    else:
        x = 10 + random.random() * 50
        y = 50 + random.random() * 20
    node = node_class.objects.create(analysis=analysis, x=x, y=y, **menu_entry.initial_kwargs)
    update_analysis(node.analysis_id)
    return JsonResponse(get_rendering_dict(node))


@require_POST
def nodes_copy(request, analysis_id):
    node_ids = json.loads(request.POST["nodes"])
    node_ids = {int(i) for i in node_ids}

    nodes = []
    edges = []

    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    nodes_qs = analysis.analysisnode_set.filter(id__in=node_ids).select_subclasses()
    topo_sorted = get_toposorted_nodes(nodes_qs)

    # Nudge the copy clear of the original - along the flow in horizontal mode, so it reads as the next node
    copy_x_offset = 80 if analysis.analysis_horizontal_mode else 10

    old_new_map = {}
    for group in topo_sorted:
        for node in group:
            if analysis_id is None:
                analysis_id = node.analysis_id

            template_node = get_node_subclass_or_404(request.user, node.id)
            parents = template_node.analysisnode_ptr.parents().filter(id__in=old_new_map).values_list('id', flat=True)

            clone_node = template_node.save_clone()
            clone_node.x += copy_x_offset
            clone_node.y += 10
            clone_node.status = NodeStatus.DIRTY
            clone_node.save()
            old_new_map[node.id] = clone_node

            clone_node.adjust_cloned_parents(old_new_map)

            for parent_id in parents:
                new_parent = old_new_map[parent_id]
                new_parent.add_child(clone_node)

                edge = clone_node.get_connection_data(new_parent)
                edges.append(edge)

            if not clone_node.is_valid:
                clone_node.count = None

            clone_node.save()
            nodes.append(get_rendering_dict(clone_node))

    update_analysis(analysis.pk)
    return JsonResponse({"nodes": nodes,
                         "edges": edges})


@require_POST
def nodes_delete(request, analysis_id):
    node_ids = json.loads(request.POST["nodes"])
    node_ids = {int(i) for i in node_ids}

    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    nodes_qs = analysis.analysisnode_set.filter(id__in=node_ids).select_subclasses()
    topo_sorted = get_toposorted_nodes(nodes_qs)

    for group in reversed(topo_sorted):
        for node in group:
            if node.id not in node_ids:
                continue  # parent we don't care about
            # Detach first
            for kid in node.analysisnode_ptr.children.select_subclasses():
                kid.remove_parent(node)
                kid.parents_changed = True
                kid.save()

            node.delete()

    update_analysis(analysis.pk)
    return JsonResponse({})


@require_POST
def set_variant_tag(request, location):
    """ Can be called from analysis or variant details page (location = A or V) """
    location = TagLocation(location)
    variant_id = request.POST['variant_id']
    tag_id = request.POST['tag_id']
    op = request.POST['op']

    # Optional
    variant_tag_id = request.POST.get('variant_tag_id')  # Pass in PK to delete
    genome_build_name = request.POST.get('genome_build_name')

    if analysis_id := request.POST.get('analysis_id'):
        node_id = request.POST.get('node_id')
        analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    else:
        analysis = None
        node_id = None

    tag = get_object_or_404(Tag, pk=tag_id)
    ret = {}  # Empty
    if op == 'add':
        if analysis:
            genome_build = analysis.genome_build
            variant_tag, created = VariantTag.objects.get_or_create(variant_id=variant_id, tag=tag,
                                                                    genome_build=genome_build, location=location,
                                                                    analysis=analysis, user=request.user)
            if node_id:
                variant_tag.node_id = node_id
                # Stamp what the node was showing, so a reviewer can later tell why the variant was in it
                node_version = NodeVersion.objects.filter(node_id=node_id,
                                                          version=F("node__version")).first()
                variant_tag.node_version = node_version
                variant_tag.node_live_data_sources = node_version.live_data_sources if node_version else {}
                variant_tag.save()
        else:
            if genome_build_name is None:
                raise ValueError("Adding requires either 'analysis_id' or 'genome_build_name'")
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

            variant_tag, created = VariantTag.objects.get_or_create(variant_id=variant_id, tag=tag,
                                                                    analysis=None, location=location,
                                                                    user=request.user,
                                                                    defaults={"genome_build": genome_build})
        if created:  # Only return new if anything created
            ret = VariantTagSerializer(variant_tag, context={"request": request}).data
    elif op == 'del':
        # Deletion of tags is for analysis (all users)
        if analysis:
            VariantTag.objects.filter(variant_id=variant_id, analysis=analysis, tag=tag).delete()
        elif variant_tag_id:
            variant_tag = VariantTag.get_for_user(request.user, pk=variant_tag_id, write=True)
            variant_tag.delete()
        else:
            raise ValueError("Deletion requires either 'analysis_id' or 'variant_tag_id'")

    if analysis:
        # Tagging can add/remove that tag's own node count, so tell the client what to draw
        ret["node_count_types"] = analysis.get_node_count_types()

    return JsonResponse(ret)


@require_POST
def set_variant_selected(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id, write=True)
    variant_id = request.POST['variant_id']
    checked = json.loads(request.POST['checked'])

    kwargs = {"variant_id": variant_id,
              "node_id": node.pk}
    if checked:
        NodeVariant.objects.get_or_create(**kwargs)
    else:
        NodeVariant.objects.filter(**kwargs).delete()

    kids_qs = AnalysisEdge.objects.filter(parent=node).values_list("child_id", flat=True)  # @UndefinedVariable
    for node in SelectedInParentNode.objects.filter(pk__in=kids_qs):
        node.queryset_dirty = True
        node.save()

    update_analysis(node.analysis.pk)
    return JsonResponse({})


def create_filter_child(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    column_name = request.POST['column_name']
    column_filter = request.POST['column_filter']

    child_node = create_filter_child_node(node, column_name, column_filter)

    data = get_rendering_dict(child_node)
    data["node_id"] = node.get_css_id()
    return JsonResponse(data)


@require_POST
def create_extra_filter_child(request, analysis_id, node_id, extra_filters):
    node = get_node_subclass_or_404(request.user, node_id, write=True)
    x, y = get_child_position(node)
    if tag_ids := TagFilter.get_tag_ids(extra_filters):
        tag_kwargs = {"mode": TagNodeMode.THIS_ANALYSIS}
        if isinstance(node, TagNode):
            # The parent's scope and cutoff are what produced the rows being filtered - keep them
            tag_kwargs = {"mode": node.mode, "tagged_within_days": node.tagged_within_days}
        tag_node = TagNode.objects.create(analysis=node.analysis,
                                          x=x,
                                          y=y,
                                          ready=False,
                                          **tag_kwargs)
        for tag_id in tag_ids:
            tag_node.tagnodetag_set.create(tag_id=tag_id)
        # Re-load so the node name picks up the tags - TagNode.tag_ids is cached from the create() above
        filter_node = TagNode.objects.get(pk=tag_node.pk)
    else:
        filter_node = BuiltInFilterNode.objects.create(analysis=node.analysis,
                                                       built_in_filter=extra_filters,
                                                       x=x,
                                                       y=y,
                                                       ready=False)
    filter_node.add_parent(node)
    filter_node.save()

    update_analysis(node.analysis.pk)
    data = get_rendering_dict(filter_node)
    data["node_id"] = node.get_css_id()
    return JsonResponse(data)


def create_selected_child(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    x, y = get_child_position(node)

    selected_node = SelectedInParentNode.objects.create(analysis=node.analysis,
                                                        x=x,
                                                        y=y,
                                                        ready=False)
    selected_node.add_parent(node)
    selected_node.save()
    update_analysis(node.analysis.pk)

    data = get_rendering_dict(selected_node)
    data["node_id"] = node.get_css_id()
    return JsonResponse(data)


def get_sample_patient_gene_disease_data(sample: Sample, ontology_version: OntologyVersion):
    data = {
        "patient_id": sample.patient_id
    }
    if sample.patient:
        all_terms = OntologyTerm.objects.filter(pk__in=sample.patient.get_ontology_term_ids())
        gene_disease_qs = ontology_version.get_gene_disease_relations_qs()
        gene_disease_terms = all_terms.filter(subject__in=gene_disease_qs).distinct()
        data["patient"] = str(sample.patient)
        data["total_terms"] = all_terms.count()
        data["terms"] = [OntologyTermSerializer(t).data for t in gene_disease_terms]
    return data


def sample_patient_gene_disease(request, sample_id, ontology_version_id):
    """ For a sample, return patient MONDO terms that are associated with gene/disease
        Used by MOI Node """
    sample = Sample.get_for_user(request.user, sample_id)
    ontology_version = get_object_or_404(OntologyVersion, pk=ontology_version_id)
    data = get_sample_patient_gene_disease_data(sample, ontology_version)
    return JsonResponse(data)


@require_POST
def analysis_reload(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    only_errors = json.loads(request.POST.get("only_errors", "false"))
    node_utils.reload_analysis_nodes(analysis.pk, only_errors=only_errors)
    return JsonResponse({})


@never_cache
def nodes_status(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    nodes = json.loads(request.GET['nodes'])
    node_counts_qs = NodeCount.objects.filter(node_version__node__in=nodes)
    node_counts = defaultdict(dict)
    # Tag counts are recalculated without bumping the node version (@see update_analysis_tag_node_counts),
    # so hand the client the last time the counts changed - that's how it knows a recount landed
    counts_modified = {}
    for node_id, version, label, count, modified in node_counts_qs.values_list(
            "node_version__node_id", "node_version__version", "label", "count", "modified"):
        key = f"{node_id}_{version}"
        node_counts[key][label] = count
        if modified and modified.isoformat() > counts_modified.get(key, ""):
            counts_modified[key] = modified.isoformat()

    node_version_qs = NodeVersion.objects.filter(node__in=nodes)
    live_data_sources = {f"{node_id}_{version}": sources
                         for node_id, version, sources in node_version_qs.values_list("node_id", "version",
                                                                                      "live_data_sources")}

    qs = analysis.analysisnode_set.filter(id__in=nodes)
    node_status_list = []
    for data in qs.values("id", "version", "status", "count", "shadow_color"):
        node_id = data["id"]
        version = data["version"]

        data["valid"] = not NodeStatus.is_error(data["status"])
        data["ready"] = NodeStatus.is_ready(data["status"])

        counts = node_counts.get(f"{node_id}_{version}", {})
        counts[BuiltInFilters.TOTAL] = data["count"]
        data["counts"] = counts
        data["counts_modified"] = counts_modified.get(f"{node_id}_{version}", "")
        # A node reading mutable tables has an advisory count - the client shows it abbreviated (#235)
        sources = live_data_sources.get(f"{node_id}_{version}") or {}
        data["deterministic"] = not sources
        data["live_data_sources"] = sources
        node_status_list.append(data)
    return JsonResponse({"node_status": node_status_list})


@never_cache
def nodes_tasks(request, analysis_id):
    """ Returns a dict of node tasks states / counts """
    NODE_NAME = 'celery@analysis_workers'
    analysis = get_analysis_or_404(request.user, analysis_id)
    inspect = app.control.inspect([NODE_NAME])
    if active := inspect.active():
        active_jobs = {a["id"] for a in active[NODE_NAME]}

        summary = Counter()
        # Current-version NodeTasks for this analysis's visible, still-loading nodes
        node_task_qs = NodeTask.objects.filter(
            node_version__node__analysis=analysis,
            node_version__node__visible=True,
            node_version__node__status__in=NodeStatus.LOADING_STATUSES,
            node_version__version=F("node_version__node__version"))
        for celery_task in node_task_qs.values_list("celery_task", flat=True):
            status = "QUEUED"
            if celery_task:
                if celery_task in active_jobs:
                    status = "ACTIVE"
                else:
                    result = AsyncResult(celery_task)
                    status = result.status
            summary[status] += 1
        data = dict(summary)
    else:
        data = {"error": "No analysis workers found!"}
    return JsonResponse(data)


@require_POST
def analysis_set_panel_size(request, analysis_id):
    """ This is set from AJAX queries, ie dragging a panel border """
    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    analysis.analysis_panel_fraction = request.POST["analysis_panel_fraction"]
    analysis.save()
    return JsonResponse({})


@require_POST
@require_superuser
def node_populate_clingen_alleles(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    node.analysis.check_can_view(request.user)
    populate_clingen_alleles_from_analysis_node.si(node_id, settings.CLINGEN_ALLELE_REGISTRY_MAX_MANUAL_REQUESTS).apply_async()
    return JsonResponse({})


@require_POST
def analysis_template_variable(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id, write=True)

    field = request.POST["field"]
    operation = request.POST["op"]

    kwargs = {"node": node, "field": field}
    if operation == 'add':
        class_name = AnalysisVariable.get_node_field_class_name(node, field)
        AnalysisVariable.objects.get_or_create(**kwargs, defaults={"class_name": class_name})
    elif operation == 'del':
        AnalysisVariable.objects.filter(**kwargs).delete()

    return JsonResponse({})


@require_POST
def analysis_template_save(request, pk):
    """ Creates a new AnalysisTemplateVersion for an AnalysisTemplate """
    analysis_template = AnalysisTemplate.get_for_user(request.user, pk, write=True)
    analysis_name_template = request.POST.get("analysis_name_template")
    activate = request.POST.get("activate") not in (None, "", "0", "false")

    try:
        previously_active = analysis_template.active
        atv = analysis_template.new_version(analysis_name_template)
        replaced_version = None
        if activate:
            atv.activate()
            if previously_active:
                replaced_version = previously_active.version
        return JsonResponse({"version": atv.version, "created": atv.created.isoformat(),
                             "active": atv.active, "replaced_version": replaced_version})
    except ValueError:
        return JsonResponse({
            "error": f"Could not create new analysis template version for '{analysis_template}'"
        })


@require_POST
def analysis_template_clone(request, pk):
    analysis_template = AnalysisTemplate.get_for_user(request.user, pk, write=False)
    new_analysis_template = analysis_template.clone(request.user)
    reload_analysis_nodes(new_analysis_template.analysis.pk)
    return JsonResponse({"analysis_template_id": new_analysis_template.pk})


def get_candidate_search_run_json(request, pk):
    candidate_search_run = CandidateSearchRun.get_for_user(request.user, pk=pk)
    serializer = CandidateSearchRunSerializer(candidate_search_run)
    return JsonResponse(serializer.data)


@require_POST
def set_candidate_status(request, candidate_id):
    candidate = Candidate.get_permission_check(candidate_id, request.user, write=True)
    status = request.POST.get("status")
    candidate.status = CandidateStatus(status)
    candidate.reviewer = request.user
    candidate.save()
    return JsonResponse({"success": True})
