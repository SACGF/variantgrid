import json
import logging
import random

from django.conf import settings
from django.db.models import Max
from django.http.response import JsonResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_POST

from analysis.models import AnalysisVariable, AnalysisTemplate, AnalysisTemplateType, NodeCount, \
    AnalysisTemplateVersion, VariantTag
from analysis.models.enums import TagLocation
from analysis.models.nodes import node_utils
from analysis.models.nodes.analysis_node import NodeStatus, AnalysisEdge, NodeVersion, AnalysisNodeAlleleSource, \
    AnalysisNode
from analysis.models.nodes.filter_child import create_filter_child_node
from analysis.models.nodes.filters.built_in_filter_node import BuiltInFilterNode
from analysis.models.nodes.filters.selected_in_parent_node import NodeVariant, SelectedInParentNode
from analysis.models.nodes.filters.venn_node import VennNode
from analysis.models.nodes.node_types import get_node_types_hash_by_class_name
from analysis.models.nodes.node_utils import reload_analysis_nodes, update_analysis, \
    get_toposorted_nodes, get_rendering_dict
from analysis.serializers import VariantTagSerializer
from analysis.views.analysis_permissions import get_analysis_or_404, get_node_subclass_or_404, \
    get_node_subclass_or_non_fatal_exception
from analysis.views.node_json_view import NodeJSONPostView
from library.django_utils import require_superuser
from ontology.models import OntologyTermRelation, OntologyTerm
from ontology.serializers import OntologyTermSerializer
from snpdb.models import Tag, BuiltInFilters, GenomeBuild, Sample
from snpdb.tasks.clingen_tasks import populate_clingen_alleles_from_allele_source


@require_POST
def clone_analysis(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    new_analysis = analysis.clone(request.user)
    reload_analysis_nodes(new_analysis.pk)

    return JsonResponse({"analysis_id": new_analysis.pk})


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
            node.x = params['x']
            node.y = params['y']
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


NODE_TYPES_HASH = None


def node_data(request, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    return JsonResponse(get_rendering_dict(node))


@require_POST
def node_create(request, analysis_id, node_type):
    global NODE_TYPES_HASH
    if NODE_TYPES_HASH is None:
        NODE_TYPES_HASH = get_node_types_hash_by_class_name()

    analysis = get_analysis_or_404(request.user, analysis_id, write=True)

    node_class = NODE_TYPES_HASH[node_type]
    node = node_class.objects.create(analysis=analysis)
    update_analysis(node.analysis_id)
    return JsonResponse(get_rendering_dict(node))


@require_POST
def nodes_copy(request, analysis_id):
    node_ids = json.loads(request.POST["nodes"])
    node_ids = set([int(i) for i in node_ids])

    nodes = []
    edges = []

    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    nodes_qs = analysis.analysisnode_set.filter(id__in=node_ids).select_subclasses()
    topo_sorted = get_toposorted_nodes(nodes_qs)

    old_new_map = {}
    for group in topo_sorted:
        for node in group:
            if analysis_id is None:
                analysis_id = node.analysis_id

            template_node = get_node_subclass_or_404(request.user, node.id)
            parents = template_node.analysisnode_ptr.parents().filter(id__in=old_new_map).values_list('id', flat=True)

            clone_node = template_node.save_clone()
            clone_node.x += 10
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

            if not clone_node.is_valid():
                clone_node.count = None

            clone_node.save()
            nodes.append(get_rendering_dict(clone_node))

    update_analysis(analysis.pk)
    return JsonResponse({"nodes": nodes,
                         "edges": edges})


@require_POST
def nodes_delete(request, analysis_id):
    node_ids = json.loads(request.POST["nodes"])
    node_ids = set([int(i) for i in node_ids])

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

    return JsonResponse(ret)


@require_POST
def set_variant_selected(request, node_id):
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


def create_filter_child(request, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    column_name = request.POST['column_name']
    column_filter = request.POST['column_filter']

    child_node = create_filter_child_node(node, column_name, column_filter)

    data = get_rendering_dict(child_node)
    data["node_id"] = node.get_css_id()
    return JsonResponse(data)


@require_POST
def create_extra_filter_child(request, node_id, extra_filters):
    node = get_node_subclass_or_404(request.user, node_id, write=True)
    x = node.x + 50 + random.randrange(-10, 10)
    y = node.y + 100 + random.randrange(-10, 10)
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


def create_selected_child(request, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    x = node.x + 50 + random.randrange(-10, 10)
    y = node.y + 100 + random.randrange(-10, 10)

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


def get_sample_patient_gene_disease_data(sample: Sample):
    data = {
        "patient_id": sample.patient_id
    }
    if sample.patient:
        all_terms = OntologyTerm.objects.filter(pk__in=sample.patient.get_ontology_term_ids())
        gene_disease_qs = OntologyTermRelation.gene_disease_relations()
        gene_disease_terms = all_terms.filter(subject__in=gene_disease_qs).distinct()
        data["patient"] = str(sample.patient)
        data["total_terms"] = all_terms.count()
        data["terms"] = [OntologyTermSerializer(t).data for t in gene_disease_terms]
    return data


def sample_patient_gene_disease(request, sample_id):
    """ For a sample, return patient MONDO terms that are associated with gene/disease
        Used by MOI Node """
    sample = Sample.get_for_user(request.user, sample_id)
    data = get_sample_patient_gene_disease_data(sample)
    return JsonResponse(data)


@require_POST
def analysis_reload(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    node_utils.reload_analysis_nodes(analysis.pk)
    return JsonResponse({})


def get_node_counts(node_id, version, total_count):
    counts = {}
    if total_count is not None:  # Sub-counts will only be here if total count is
        counts[BuiltInFilters.TOTAL] = total_count
        node_version = NodeVersion.objects.get_or_create(node_id=node_id, version=version)[0]
        for nc in NodeCount.objects.filter(node_version=node_version):
            counts[nc.label] = nc.count
    return counts


def nodes_status(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    nodes = json.loads(request.GET['nodes'])
    qs = analysis.analysisnode_set.filter(id__in=nodes)
    node_status_list = []
    for data in qs.values("id", "version", "status", "count", "shadow_color"):
        node_id = data["id"]
        version = data["version"]

        data["valid"] = not NodeStatus.is_error(data["status"])
        data["ready"] = NodeStatus.is_ready(data["status"])
        data["counts"] = get_node_counts(node_id, version, data["count"])
        node_status_list.append(data)
    return JsonResponse({"node_status": node_status_list})


@require_POST
def analysis_set_panel_size(request, analysis_id):
    """ This is set from AJAX queries, ie dragging a panel border """
    analysis = get_analysis_or_404(request.user, analysis_id, write=True)
    analysis.analysis_panel_fraction = request.POST["analysis_panel_fraction"]
    analysis.save()
    return JsonResponse({})


@require_POST
@require_superuser
def node_populate_clingen_alleles(request, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    an_as, _ = AnalysisNodeAlleleSource.objects.get_or_create(node=node)
    populate_clingen_alleles_from_allele_source.si(an_as.pk, settings.CLINGEN_ALLELE_REGISTRY_MAX_MANUAL_REQUESTS).apply_async()
    return JsonResponse({})


@require_POST
def analysis_template_variable(request, node_id):
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

    # Mark all previous as inactive
    analysis_template.analysistemplateversion_set.all().update(active=False)

    analysis_snapshot = analysis_template.analysis.clone()
    analysis_snapshot.visible = False
    analysis_snapshot.template_type = AnalysisTemplateType.SNAPSHOT
    analysis_snapshot.save()

    sample_gene_list = analysis_snapshot.analysisnode_set.filter(analysisvariable__field='sample_gene_list',
                                                                 analysisvariable__class_name='genes.SampleGeneList')
    requires_sample_gene_list = sample_gene_list.exists()

    data = analysis_template.analysistemplateversion_set.all().aggregate(max_version=Max("version"))
    current_max_version = data.get("max_version") or 0
    version = current_max_version + 1
    AnalysisTemplateVersion.objects.create(template=analysis_template,
                                           version=version,
                                           analysis_name_template=analysis_name_template,
                                           analysis_snapshot=analysis_snapshot,
                                           active=True,
                                           requires_sample_gene_list=requires_sample_gene_list)
    return JsonResponse({"version": version})
