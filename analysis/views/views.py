import json
from collections import Counter, defaultdict

from celery.contrib.abortable import AbortableAsyncResult
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied, EmptyResultSet
from django.db.models import Count
from django.forms.models import inlineformset_factory
from django.http.response import HttpResponse, HttpResponseRedirect, JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.test.client import RequestFactory
from django.urls.base import reverse
from django.utils import timezone
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST
from django.views.decorators.vary import vary_on_cookie
import importlib
import inspect
import logging

from htmlmin.decorators import not_minified_response
from lazy import lazy

from analysis import forms
from analysis.analysis_templates import populate_analysis_from_template_run
from analysis.exceptions import NonFatalNodeError, NodeOutOfDateException
from analysis.forms import SelectGridColumnForm, UserTrioWizardForm, VCFLocusFilterForm, InputSamplesForm, \
    AnalysisChoiceForm, AnalysisTemplateTypeChoiceForm
from analysis.graphs.column_boxplot_graph import ColumnBoxplotGraph
from analysis.grids import VariantGrid
from analysis.models import AnalysisNode, NodeGraphType, VariantTag, TagNode, AnalysisVariable, AnalysisTemplate, \
    AnalysisTemplateRun, AnalysisLock, Analysis
from analysis.models.enums import SNPMatrix, MinimisationResultType, NodeStatus, TrioSample
from analysis.models.mutational_signatures import MutationalSignature
from analysis.models.nodes import node_utils
from analysis.models.nodes.analysis_node import NodeVCFFilter, AnalysisClassification, NodeTask
from analysis.models.nodes.node_counts import get_node_count_colors, get_node_counts_mine_and_available
from analysis.models.nodes.node_types import get_node_types_hash
from analysis.models.nodes.sources.cohort_node import CohortNodeZygosityFiltersCollection, CohortNodeZygosityFilter
from analysis.serializers import AnalysisNodeSerializer
from analysis.views.analysis_permissions import get_analysis_or_404, get_node_subclass_or_404, \
    get_node_subclass_or_non_fatal_exception
from analysis.views.nodes.node_view import NodeView
from annotation.models.models import MutationalSignatureInfo
from classification.views.views import create_classification_object, CreateClassificationForVariantView
from library import pandas_utils
from library.constants import WEEK_SECS, HOUR_SECS
from library.database_utils import run_sql
from library.django_utils import add_save_message, get_field_counts, set_form_read_only
from library.guardian_utils import is_superuser
from library.utils import full_class_name, defaultdict_to_dict
from pedigree.models import Pedigree
from snpdb.graphs import graphcache
from snpdb.models import UserSettings, Sample, \
    Cohort, CohortSample, ImportStatus, VCF, get_igv_data, Trio, Variant, GenomeBuild
from variantgrid.celery import app
import numpy as np
import pandas as pd


def analysis_list(request):
    form = forms.CreateAnalysisForm(request.POST or None, user=request.user)
    if request.method == "POST":
        if form.is_valid():
            analysis = form.save()
            return redirect(f"/analysis/{analysis.id}")
        add_save_message(request, False, "Analysis")

    context = {"create_analysis_form": form,
               "analysis_choice_form": AnalysisChoiceForm()}
    return render(request, 'analysis/analyses.html', context)


def analysis_templates(request):
    form = forms.CreateAnalysisTemplateForm(request.POST or None, user=request.user)
    if request.method == "POST":
        if form.is_valid():
            analysis_template = form.save()
            return redirect(analysis_template.analysis)

    context = {"create_analysis_template_form": form,
               "analysis_template_choice_form": AnalysisTemplateTypeChoiceForm()}
    return render(request, 'analysis/analysis_templates.html', context)


def get_analysis_settings(user, analysis):
    user_settings = UserSettings.get_for_user(user)
    igv_data = get_igv_data(user, genome_build=analysis.genome_build)
    if analysis.canonical_transcript_collection:
        canonical_transcript_collection = str(analysis.canonical_transcript_collection)
    else:
        canonical_transcript_collection = ""

    analysis_settings = {
        "annotation_version": analysis.annotation_version_id,
        "node_count_types": analysis.get_node_count_types(),
        "canonical_transcript_collection": canonical_transcript_collection,
        "show_igv_links": analysis.show_igv_links,
        "igv_data": igv_data,
        "open_variant_details_in_new_window": user_settings.variant_link_in_analysis_opens_new_tab,
        "genome_build": str(analysis.genome_build),
    }
    return analysis_settings


def view_analysis(request, analysis_id, active_node_id=0):
    analysis = get_analysis_or_404(request.user, analysis_id)

    nodes = analysis.analysisnode_set.filter(visible=True).select_subclasses()
    node_classes_kwargs = {}
    if analysis.lock_input_sources:
        node_classes_kwargs = {"source_nodes": False}

    node_help_dict = {}
    user_settings = UserSettings.get_for_user(request.user)
    if user_settings.tool_tips:
        node_help_dict = {label: subclass.get_help_text() for label, subclass in get_node_types_hash().items()}

    analysis_variables = [[av.node_id, av.field] for av in AnalysisVariable.objects.filter(node__analysis=analysis)]
    analysis_tags_node = TagNode.get_analysis_tags_node(analysis)

    context = {"node_classes_form": forms.AnalysisNodeClassesForm(**node_classes_kwargs),
               "nodes": nodes,
               "node_count_colors": get_node_count_colors("color"),
               "analysis": analysis,
               "analysis_settings": get_analysis_settings(request.user, analysis),
               "analysis_tags_node": analysis_tags_node,
               "active_node_id": active_node_id,
               "node_help": node_help_dict,
               "analysis_variables": analysis_variables,
               "has_write_permission": analysis.can_write(request.user),
               "warnings": analysis.get_toolbar_warnings(request.user),
               "ANALYSIS_DUAL_SCREEN_MODE_FEATURE_ENABLED": settings.ANALYSIS_DUAL_SCREEN_MODE_FEATURE_ENABLED}
    return render(request, 'analysis/analysis.html', context)


def view_active_node(analysis, active_node=None):
    kwargs = {"analysis_id": analysis.pk}
    url = reverse('analysis', kwargs=kwargs)
    if active_node:
        kwargs["active_node_id"] = active_node.pk
        url = reverse('analysis_node', kwargs=kwargs)
    return redirect(url)


@require_POST
def create_analysis_from_template(request, genome_build_name):
    data = request.POST.dict()
    tag_uuid = data.pop("tag_uuid")
    analysis_template_key = f"{tag_uuid}-analysis_template"
    analysis_template_name = data.pop(analysis_template_key)
    analysis_template = AnalysisTemplate.get_for_user(request.user, analysis_template_name)

    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    template_run = AnalysisTemplateRun.create(analysis_template, genome_build, user=request.user)
    template_run.populate_arguments(data)
    populate_analysis_from_template_run(template_run)

    return view_active_node(template_run.analysis, None)


def trio_wizard(request, cohort_id, sample1_id, sample2_id, sample3_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    if cohort.import_status != ImportStatus.SUCCESS:
        import_status = cohort.get_import_status_display()
        msg = f"Can't create analysis for {cohort} of status {import_status}"
        raise PermissionDenied(msg)

    sample_1 = Sample.get_for_user(request.user, sample1_id)
    sample_2 = Sample.get_for_user(request.user, sample2_id)
    sample_3 = Sample.get_for_user(request.user, sample3_id)

    patient_description_results = []
    for sample in [sample_1, sample_2, sample_3]:
        description = ''
        results = []
        try:
            description = sample.patient.phenotype
            results = sample.patient.patient_text_phenotype.phenotype_description.get_results()
        except:
            pass
        patient_description_results.append([description, results])

    form = UserTrioWizardForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            mother_affected = form.cleaned_data['mother_affected']
            father_affected = form.cleaned_data['father_affected']
            sample_1_person = form.cleaned_data['sample_1']
            sample_2_person = form.cleaned_data['sample_2']
            sample_3_person = form.cleaned_data['sample_3']

            samples = [sample_1, sample_2, sample_3]
            people = [sample_1_person, sample_2_person, sample_3_person]

            mother_cs = None
            father_cs = None
            proband_cs = None

            def get_cohort_sample(sample):
                return cohort.cohortsample_set.get(sample=sample)

            for s, p in zip(samples, people):
                if p == TrioSample.FATHER:
                    father_cs = get_cohort_sample(s)
                elif p == TrioSample.MOTHER:
                    mother_cs = get_cohort_sample(s)
                elif p == TrioSample.PROBAND:
                    proband_cs = get_cohort_sample(s)

            sample_names = "/".join((mother_cs.name, father_cs.name, proband_cs.name))
            trio_name = f"{sample_names} from {cohort}"
            trio, _ = Trio.objects.get_or_create(cohort=cohort,
                                                 mother=mother_cs,
                                                 mother_affected=mother_affected,
                                                 father=father_cs,
                                                 father_affected=father_affected,
                                                 proband=proband_cs,
                                                 defaults={"name": trio_name})
            return redirect(trio)

    context = {"cohort": cohort,
               "sample_1": sample_1,
               "sample_2": sample_2,
               "sample_3": sample_3,
               "form": form,
               "patient_description_results": patient_description_results}
    return render(request, 'analysis/trio_wizard.html', context)


def analysis_editor_and_grid(request, analysis_id, stand_alone=False):
    analysis = get_analysis_or_404(request.user, analysis_id)
    context = {"select_grid_column_form": SelectGridColumnForm(),
               "stand_alone": stand_alone,
               "analysis": analysis}
    return render(request, 'analysis/analysis_editor_and_grid.html', context)


def stand_alone_analysis_editor_and_grid(request, analysis_id):
    return analysis_editor_and_grid(request, analysis_id, True)


def analysis_templates_list(request, pk):
    analysis_template = AnalysisTemplate.get_for_user(request.user, pk)
    context = {"analysis_template": analysis_template,
               "analysis_template_versions": analysis_template.analysistemplateversion_set.order_by("-pk")}
    return render(request, 'analysis/analysis_templates_list.html', context)


def get_node_views_by_class():
    node_views_by_class = {}
    module = importlib.import_module("analysis.views.nodes")
    for _, klass in inspect.getmembers(module, inspect.isclass):
        if issubclass(klass, NodeView):
            if klass.model:
                node_views_by_class[klass.model] = klass.as_view()
    return node_views_by_class


NODE_DISPATCHER = get_node_views_by_class()


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_view(request, analysis_version, node_id, node_version, extra_filters):
    """ So we don't fill up urls with lots of different views, just come here and dispatch
        to subclasses of NodeView in analysis.views.nodes based on the model field """
    node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    view = NODE_DISPATCHER[node.__class__]
    return view(request, pk=node_id, version_id=node_version, extra_filters=extra_filters)


def get_node_sql(grid):
    # Temporarily disabling SQL formatting as it's really slow.
    request = RequestFactory().get('/fake')
    grid_sql = grid.get_sql_params_and_columns(request)[0]
    # grid_sql = sqlparse.format(grid_sql, reindent=True, keyword_case='upper')

    grid.fields = ['id']
    node_sql_ = grid.get_sql_params_and_columns(request)[0]
    # node_sql_ = sqlparse.format(node_sql_, reindent=True, keyword_case='upper')

    return node_sql_, grid_sql


@not_minified_response
@cache_page(WEEK_SECS)
@vary_on_cookie
def node_debug(request, analysis_version, node_id, node_version, extra_filters):
    node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    model_name = node._meta.label
    node_serializers = AnalysisNodeSerializer.get_node_serializers()
    serializer_klass = node_serializers.get(model_name, AnalysisNodeSerializer)
    serializer = serializer_klass(node, context={"request": request})

    context = {"node": node,
               "node_data": dict(sorted(serializer.data.items()))}
    if node.valid:
        grid = VariantGrid(request.user, node, extra_filters)
        try:
            node_sql_, grid_sql = get_node_sql(grid)
            context['node_sql'] = node_sql_
            context['grid_sql'] = grid_sql
        except EmptyResultSet:
            pass
    return render(request, "analysis/node_editors/grid_editor_debug_tab.html", context)


def node_doc(request, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    has_write_permission = node.analysis.can_write(request.user)
    form = forms.AnalysisNodeForm(request.POST or None, instance=node)
    if not has_write_permission:
        set_form_read_only(form)

    if request.method == "POST":
        node.analysis.check_can_write(request.user)
        if form.is_valid():
            # Doesn't set "queryset_dirty" so won't cause expensive reloads
            node = form.save()

    context = {"form": form,
               "node": node,
               "has_write_permission": has_write_permission}
    return render(request, "analysis/node_editors/grid_editor_doc_tab.html", context)


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_data_grid(request, analysis_version, node_id, node_version, extra_filters):
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        return HttpResponseRedirect(reverse("node_load", kwargs={"node_id": node_id}))

    context = {"analysis_version": analysis_version,
               "node_id": node_id,
               "node_version": node_version,
               "extra_filters": extra_filters,
               "bams_dict": node.get_bams_dict()}
    return render(request, 'analysis/node_data/node_data_grid.html', context)


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_column_summary(request, analysis_version, node_id, node_version, extra_filters, grid_column_name, significant_figures):
    node = get_node_subclass_or_404(request.user, node_id, version=node_version)

    grid = VariantGrid(request.user, node, extra_filters)
    cm = grid.get_column_colmodel(grid_column_name)
    variant_column = cm['name']
    sorttype = cm.get('sorttype')
    quantitative = sorttype in ['float', 'int']

    context = {"node_id": node.id,
               "node_version": node.version,
               "extra_filters": extra_filters}

    if quantitative:
        label = cm['label']

        poll_url = reverse(column_summary_boxplot, kwargs={"node_id": node_id, "label": label, "variant_column": variant_column})
        context["poll_url"] = poll_url
        template = 'analysis/node_data/node_data_graph.html'
        return render(request, template, context)

    context['grid_column_name'] = grid_column_name
    context['significant_figures'] = int(significant_figures)
    template = 'analysis/node_data/node_data_column_summary_grid.html'
    return render(request, template, context)


def get_snp_matrix_counts(user: User, node_id, version_id):
    node = get_node_subclass_or_404(user, node_id, version=version_id)
    qs = node.get_queryset().filter(locus__ref__length=1, alt__length=1)
    count_qs = qs.values_list('locus__ref__seq', 'alt__seq').distinct().annotate(Count('id'))

    bases = list('ACGT')
    counts_df = pd.DataFrame(index=bases, columns=bases, dtype='i').fillna(0)

    for ref, alt, count in count_qs:
        if alt == Variant.REFERENCE_ALT:
            alt = ref
        counts_df[ref][alt] = count
    return counts_df


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_snp_matrix(request, node_id, node_version, conversion, significant_figures):
    counts_df = get_snp_matrix_counts(request.user, node_id, node_version)
    total = counts_df.sum().sum()
    ti = tv = ti_tv_ratio = None
    if total:
        ti = 0
        ti += counts_df['A']['G']
        ti += counts_df['G']['A']
        ti += counts_df['C']['T']
        ti += counts_df['T']['C']
        tv = total - ti
        ti_tv_ratio = ti / tv

    df = None
    if conversion == SNPMatrix.TOTAL_PERCENT:
        df = pandas_utils.get_total_percent_dataframe(counts_df)
    elif conversion == SNPMatrix.ROWS_PERCENT:
        df = pandas_utils.get_rows_percent_dataframe(counts_df)
    elif conversion == SNPMatrix.COLS_PERCENT:
        df = pandas_utils.get_columns_percent_dataframe(counts_df)

    conversion_description = SNPMatrix(conversion).label

    context = {"node_id": node_id,
               "node_version": node_version,
               'counts_df': counts_df,
               'conversion_description': conversion_description,
               'other_df': df,
               'significant_figures': int(significant_figures),
               'ti': ti,
               'tv': tv,
               'ti_tv_ratio': ti_tv_ratio,
               "extra_filters": None}
    template = 'analysis/node_data/node_data_snp_matrix.html'
    return render(request, template, context)


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_data_graph(request, analysis_version, node_id, node_version, graph_type_id, cmap):
    context = {"node_id": node_id,
               "node_version": node_version,
               "extra_filters": None}

    node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    poll_url = reverse(node_graph, kwargs={"node_id": node.id, "graph_type_id": graph_type_id, "cmap": cmap})
    context["poll_url"] = poll_url
    template = 'analysis/node_data/node_data_graph.html'
    return render(request, template, context)


@cache_page(HOUR_SECS)
@vary_on_cookie
def node_async_wait(request, analysis_version, node_id, node_version, extra_filters):
    node = get_node_subclass_or_404(request.user, node_id)

    context = {"analysis_version": analysis_version,
               "node_id": node_id,
               "node_version": node_version,
               "node": node,
               "extra_filters": extra_filters}

    template = 'analysis/node_data/node_async_wait.html'
    return render(request, template, context)


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_errors(request, analysis_version, node_id, node_version, extra_filters):
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        return HttpResponseRedirect(reverse("node_load", kwargs={"node_id": node_id}))

    context = {"analysis_version": analysis_version,
               "node_id": node_id,
               "node_version": node_version,
               "node": node,
               "errors": node.get_errors(flat=True),
               "extra_filters": extra_filters,
               "status": node.status}

    template = 'analysis/node_data/node_errors.html'
    return render(request, template, context)


def node_load(request, node_id):
    """ loads main grid area, and triggers loading of node editor
        we can't cache this as the node version may have changed since last visit """

    try:
        errors = []
        node = get_node_subclass_or_non_fatal_exception(request.user, node_id)
        extra_filters = request.GET.get("extra_filters", "default")
        kwargs = {"analysis_version": node.analysis.version,
                  "node_id": node.id,
                  "node_version": node.version,
                  "extra_filters": extra_filters}
        if NodeStatus.is_error(node.status):
            view_name = "node_errors"
        elif NodeStatus.is_ready(node.status):
            view_name = "node_data_grid"
        else:
            view_name = "node_async_wait"
        url = reverse(view_name, kwargs=kwargs)
        if node.cloned_from:
            url += f"?cloned_from={node.cloned_from}"
        return HttpResponseRedirect(url)
    except NonFatalNodeError as e:
        errors = [str(e)]

    context = {"errors": errors}
    template = 'analysis/node_data/node_errors.html'
    return render(request, template, context)


@require_POST
def node_cancel_load(request, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    if node_task := NodeTask.objects.filter(node=node, version=node.version).first():
        if node_task.celery_task:
            logging.debug("TODO: Cancelling task %s", node_task.celery_task)
            app.control.revoke(node_task.celery_task, terminate=True)  # @UndefinedVariable

            result = AbortableAsyncResult(node_task.celery_task)
            result.abort()

        if node_task.db_pid:
            run_sql("select pg_cancel_backend(%s)", [node_task.db_pid])
    else:
        logging.error("No task set for node %s", node_id)

    node.status = NodeStatus.CANCELLED
    node.save()

    return HttpResponse()


def node_graph(request, node_id, graph_type_id, cmap):
    get_node_subclass_or_404(request.user, node_id)  # Permission check
    node_graph_type = NodeGraphType.objects.get(pk=graph_type_id)
    cached_graph = graphcache.async_graph(node_graph_type.graph_class, cmap, node_id)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def column_summary_boxplot(request, node_id, label, variant_column):
    get_node_subclass_or_404(request.user, node_id)  # Permission check
    graph_class_name = full_class_name(ColumnBoxplotGraph)
    cached_graph = graphcache.async_graph(graph_class_name, node_id, label, variant_column)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def cohort_zygosity_filters(request, cohort_node_id, cohort_id):
    """ Called from Cohort Node editor - cohort is what's currently selected,
        not what is saved in node """
    cohort_node = get_node_subclass_or_404(request.user, cohort_node_id)
    cohort = Cohort.get_for_user(request.user, cohort_id)

    cnzfc = CohortNodeZygosityFiltersCollection.get_for_node_and_cohort(cohort_node, cohort)

    class CohortSampleNameReadOnlyTextInput(forms.TextInput):
        input_type = 'text'

        def render(self, name, value, attrs=None, renderer=None):
            cs = CohortSample.objects.get(pk=value)
            return cs.sample.name

    CNZFFormSet = inlineformset_factory(CohortNodeZygosityFiltersCollection,
                                        CohortNodeZygosityFilter,
                                        can_delete=False,
                                        fields=['cohort_sample', 'show_in_grid', 'zygosity_ref', 'zygosity_het', 'zygosity_hom', 'zygosity_none'],
                                        widgets={'cohort_sample': CohortSampleNameReadOnlyTextInput(attrs={'readonly':'readonly'})},
                                        extra=0)

    formset = CNZFFormSet(request.POST or None, instance=cnzfc)
    context = {'formset': formset}

    template = 'analysis/node_editors/cohort_zygosity_filters.html'
    return render(request, template, context)


def vcf_locus_filters(request, node_id, vcf_id):
    node = get_node_subclass_or_404(request.user, node_id)
    if vcf_id:
        vcf = VCF.get_for_user(request.user, vcf_id)
    else:
        vcf = None

    context = {"vcf": vcf}
    if vcf:
        vcf_filter_descriptions = {"PASS": "All filters passed"}
        set_filters = {}
        for npf in NodeVCFFilter.filter_for_node(node, vcf):
            if npf.vcf_filter:
                filter_id = npf.vcf_filter.filter_id
            else:
                filter_id = "PASS"
            set_filters[filter_id] = True
        existing_filter_settings = {"PASS": "PASS" in set_filters}

        for vcf_filter in vcf.vcffilter_set.all():
            filter_id = vcf_filter.filter_id
            vcf_filter_descriptions[filter_id] = vcf_filter.description
            existing_filter_settings[filter_id] = filter_id in set_filters

        context["has_filters"] = vcf.has_filters
        context["has_filters_set"] = bool(set_filters)
        context["vlf_form"] = VCFLocusFilterForm(vcf_filters=existing_filter_settings)
        context["vlf_descriptions"] = vcf_filter_descriptions

    return render(request, 'analysis/node_editors/vcf_locus_filters.html', context)


def sample_vcf_locus_filters(request, node_id, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    return vcf_locus_filters(request, node_id, sample.vcf.pk)


def cohort_vcf_locus_filters(request, node_id, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.vcf:
        vcf_id = cohort.vcf.pk
    else:
        vcf_id = None

    return vcf_locus_filters(request, node_id, vcf_id)


def pedigree_vcf_locus_filters(request, node_id, pedigree_id):
    pedigree = Pedigree.get_for_user(request.user, pedigree_id)
    if pedigree.cohort.vcf:
        vcf_id = pedigree.cohort.vcf.pk
    else:
        vcf_id = None

    return vcf_locus_filters(request, node_id, vcf_id)


def view_analysis_settings(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    analysis_settings = get_analysis_settings(request.user, analysis)

    form = forms.CreateAnalysisTemplateForm(request.POST or None, user=request.user, analysis=analysis)
    if request.method == "POST":
        if form.is_valid():
            analysis_template = form.save()
            return JsonResponse({"analysis_id": analysis_template.analysis_id})

    context = {"analysis": analysis,
               "create_analysis_template_form": form,
               "new_analysis_settings": analysis_settings,
               "has_write_permission": analysis.can_write(request.user),
               "can_unlock": analysis.can_unlock(request.user)}
    return render(request, 'analysis/analysis_settings.html', context)


def analysis_settings_details_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    old_annotation_version = analysis.annotation_version
    form = forms.AnalysisForm(request.POST or None, user=request.user, instance=analysis)
    has_write_permission = analysis.can_write(request.user)
    if not has_write_permission:
        set_form_read_only(form)

    reload_analysis = False
    if request.method == "POST":
        analysis.check_can_write(request.user)
        if form.has_changed:
            valid = form.is_valid()
            if valid:
                analysis = form.save()
                analysis.save()
                if reload_analysis := (old_annotation_version != analysis.annotation_version):
                    node_utils.reload_analysis_nodes(analysis.pk)

            add_save_message(request, valid, "Analysis Settings")

    for error in analysis.get_errors():
        messages.add_message(request, messages.ERROR, error)

    for warning in analysis.get_warnings():
        messages.add_message(request, messages.WARNING, warning)

    analysis_settings = get_analysis_settings(request.user, analysis)
    context = {"analysis": analysis,
               "form": form,
               "new_analysis_settings": analysis_settings,
               "has_write_permission": has_write_permission,
               "reload_analysis": reload_analysis}
    return render(request, 'analysis/analysis_settings_details_tab.html', context)


def analysis_settings_node_counts_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    return _analysis_settings_node_counts_tab(request, analysis, has_write_permission=analysis.can_write(request.user))


def _analysis_settings_node_counts_tab(request, analysis, pass_analysis_settings=True, has_write_permission=True):
    if request.method == "POST":
        if has_write_permission is False:
            raise PermissionDenied()
        node_counts_str = request.POST.get("node_counts")
        if node_counts_str is not None:
            node_counts_array = node_counts_str.split(',')
            analysis.set_node_count_types(node_counts_array)
        add_save_message(request, True, "Node Counts")

    my_node_counts_list, available_node_counts_list = get_node_counts_mine_and_available(analysis)
    context = {"my_node_counts_list": my_node_counts_list,
               "available_node_counts_list": available_node_counts_list,
               "has_write_permission": has_write_permission}

    if pass_analysis_settings:
        analysis_settings = get_analysis_settings(request.user, analysis)
        context["new_analysis_settings"] = analysis_settings

    return render(request, 'analysis/analysis_settings_node_counts_tab.html', context)


def analysis_settings_template_run_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)

    node_variables = defaultdict(list)
    for node in analysis.analysisnode_set.filter(analysisvariable__isnull=False).distinct().order_by("y"):
        for av in node.analysisvariable_set.all().order_by("field"):
            node_variables[node].append(av)

    context = {"analysis_template_run": analysis.analysistemplaterun,
               "node_variables": defaultdict_to_dict(node_variables)}
    return render(request, 'analysis/analysis_settings_template_run_tab.html', context)


@require_POST
def analysis_settings_lock(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    if not analysis.can_unlock(request.user):  # check_can_write returns false if locked
        raise PermissionDenied(f"You do not have write access to {analysis.pk}")
    lock = json.loads(request.POST["lock"])
    AnalysisLock.objects.create(analysis=analysis, locked=lock, user=request.user, date=timezone.now())
    # Bump version to expire cache
    analysis.version += 1
    analysis.save()
    return redirect(analysis)  # Reload


def analysis_input_samples(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    input_nodes = analysis.analysisnode_set.filter(analysisnode_parent__isnull=True)
    input_nodes = input_nodes.select_subclasses()

    node_inputs = []
    for node in input_nodes:
        if samples := node.get_samples_from_node_only_not_ancestors():
            node_input_data = {
                "class": node.get_class_name(),
                "name": node.name,
                "samples": samples,
            }
            for field in ["trio", "pedigree", "cohort"]:
                try:
                    node_input_data[field] = getattr(node, field)
                except AttributeError:
                    pass

            node_inputs.append(node_input_data)

    context = {"analysis": analysis,
               "node_inputs": node_inputs}
    return render(request, 'analysis/analysis_input_samples.html', context)


def node_method_description(request, node_id, node_version):
    node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    nodes = AnalysisNode.depth_first(node)

    context = {"node": node,
               "nodes": nodes}
    return render(request, 'analysis/node_method_description.html', context)


@user_passes_test(is_superuser)
def view_analysis_issues(request):
    all_nodes = AnalysisNode.objects.all()
    field_counts = get_field_counts(all_nodes, "status")
    summary_data = Counter()
    for status, count in field_counts.items():
        summary = NodeStatus.get_summary_state(status)
        summary_data[summary] += count

    field_counts = {NodeStatus(k).label: v for k, v in field_counts.items()}
    context = {"nodes_status_summary": summary_data,
               "field_counts": field_counts}
    return render(request, 'analysis/view_analysis_issues.html', context)


def view_mutational_signature(request, pk):
    mutational_signature = get_object_or_404(MutationalSignature, pk=pk)
    # Ensure you can see this sample - sample
    mutational_signature.sample.check_can_view(request.user)
    results_qs = mutational_signature.mutationalsignatureminimisationresult_set.all()
    bootstrapped_results_qs = results_qs.filter(result_type=MinimisationResultType.BOOTSTRAPPED)

    sorted_index = np.flip(np.argsort(mutational_signature.mean), axis=0)
    bootstrapped_data = []
    for r in bootstrapped_results_qs.order_by("iteration"):
        bootstrapped_data.append(r.solution_array)
    sorted_bootstrap_transposed = np.array(bootstrapped_data).T[sorted_index]
    sorted_bootstrap_std = sorted_bootstrap_transposed.std(axis=1)
    sorted_mean = np.array(mutational_signature.mean)[sorted_index]
    labels = [f"Sig {i}" for i in sorted_index + 1]  # signatures are 1 based
    signature_info = list(MutationalSignatureInfo.objects.all().order_by("signature_id"))
    sorted_signature_info = [signature_info[i] for i in sorted_index]
    sorted_data = zip(labels, sorted_mean, sorted_bootstrap_std, sorted_signature_info)

    # Data Fit plot
    mut_count_qs = mutational_signature.mutationalsignaturemutationcount_set.all().order_by("pk")
    mutation_types = []
    mutation_type_counts = []

    for mutation_type, count in mut_count_qs.values_list("mutation_type", "count"):
        mutation_types.append(mutation_type)
        mutation_type_counts.append(count)

    context = {"mutational_signature": mutational_signature,
               "labels": labels,
               "sorted_bootstrap_transposed": sorted_bootstrap_transposed.tolist(),
               "mutation_types": mutation_types,
               "mutation_type_counts": mutation_type_counts,
               "sorted_data": sorted_data}
    return render(request, 'analysis/view_mutational_signature.html', context)


class CreateClassificationForVariantTagView(CreateClassificationForVariantView):
    template_name = "analysis/create_classification_for_variant_tag.html"

    def _get_variant(self):
        return self.variant_tag.variant

    def _get_genome_build(self) -> GenomeBuild:
        return self.variant_tag.genome_build

    def _get_form_post_url(self) -> str:
        if self.variant_tag.analysis:
            return reverse("create_classification_for_analysis", kwargs={"analysis_id": self.variant_tag.analysis.pk})
        else:
            # Just for variant
            return super()._get_form_post_url()

    def _get_sample_form(self):
        # If we have a node with input samples, use that. Then fall back on all samples in analysis.
        # Otherwise fall back on default (all samples in DB visible to user)
        samples = None
        if self.variant_tag.analysis:
            if self.variant_tag.node:
                samples = self.variant_tag.node.get_subclass().get_samples()

            if not samples:
                samples = self.variant_tag.analysis.get_samples()

        if samples:
            form = InputSamplesForm(samples=samples)
            form.fields['sample'].required = False
        else:
            form = super()._get_sample_form()
        return form

    @lazy
    def variant_tag(self):
        variant_tag_id = self.kwargs["variant_tag_id"]
        variant_tag = VariantTag.get_for_user(self.request.user, pk=variant_tag_id)
        return variant_tag

    def get_context_data(self, *args, **kwargs):
        try:
            context = super().get_context_data(*args, **kwargs)
            context["variant_tag"] = self.variant_tag

            if not self.variant_tag.can_write(self.request.user):
                if self.variant_tag.analysis:
                    read_only_message = "You have read-only access to this analysis. You can create a " \
                                        "classification but it will not be linked to the analysis and the " \
                                        f"{settings.TAG_REQUIRES_CLASSIFICATION} tag will not be deleted."
                else:
                    read_only_message = "You have read-only access to this tag. You can create a classification " \
                                        f"but the {settings.TAG_REQUIRES_CLASSIFICATION} tag will not be deleted."
                messages.add_message(self.request, messages.WARNING, read_only_message)
        except VariantTag.DoesNotExist:
            variant_tag_id = self.kwargs["variant_tag_id"]
            msg = f"The VariantTag ({variant_tag_id}) does not exist. It may have been deleted or already classified."
            context = {"error_message": msg}
        return context


@require_POST
def create_classification_for_analysis(request, analysis_id):
    classification = create_classification_object(request)
    analysis = Analysis.get_for_user(request.user, pk=analysis_id)

    if analysis.can_write(request.user):
        AnalysisClassification.objects.create(analysis=analysis, classification=classification)

        # Remove "Requires classification" tag
        VariantTag.objects.filter(variant=classification.variant, analysis=analysis,
                                  tag_id=settings.TAG_REQUIRES_CLASSIFICATION).delete()
    return redirect(classification)
