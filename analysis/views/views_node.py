import importlib
import inspect
import logging

import pandas as pd
from celery.contrib.abortable import AbortableAsyncResult
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import EmptyResultSet
from django.db.models import Count
from django.forms.models import inlineformset_factory
from django.http.response import HttpResponse, HttpResponseRedirect
from django.shortcuts import redirect, render
from django.urls.base import reverse
from django.views.decorators.cache import cache_page, never_cache
from django.views.decorators.http import require_POST
from django.views.decorators.vary import vary_on_cookie
from htmlmin.decorators import not_minified_response

from analysis import forms
from analysis.exceptions import NodeOutOfDateException, NonFatalNodeError
from analysis.forms import VCFLocusFilterForm
from analysis.graphs.column_boxplot_graph import ColumnBoxplotGraph
from analysis.grids import VariantGrid
from analysis.models import AnalysisNode, NodeGraphType
from analysis.models.enums import NodeStatus, SNPMatrix
from analysis.models.nodes.analysis_node import NodeTask, NodeVCFFilter
from analysis.models.nodes.node_counts import get_extra_filters_count, is_extra_filter
from analysis.models.nodes.sources.cohort_node import (
    CohortNodeZygosityFilter,
    CohortNodeZygosityFiltersCollection,
)
from analysis.serializers import AnalysisNodeSerializer
from analysis.views.analysis_permissions import (
    get_node_subclass_or_404,
    get_node_subclass_or_non_fatal_exception,
)
from analysis.views.nodes.node_view import NodeView
from library import pandas_utils
from library.constants import HOUR_SECS, WEEK_SECS
from library.django_utils import set_form_read_only
from library.utils import full_class_name
from library.utils.database_utils import queryset_to_sql, render_empty_result_set_sql, run_sql
from pedigree.models import Pedigree
from snpdb.graphs import graphcache
from snpdb.models import VCF, Cohort, CohortSample, UserSettings, Variant
from variantgrid.celery import app


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
def node_view(request, analysis_id, analysis_version, node_id, node_version, extra_filters):
    """ So we don't fill up urls with lots of different views, just come here and dispatch
        to subclasses of NodeView in analysis.views.nodes based on the model field

        We use analysis version to be able to expire cache if the custom columns etc change """
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        return HttpResponseRedirect(reverse("node_load",
                                            kwargs={"analysis_id": analysis_id, "node_id": node_id}))
    view = NODE_DISPATCHER[node.__class__]
    return view(request, pk=node_id, version_id=node_version, extra_filters=extra_filters)


def get_node_sql(grid):
    """ (the node's own query, the grid query built on top of it) - the debug tab shows both, so you
        can see what the display columns cost over the rows the node selects """
    # Temporarily disabling SQL formatting as it's really slow.
    # A node that filters on an empty set (eg comp het with no 2-hit genes) would otherwise
    # short circuit with EmptyResultSet - we still want to show the SQL in the debug tab
    with render_empty_result_set_sql():
        qs = grid.get_initial_queryset()
        grid_sql = queryset_to_sql(qs.values(*grid.value_columns()))
        # grid_sql = sqlparse.format(grid_sql, reindent=True, keyword_case='upper')

        node_sql = queryset_to_sql(qs.values("id"))
        # node_sql_ = sqlparse.format(node_sql_, reindent=True, keyword_case='upper')

    return node_sql, grid_sql


@not_minified_response
@cache_page(WEEK_SECS)
@vary_on_cookie
def node_debug(request, analysis_id, analysis_version, node_id, node_version, extra_filters):
    """ We use analysis version to be able to expire cache if the custom columns etc change """
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        return HttpResponseRedirect(reverse("node_load",
                                            kwargs={"analysis_id": analysis_id, "node_id": node_id}))
    model_name = node._meta.label
    node_serializers = AnalysisNodeSerializer.get_node_serializers()
    serializer_klass = node_serializers.get(model_name, AnalysisNodeSerializer)
    serializer = serializer_klass(node, context={"request": request})

    context = {"node": node,
               "node_data": dict(sorted(serializer.data.items()))}
    if node.valid:
        grid = VariantGrid(request, node, extra_filters)
        try:
            node_sql_, grid_sql = get_node_sql(grid)
            context['node_sql'] = node_sql_
            context['grid_sql'] = grid_sql
        except EmptyResultSet:
            pass
    return render(request, "analysis/node_editors/grid_editor_debug_tab.html", context)


@not_minified_response
# @cache_page(WEEK_SECS)
# @vary_on_cookie
def node_audit_log(request, analysis_id, analysis_version, node_id, node_version, extra_filters):
    """ We use analysis version to be able to expire cache if the custom columns etc change """
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        return HttpResponseRedirect(reverse("node_load",
                                            kwargs={"analysis_id": analysis_id, "node_id": node_id}))
    context = {
        "node": node
    }
    return render(request, "analysis/node_editors/grid_editor_audit_log_tab.html", context)


def node_doc(request, analysis_id, node_id):
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
def node_data_grid(request, analysis_id, analysis_version, node_id, node_version, extra_filters):
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        kwargs = {
            "analysis_id": analysis_id,
            "node_id": node_id,
        }
        return HttpResponseRedirect(reverse("node_load", kwargs=kwargs))

    # Use the count for the view being shown - if an extra filter (eg clinvar, a tag selection) is
    # selected, that filtered count (not the whole node's) decides auto-load/sorting.
    # Mirrors VariantGrid._grid_row_count().
    if is_extra_filter(extra_filters):
        grid_row_count = get_extra_filters_count(node, extra_filters)
    else:
        grid_row_count = node.count

    max_variants = (UserSettings.get_for_user(request.user).node_grid_auto_load_max_variants
                    or settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS)
    grid_auto_load = (max_variants is None) or (grid_row_count is not None and grid_row_count < max_variants)

    max_variants_display = None
    if max_variants is not None:
        # eg 50000 -> "50k", 50500 -> "50.5k"
        max_variants_display = f"{max_variants / 1000:g}k" if max_variants >= 1000 else str(max_variants)

    # Sorting is disabled on large grids (see issue #1651) - warn the user so they filter first
    grid_sort_max_variants = settings.ANALYSIS_GRID_SORT_MAX_ROWS
    grid_sorting_disabled = grid_row_count is None or grid_row_count >= grid_sort_max_variants

    context = {
        "analysis_id": analysis_id,
        "analysis_version": analysis_version,
        "node_id": node_id,
        "node_version": node_version,
        "extra_filters": extra_filters,
        "bams_dict": node.get_bams_dict(),
        "node": node,
        "grid_row_count": grid_row_count,
        "grid_auto_load": grid_auto_load,
        "grid_auto_load_max_variants_display": max_variants_display,
        "grid_sorting_disabled": grid_sorting_disabled,
        "grid_sort_max_variants": grid_sort_max_variants,
    }
    return render(request, 'analysis/node_data/node_data_grid.html', context)


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_column_summary(request, analysis_id, analysis_version, node_id, node_version, extra_filters, grid_column_name, significant_figures):
    node = get_node_subclass_or_404(request.user, node_id, version=node_version)

    grid = VariantGrid(request, node, extra_filters)
    rich_column = grid.column(grid_column_name)
    variant_column = rich_column.name
    column_filter = rich_column.column_filter
    quantitative = column_filter is not None and column_filter.type in ['float', 'int']

    context = {
        "analysis_id": analysis_id,
        "node_id": node.id,
        "node_version": node.version,
        "extra_filters": extra_filters
    }

    if quantitative:
        label = rich_column.label

        poll_url = reverse(column_summary_boxplot, kwargs={"analysis_id": analysis_id, "node_id": node_id,
                                                           "label": label, "variant_column": variant_column})
        context["poll_url"] = poll_url
        template = 'analysis/node_data/node_data_graph.html'
        return render(request, template, context)

    context['grid_column_name'] = grid_column_name
    context['significant_figures'] = int(significant_figures)
    template = 'analysis/node_data/node_data_column_summary_grid.html'
    return render(request, template, context)


def get_snp_matrix_counts(user: User, node_id, version_id):
    node = get_node_subclass_or_404(user, node_id, version=version_id)
    qs = node.get_queryset().filter(Variant.get_snp_q())
    count_qs = qs.values_list('locus__ref__seq', 'alt__seq').distinct().annotate(Count('id'))

    bases = list('ACGT')
    counts_df = pd.DataFrame(index=bases, columns=bases, dtype='i').fillna(0)

    for ref, alt, count in count_qs:
        if alt == Variant.REFERENCE_ALT:
            alt = ref
        counts_df.loc[alt, ref] = count
    return counts_df


@cache_page(WEEK_SECS)
@vary_on_cookie
def node_snp_matrix(request, analysis_id, node_id, node_version, conversion, significant_figures):
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
def node_data_graph(request, analysis_id, analysis_version, node_id, node_version, graph_type_id, cmap):
    context = {"node_id": node_id,
               "node_version": node_version,
               "extra_filters": None}

    node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    poll_url = reverse(node_graph, kwargs={"analysis_id": analysis_id, "node_id": node.id,
                                           "graph_type_id": graph_type_id, "cmap": cmap})
    context["poll_url"] = poll_url
    template = 'analysis/node_data/node_data_graph.html'
    return render(request, template, context)


@cache_page(HOUR_SECS)
@vary_on_cookie
def node_async_wait(request, analysis_id, analysis_version, node_id, node_version, extra_filters):
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
def node_errors(request, analysis_id, analysis_version, node_id, node_version, extra_filters):
    try:
        node = get_node_subclass_or_404(request.user, node_id, version=node_version)
    except NodeOutOfDateException:
        return HttpResponseRedirect(reverse("node_load", kwargs={"node_id": node_id}))

    context = {"analysis_id": analysis_id,
               "analysis_version": analysis_version,
               "node_id": node_id,
               "node_version": node_version,
               "node": node,
               "errors": node.get_errors(flat=True),
               "extra_filters": extra_filters,
               "status": node.status}

    template = 'analysis/node_data/node_errors.html'
    return render(request, template, context)


@never_cache
def node_load(request, analysis_id, node_id):
    """ loads main grid area, and triggers loading of node editor
        we can't cache this as the node version may have changed since last visit """

    try:
        errors = []
        node = get_node_subclass_or_non_fatal_exception(request.user, node_id)
        extra_filters = request.GET.get("extra_filters", "default")
        kwargs = {
            "analysis_id": node.analysis_id,
            "analysis_version": node.analysis.version,
            "node_id": node.id,
            "node_version": node.version,
            "extra_filters": extra_filters
        }
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
def node_cancel_load(request, analysis_id, node_id):
    node = get_node_subclass_or_404(request.user, node_id)
    if node_task := NodeTask.objects.filter(node_version__node=node, node_version__version=node.version).first():
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


def node_graph(request, analysis_id, node_id, graph_type_id, cmap):
    """ This is used in node_data_graph """
    get_node_subclass_or_404(request.user, node_id)  # Permission check
    node_graph_type = NodeGraphType.objects.get(pk=graph_type_id)
    cached_graph = graphcache.async_graph(node_graph_type.graph_class, cmap, node_id)
    return redirect(cached_graph)


def column_summary_boxplot(request, analysis_id, node_id, label, variant_column):
    """ This is used in node_column_summary """
    get_node_subclass_or_404(request.user, node_id)  # Permission check
    graph_class_name = full_class_name(ColumnBoxplotGraph)
    cached_graph = graphcache.async_graph(graph_class_name, node_id, label, variant_column)
    return redirect(cached_graph)


def cohort_zygosity_filters(request, analysis_id, node_id, cohort_id):
    """ Called from Cohort Node editor - cohort is what's currently selected,
        not what is saved in node """
    cohort_node = get_node_subclass_or_404(request.user, node_id)
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
    context = {'formset': formset,
               'cohort': cohort}

    template = 'analysis/node_editors/cohort_zygosity_filters.html'
    return render(request, template, context)


def _render_vcf_locus_filters(request, node, vcf):
    """ One VCF's own FILTER codes, plus the node level PASS - only PASS means the same thing in
        every VCF, so a code is never offered from or resolved into another one """
    context = {"vcf": vcf}
    if vcf:
        set_filter_ids = {filter_id for _, filter_id in NodeVCFFilter.get_vcf_filter_ids(node, vcf)}
        pass_set = NodeVCFFilter.has_pass(node)

        vcf_filter_descriptions = {"PASS": "All filters passed"}
        existing_filter_settings = {"PASS": pass_set}
        for vcf_filter in vcf.vcffilter_set.all():
            filter_id = vcf_filter.filter_id
            vcf_filter_descriptions[filter_id] = vcf_filter.description
            existing_filter_settings[filter_id] = filter_id in set_filter_ids

        context["has_filters"] = vcf.has_filters
        context["has_filters_set"] = pass_set or bool(set_filter_ids)
        context["vlf_form"] = VCFLocusFilterForm(vcf_filters=existing_filter_settings)
        context["vlf_descriptions"] = vcf_filter_descriptions

    return render(request, 'analysis/node_editors/vcf_locus_filters.html', context)


def vcf_locus_filters(request, analysis_id, node_id, vcf_id):
    node = get_node_subclass_or_404(request.user, node_id)
    vcf = VCF.get_for_user(node.analysis.user, vcf_id) if vcf_id else None
    return _render_vcf_locus_filters(request, node, vcf)


def cohort_vcf_locus_filters(request, analysis_id, node_id, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.vcf:
        vcf_id = cohort.vcf.pk
    else:
        vcf_id = None

    return vcf_locus_filters(request, analysis_id, node_id, vcf_id)


def pedigree_vcf_locus_filters(request, analysis_id, node_id, pedigree_id):
    pedigree = Pedigree.get_for_user(request.user, pedigree_id)
    if pedigree.cohort.vcf:
        vcf_id = pedigree.cohort.vcf.pk
    else:
        vcf_id = None

    return vcf_locus_filters(request, analysis_id, node_id, vcf_id)


def node_method_description(request, analysis_id, node_id, node_version):
    # Deliberately don't check node_version - always describe the latest version (#789)
    node = get_node_subclass_or_404(request.user, node_id)
    nodes = AnalysisNode.depth_first(node)

    context = {"node": node,
               "nodes": nodes}
    return render(request, 'analysis/node_method_description.html', context)
