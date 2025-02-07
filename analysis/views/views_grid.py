import logging
import time
from urllib.parse import urlencode

from django.contrib.postgres.aggregates.general import StringAgg
from django.core.cache import cache
from django.http.response import StreamingHttpResponse, HttpResponseRedirect
from django.shortcuts import redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from analysis import grids
from analysis.grid_export import node_grid_get_export_iterator
from analysis.models import AnalysisNode
from analysis.tasks.analysis_grid_export_tasks import export_cohort_to_downloadable_file, \
    export_sample_to_downloadable_file, get_grid_downloadable_file_params_hash
from analysis.views.analysis_permissions import get_node_subclass_or_non_fatal_exception
from analysis.views.node_json_view import NodeJSONGetView, NodeJSONViewMixin
from library.constants import WEEK_SECS
from snpdb.models import Sample, Cohort, CachedGeneratedFile
from snpdb.models.models_variant import Variant

_NODE_GRID_ALLOWED_PARAMS = {
    '_filters',
    '_search',
    'ccc_id',
    'ccc_version_id',
    'extra_filters',
    'filters',
    'node_id',
    'page',
    'rows',
    'sidx',
    'sord',
    'version_id',
    'zygosity_samples_hash',
}


def _add_allowed_node_grid_params(url: str, params: dict) -> str:
    cleaned_params = {}
    for key, value in params.items():
        if key in _NODE_GRID_ALLOWED_PARAMS:
            cleaned_params[key] = value
        else:
            logging.warning(f"Node redirect had disallowed GET param: %s", key)
    return f"{url}?" + urlencode(cleaned_params)


@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridHandler(NodeJSONViewMixin):
    def get(self, request, *args, **kwargs):
        """ This can be a really expensive operation (ie a few mins)
            And users can sometimes click multiple times, causing the DB to get slow running the same query
            multiple times, interfering with itself - so make a per-user lock, and redirect any further calls
            which should hopefully hit the cache next time
        """
        LOCK_EXPIRE = 60 * 10  # 10 mins
        node = self._get_node(request)
        url = reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id})
        url = _add_allowed_node_grid_params(url, request.GET.dict())
        lock_id = f"{url}_{request.user}"
        if cache.add(lock_id, "true", LOCK_EXPIRE):  # Acquire lock
            try:
                logging.info("Got the lock...")
                response = self.get_response(request, *args, **kwargs)
            finally:
                cache.delete(lock_id)  # release lock
        else:
            logging.info("Don't have lock - going to sleep then retry...")
            time.sleep(2)
            response = HttpResponseRedirect(url)
        return response

    def _get_node(self, request, **kwargs) -> AnalysisNode:
        return _node_from_request(request)

    def _get_redirect(self, request, node):
        """ If node can be represented by another, use that to re-use cache """
        ret = None
        grid_node_id, grid_node_version = node.get_grid_node_id_and_version()
        if grid_node_id != node.pk:
            node_grid_handler = reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id})
            params = request.GET.dict()
            params.update({"node_id": grid_node_id, "version_id": grid_node_version})
            url = _add_allowed_node_grid_params(node_grid_handler, params)
            ret = HttpResponseRedirect(url)
        return ret

    def _get_data(self, request, node, **kwargs):
        grid = _variant_grid_from_request(request, node)
        return grid.get_data(request)


@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridConfig(NodeJSONGetView):
    def _get_node(self, request, **kwargs) -> AnalysisNode:
        return get_node_subclass_or_non_fatal_exception(request.user, kwargs["node_id"], version=kwargs["node_version"])

    def _get_data(self, request, node, **kwargs):
        errors = node.get_errors(flat=True)

        if errors:
            ret = {"errors": errors}
        else:
            grid = grids.VariantGrid(request.user, node, kwargs["extra_filters"])
            ret = grid.get_config(as_json=False)
        return ret


def cohort_grid_export(request, cohort_id, export_type):
    EXPORT_TYPES = {"csv", "vcf"}
    Cohort.get_for_user(request.user, cohort_id)  # Permission check
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")

    params_hash = get_grid_downloadable_file_params_hash(cohort_id, export_type)
    task = export_cohort_to_downloadable_file.si(cohort_id, export_type)
    cgf = CachedGeneratedFile.get_or_create_and_launch("export_cohort_to_downloadable_file", params_hash, task)
    if cgf.exception:
        raise ValueError(cgf.exception)
    return redirect(cgf)


def sample_grid_export(request, sample_id, export_type):
    EXPORT_TYPES = {"csv", "vcf"}
    Sample.get_for_user(request.user, sample_id)  # Permission check
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")

    params_hash = get_grid_downloadable_file_params_hash(sample_id, export_type)
    task = export_sample_to_downloadable_file.si(sample_id, export_type)
    cgf = CachedGeneratedFile.get_or_create_and_launch("export_sample_to_downloadable_file", params_hash, task)
    if cgf.exception:
        raise ValueError(cgf.exception)
    return redirect(cgf)


def node_grid_export(request, analysis_id):
    export_type = request.GET["export_type"]
    use_canonical_transcripts = request.GET.get("use_canonical_transcripts")

    node = _node_from_request(request)
    canonical_transcript_collection = None
    if use_canonical_transcripts:
        # Whether to use it or not is set server-side. Just use client to see what they wanted
        if ctc := node.analysis.canonical_transcript_collection:
            canonical_transcript_collection = ctc
        else:
            logging.warning("Grid request had 'use_canonical_transcripts' but analysis did not.")

    variant_tags_qs = Variant.objects.filter(varianttag__analysis=node.analysis)
    variant_tags_qs = variant_tags_qs.annotate(tags=StringAgg("varianttag__tag", delimiter=', ', distinct=True))
    variant_tags_dict = dict(variant_tags_qs.values_list("id", "tags"))

    filename, file_iterator = node_grid_get_export_iterator(request, node, export_type,
                                                            canonical_transcript_collection, variant_tags_dict)
    return _get_streaming_response(filename, file_iterator)


def _get_streaming_response(filename, file_iterator):
    response = StreamingHttpResponse(file_iterator(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


def _node_from_request(request) -> AnalysisNode:
    node_id = request.GET["node_id"]
    version_id = int(request.GET["version_id"])
    return get_node_subclass_or_non_fatal_exception(request.user, node_id, version=version_id)


def _variant_grid_from_request(request, node: AnalysisNode, **kwargs):
    extra_filters = request.GET.get("extra_filters")
    return grids.VariantGrid(request.user, node, extra_filters, **kwargs)
