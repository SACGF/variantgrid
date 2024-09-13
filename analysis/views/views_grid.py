import logging
import time
from urllib.parse import urlencode

from django.contrib.postgres.aggregates.general import StringAgg
from django.core.cache import cache
from django.http.response import StreamingHttpResponse, HttpResponseRedirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from analysis import grids
from analysis.analysis_templates import get_sample_analysis, get_cohort_analysis
from analysis.grid_export import node_grid_get_export_iterator
from analysis.models import AnalysisNode, AnalysisTemplate, SampleNode
from analysis.views.analysis_permissions import get_node_subclass_or_non_fatal_exception
from analysis.views.node_json_view import NodeJSONGetView, NodeJSONViewMixin
from library.constants import WEEK_SECS
from library.utils import name_from_filename
from snpdb.models import Sample, Cohort
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

    cohort = Cohort.get_for_user(request.user, cohort_id)
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")

    analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
    analysis = get_cohort_analysis(cohort, analysis_template)
    node = analysis.analysisnode_set.get_subclass(output_node=True)  # Should only be 1
    basename = "_".join([name_from_filename(cohort.name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(cohort.genome_build)])

    filename, file_iterator = node_grid_get_export_iterator(request, node, export_type, basename=basename)
    return _get_streaming_response(filename, file_iterator)


def sample_grid_export(request, sample_id, export_type):
    EXPORT_TYPES = {"csv", "vcf"}

    sample = Sample.get_for_user(request.user, sample_id)
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")

    analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_SAMPLE")
    analysis = get_sample_analysis(sample, analysis_template)
    node = SampleNode.objects.get(analysis=analysis, output_node=True)  # Should only be 1
    basename = "_".join([name_from_filename(sample.name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(sample.genome_build)])
    filename, file_iterator = node_grid_get_export_iterator(request, node, export_type, basename=basename)
    return _get_streaming_response(filename, file_iterator)


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
