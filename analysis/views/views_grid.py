import logging
import time
from urllib.parse import urlencode

from django.core.cache import cache
from django.http.response import HttpResponseRedirect
from django.shortcuts import redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from analysis import grids
from analysis.models import AnalysisNode
from analysis.tasks.analysis_grid_export_tasks import (
    NODE_EXPORT_GENERATOR,
    export_cohort_to_downloadable_file,
    export_node_to_downloadable_file,
    export_sample_to_downloadable_file,
    get_grid_downloadable_file_params_hash,
    get_node_grid_downloadable_file_params_hash,
)
from analysis.exceptions import NonFatalNodeError
from analysis.views.analysis_permissions import (
    get_node_subclass_or_404,
    get_node_subclass_or_non_fatal_exception,
)
from analysis.views.node_json_view import NodeJSONGetView, NodeJSONViewMixin
from library.constants import WEEK_SECS
from library.django_utils.jqgrid_datatable_adapter import datatable_data, datatable_definition
from library.django_utils.major_operation import TooManyMajorOperationsError, major_operation
from library.utils.hash_utils import sha256sum_str
from snpdb.models import CachedGeneratedFile, Cohort, Sample

EXPORT_TYPES = {"csv", "vcf"}

# What a node grid request may carry, so the cache key and the export params hash are built off a
# known set. The DataTables client sends start/length/order (@see translate_datatable_params turning
# them into rows/page/sidx/sord); the engine's own rows/page/sidx/sord names are here too because
# node_grid_export builds its params from this list and the export path reads them off the request.
_NODE_GRID_ALLOWED_PARAMS = {
    '_search',
    'ccc_id',
    'ccc_version_id',
    'extra_filters',
    'filters',
    'length',
    'node_id',
    'order[0][column]',
    'order[0][dir]',
    'page',
    'rows',
    'sidx',
    'sord',
    'start',
    'version_id',
    'zygosity_samples_hash',
}


def _add_allowed_node_grid_params(url: str, params: dict) -> str:
    cleaned_params = {}
    for key, value in params.items():
        if key in _NODE_GRID_ALLOWED_PARAMS:
            cleaned_params[key] = value
        else:
            logging.warning("Node redirect had disallowed GET param: %s", key)
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
        try:
            node = self._get_node(request)
        except NonFatalNodeError:
            # Let get_response raise it again and format the non-fatal JSON the client expects
            return self.get_response(request, *args, **kwargs)
        url = reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id})
        url = _add_allowed_node_grid_params(url, request.GET.dict())
        lock_id = sha256sum_str(f"{url}_{request.user}")
        if cache.add(lock_id, "true", LOCK_EXPIRE):  # Acquire lock
            try:
                logging.info("Got the lock...")
                # Cap concurrent expensive queries per-user (distinct nodes bypass the per-node lock above)
                with major_operation(request.user, "node_grid"):
                    response = self.get_response(request, *args, **kwargs)
            except TooManyMajorOperationsError:
                logging.info("Too many major operations - going to sleep then retry...")
                time.sleep(2)
                response = HttpResponseRedirect(url)
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
        # Don't build queryset if invalid (stale q-dict cache)
        if errors := node.get_errors(flat=True):
            return {"errors": errors, "data": [], "recordsTotal": 0, "recordsFiltered": 0}
        grid = _variant_grid_from_request(request, node)
        return datatable_data(request, grid)


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
            # Always deferred: whether to fetch rows at all is the page's call (the row count
            # placeholder, and the grid tab being hidden) - node_data_grid.html triggers the first load
            # The fields/operations go down for the FilterNode editor's builder, but the grid
            # itself gets no "Filter grid..." button - FilterNode is how you filter an analysis
            ret = datatable_definition(grid, defer_loading=True, filter_builder_toolbar=False)
            # The per-request state the data endpoint needs, which the page sends as its ajax params
            ret["postData"] = grid.extra_config["postData"]
        return ret


def _check_export_type(export_type: str):
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")


def _source_grid_export(request, model, obj_id: int, export_type: str, export_task, generator_name: str):
    """ Launches (or joins) a Celery export of a source node's whole grid and redirects to the
        CachedGeneratedFile poll URL. generator_name is stored against the cached file, so it stays put """
    model.get_for_user(request.user, obj_id)  # Permission check
    _check_export_type(export_type)

    params_hash = get_grid_downloadable_file_params_hash(obj_id, export_type)
    task = export_task.si(obj_id, export_type)
    cgf = CachedGeneratedFile.get_or_create_and_launch(generator_name, params_hash, task)
    if cgf.exception:
        raise ValueError(cgf.exception)
    return redirect(cgf)


def cohort_grid_export(request, cohort_id, export_type):
    return _source_grid_export(request, Cohort, cohort_id, export_type,
                               export_cohort_to_downloadable_file, "export_cohort_to_downloadable_file")


def sample_grid_export(request, sample_id, export_type):
    return _source_grid_export(request, Sample, sample_id, export_type,
                               export_sample_to_downloadable_file, "export_sample_to_downloadable_file")


def node_grid_export(request, analysis_id):
    """ Launches (or joins) a Celery export and redirects to the CachedGeneratedFile poll URL, so a big
        node isn't bound by the gunicorn request timeout (issue #1257) """
    export_type = request.GET["export_type"]
    _check_export_type(export_type)

    # Always export the latest version (node.version below) - deliberately don't check the
    # client's version_id, which goes stale whenever the node is updated, eg a tag delete (#789)
    node = get_node_subclass_or_404(request.user, request.GET["node_id"])
    canonical_transcript_collection_id = None
    if request.GET.get("use_canonical_transcripts"):
        # Whether to use it or not is set server-side. Just use client to see what they wanted
        if ctc := node.analysis.canonical_transcript_collection:
            canonical_transcript_collection_id = ctc.pk
        else:
            logging.warning("Grid request had 'use_canonical_transcripts' but analysis did not.")

    grid_params = {k: v for k, v in request.GET.items() if k in _NODE_GRID_ALLOWED_PARAMS}
    # node.version in task_args is authoritative - the client's (possibly stale) version_id would
    # otherwise split the params_hash so identical exports wouldn't share a cached file
    grid_params.pop("version_id", None)
    task_args = (node.pk, node.version, request.user.pk, export_type, canonical_transcript_collection_id, grid_params)
    params_hash = get_node_grid_downloadable_file_params_hash(*task_args)
    task = export_node_to_downloadable_file.si(*task_args)
    cgf = CachedGeneratedFile.get_or_create_and_launch(NODE_EXPORT_GENERATOR, params_hash, task)
    # A failed export reports its exception through the poll URL, so the client can offer a retry
    return redirect(cgf)


def _node_from_request(request) -> AnalysisNode:
    node_id = request.GET["node_id"]
    version_id = int(request.GET["version_id"])
    return get_node_subclass_or_non_fatal_exception(request.user, node_id, version=version_id)


def _variant_grid_from_request(request, node: AnalysisNode, **kwargs):
    extra_filters = request.GET.get("extra_filters")
    return grids.VariantGrid(request.user, node, extra_filters, **kwargs)
