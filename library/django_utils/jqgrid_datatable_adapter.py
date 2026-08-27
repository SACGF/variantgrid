""" Serves the jqGrid server classes (library.jqgrid) to the standard DataTables client.

The grid classes stay the engine: column building, querysets, sorting, counting, filtering and the
server side formatters shared with the CSV/VCF export path are all unchanged. This module only
translates protocols - a DataTables request in, a jqGrid request out; a jqGrid data envelope back,
a DataTables one out - so a converted page runs the same grid class the exports do.

@see claude/plans/variant_grid_to_datatables_plan.md
"""
import json
import logging
from typing import Any, Optional

from django.core.exceptions import PermissionDenied
from django.http import HttpRequest, HttpResponse, QueryDict
from django.urls.base import resolve, reverse
from django.views.generic.base import View

from library.django_utils.jqgrid_view import JQGridViewOp, grid_export_request
from library.jqgrid.jqgrid import json_encode
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.utils import JsonObjType, nice_class_name
from snpdb.views.datatable_view import DatabaseTableView

logger = logging.getLogger(__name__)

# jqGrid colmodel 'formatter' name -> DataTables client renderer
# @see variantgrid/static_files/default_static/js/variantgrid_formats.js
JQGRID_FORMATTER_TO_CLIENT_RENDERER = {
    "clinvarLink": "VariantGridFormat.clinvarLink",
    "cosmicLink": "VariantGridFormat.cosmicLink",
    "detailsLink": "VariantGridFormat.detailsLink",
    "formatClinGenAlleleId": "VariantGridFormat.clinGenAlleleId",
    "formatDBSNP": "VariantGridFormat.dbsnp",
    "formatMasterMindMMID3": "VariantGridFormat.masterMind",
    "formatMavedbUrnLinks": "VariantGridFormat.mavedbUrn",
    "formatOntologyTerms": "VariantGridFormat.ontologyTerms",
    "formatPubMed": "VariantGridFormat.pubMed",
    "geneSymbolLink": "VariantGridFormat.geneSymbolLink",
    "geneSymbolNewWindowLink": "VariantGridFormat.geneSymbolNewWindowLink",
    "gnomadFilteredFormatter": "VariantGridFormat.gnomadFiltered",
    "linkFormatter": "VariantGridFormat.link",
    "omimLink": "VariantGridFormat.omimLink",
    "tagsFormatter": "VariantGridFormat.tags",
    "tagsGlobalFormatter": "VariantGridFormat.tagsGlobal",
    "unitAsPercentFormatter": "VariantGridFormat.unitAsPercent",
}


def datatable_columns_from_colmodels(colmodels: list[dict]) -> list[JsonObjType]:
    """ One DataTables column per colmodel, in colmodel order - hidden columns included so a client
        column index round-trips to the colmodel it came from (@see JqGridDatatableView.translate_params) """
    columns = []
    for cm in colmodels:
        name = cm["name"]
        css_classes = [c for c in [cm.get("classes"), f"dt-{name}"] if c]
        column: JsonObjType = {
            "data": name,
            "label": cm.get("label", name),
            "orderable": bool(cm.get("sortable", True)),
            "visible": not cm.get("hidden", False),
            "className": " ".join(css_classes),
            "render": JQGRID_FORMATTER_TO_CLIENT_RENDERER.get(cm.get("formatter")),
        }
        if width := cm.get("width"):
            column["width"] = width
        if header_title := cm.get("headerTitle"):
            column["headerTitle"] = header_title
        if formatter_kwargs := cm.get("formatter_kwargs"):
            column["renderKwargs"] = formatter_kwargs
        columns.append(column)
    return columns


def datatable_order_from_config(config: dict, colmodels: list[dict]) -> Optional[list]:
    """ jqGrid sortname/sortorder -> DataTables 'order'. None where the grid has no sortable
        initial column (sortname defaults to 'pk', which is not a grid column) """
    sortname = config.get("sortname")
    if not sortname:
        return None
    for i, cm in enumerate(colmodels):
        if sortname in (cm.get("name"), cm.get("index")):
            if not cm.get("sortable", True):
                return None
            return [[i, config.get("sortorder", "asc")]]
    return None


class JqGridDatatableView(View):
    """ Drop-in replacement for JQGridView on grids whose page has been converted to DataTables.
        Same constructor conventions (URL kwargs + an 'extra_filters' JSON param), same URL shape -
        the 'op' kwarg keeps working so 'download' still streams the server side CSV.

        Add to urls.py like:
            perm_path('all_variants/grid/<genome_build_name>/<slug:op>/',
                      JqGridDatatableView.as_view(grid=AllVariantsGrid, csv_download=True),
                      name='all_variants_grid'),
    """

    grid: Optional[type[Any]] = None  # JqGrid subclass, usually JqGridUserRowConfig
    csv_download = False
    # Rows per page persists to UserGridConfig under the grid's caption, and the client sends stable,
    # minimal params so a cached data response (@see NodeGridHandler) keeps its cache key
    cache_stable_params = True
    scroll_x = True
    defer_loading = False

    @staticmethod
    def create_grid_from_request(request, grid_klass, **kwargs):
        kwargs["user"] = request.user
        # extra_filters are commonly used to filter by something custom via JS
        if extra_filters := request.GET.get("extra_filters"):
            kwargs["extra_filters"] = json.loads(extra_filters)
        return grid_klass(**kwargs)

    def _load_grid(self, request, **kwargs):
        if self.grid is None:
            msg = f"{nice_class_name(self)}.grid not set"
            raise ValueError(msg)
        return self.create_grid_from_request(request, self.grid, **kwargs)

    def get(self, request, *args, **kwargs):
        op = kwargs.get("op")
        grid_kwargs = {k: v for k, v in kwargs.items() if k != "op"}
        grid = self._load_grid(request, **grid_kwargs)

        if op == JQGridViewOp.DOWNLOAD:
            if not self.csv_download:
                msg = f"CSV download requested but 'csv_download' not set on {nice_class_name(self)}"
                raise PermissionDenied(msg)
            return grid_export_request(request, grid, self._csv_name(grid))

        if request.GET.get("dataTableDefinition"):
            return self._json_response(self.json_definition(request, grid, **kwargs))

        return self._json_response(self.get_data(request, grid))

    @staticmethod
    def _json_response(data: JsonObjType) -> HttpResponse:
        # json_encode is strict (allow_nan=False) - NaN in a float annotation column would otherwise
        # go down as a bare NaN token and blow up JSON.parse in the browser
        return HttpResponse(json_encode(data), content_type="application/json")  # pylint: disable=http-response-with-content-type-json

    def post(self, request, *args, **kwargs):
        """ DataTables can be configured to POST - the grid engine only ever reads request.GET """
        request.GET = request.POST
        return self.get(request, *args, **kwargs)

    @staticmethod
    def _csv_name(grid) -> str:
        try:
            return grid.csv_name
        except AttributeError:
            try:
                return nice_class_name(grid.model)
            except AttributeError:
                return nice_class_name(grid)

    def _download_url(self, request, **kwargs) -> Optional[str]:
        if not self.csv_download:
            return None
        if "op" not in kwargs:
            return None  # URL has no op segment, so there's nothing to reverse a download onto
        url_name = resolve(request.path_info).url_name
        return reverse(url_name, kwargs={**kwargs, "op": JQGridViewOp.DOWNLOAD})

    def json_definition(self, request: HttpRequest, grid, **kwargs) -> JsonObjType:
        """ The table definition DataTableDefinition builds the table from. Computed per request -
            hide_non_admin columns and UserGridConfig rows are both per user """
        config = grid.get_config(as_json=False)
        colmodels = grid.get_colmodels(remove_server_side_only=True)

        data: JsonObjType = {
            "columns": datatable_columns_from_colmodels(colmodels),
            "order": datatable_order_from_config(config, colmodels),
            "scrollX": self.scroll_x,
            "searchBoxEnabled": False,
            "downloadCsvButtonEnabled": False,  # server side streaming download instead, see downloadUrl
            "downloadUrl": self._download_url(request, **kwargs),
            "csvName": self._csv_name(grid),
            "ajaxType": "GET",  # keeps @cache_page on the data endpoint working
            "cacheStableParams": self.cache_stable_params,
            "deferLoading": self.defer_loading,
            "approximateCount": True,
            "extra": grid.get_datatable_extra(),
        }
        if isinstance(grid, JqGridUserRowConfig):
            # Rows per page persists per user - the DataTables client POSTs set_user_row_config on change
            data["gridName"] = grid.get_caption()
            data["pageLength"] = config["rowNum"]
            data["lengthMenu"] = config["rowList"]
        return data

    def translate_params(self, request: HttpRequest, grid) -> QueryDict:
        """ DataTables request params -> the jqGrid ones the grid engine reads off request.GET.
            Column filtering ('filters'/'_search') and any page specific params pass through untouched. """
        params = request.GET.copy()

        length = int(params.get("length") or 0)
        start = int(params.get("start") or 0)
        params["rows"] = length  # 0 means no pagination
        params["page"] = (start // length) + 1 if length else 1

        if (column_index := params.get("order[0][column]")) not in (None, ""):
            colmodels = grid.get_colmodels()
            try:
                cm = colmodels[int(column_index)]
            except (ValueError, IndexError):
                logger.warning("DataTables order column %s is not in %s's colmodel",
                               column_index, nice_class_name(grid))
            else:
                # 'index' is the server side column name. Genotype columns pack their sort info into it
                # (e.g. 'cohortgenotype_134:1:samples_zygosity') and VariantGrid._sort_items decodes it
                params["sidx"] = cm.get("index", cm["name"])
                params["sord"] = params.get("order[0][dir]", "asc")
        return params

    def get_data(self, request: HttpRequest, grid) -> JsonObjType:
        original_get = request.GET
        request.GET = self.translate_params(request, grid)
        try:
            grid_data = grid.get_data(request)
        finally:
            request.GET = original_get

        records = grid_data["records"]
        data: JsonObjType = {
            "recordsTotal": records,
            "recordsFiltered": records,
            "data": [self.prepare_row(row) for row in grid_data["rows"]],
        }
        # An estimate rather than a COUNT(*) - the pager shows it as "~N" (@see AllVariantsGrid)
        if approximate_records := grid_data.get("approximate_records"):
            data["approximateRecords"] = approximate_records
        # The client strips 'draw' so the request URL stays cacheable, and restores it on the response.
        # Echoing a 0 here would look stale to DataTables and the draw would be discarded
        if draw := request.GET.get("draw"):
            data["draw"] = int(draw)
        return data

    @staticmethod
    def prepare_row(row: dict) -> JsonObjType:
        """ Rows keep their .values() field name keys, so they line up with the columns' 'data' keys """
        return {k: DatabaseTableView.limit_value_size(DatabaseTableView.sanitize_value(v))
                for k, v in row.items()}
