import enum
import itertools
import logging
import operator
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property, reduce
from typing import Any, Generic, Optional, TypeVar, Union

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models import F, OrderBy, Q, QuerySet
from django.http import HttpRequest, QueryDict, StreamingHttpResponse
from django.urls import reverse
from kombu.utils import json

from library.django_utils.filter_rules import filter_operations, filter_rules_from_params
from library.django_utils.grid_export import csv_streaming_response, grid_export_csv
from library.log_utils import report_exc_info
from library.utils import JsonDataType, JsonObjType, full_class_name, nice_class_name, pretty_label
from snpdb.models import UserGridConfig
from snpdb.views.datatable_mixins import JSONResponseView

logger = logging.getLogger(__name__)

# Client asks for the server side CSV by adding this to the table's own ajax params
DATATABLE_CSV_PARAM = "dataTableCsv"


class SortOrder(enum.Enum):
    ASC = 'asc'
    DESC = 'desc'


@dataclass(frozen=True)
class CellData:
    """
    Parameter to be passed to server side renders,
    call .value to get the single column, otherwise can inspect columns
    """
    all_data: dict[str, Any]
    key: Optional[str]

    @property
    def value(self):
        if not self.key:
            raise ValueError("Can't call value for RichColumn that didn't specify key")
        return self.all_data[self.key]

    def __getitem__(self, item):
        return self.all_data[item]

    def __contains__(self, item) -> bool:
        return item in self.all_data

    def get_nested_json(self, key, sub_key):
        if data := self.all_data.get(key):
            return data.get(sub_key)

    def get(self, key: Any, default: Optional[Any] = None) -> Any:
        return self.all_data.get(key, default)


class RichColumn:
    """
    A column to be presented on a DataTable
    """

    @staticmethod
    def client_renderer_combine(formatters: list[str]) -> str:
        json_str = json.dumps(formatters)
        return f'TableFormat.combine.bind(null, {json_str}, null)'  # 2nd null is settings

    @staticmethod
    def client_renderer_repeat(settings: Optional[dict[str, Any]]) -> str:
        return f'TableFormat.repeat.bind(null, {settings})'

    @staticmethod
    def choices_client_renderer(choices):
        json_data = {}
        for choice in choices:
            json_data[choice[0]] = choice[1]
        json_str = json.dumps(json_data)
        return f'TableFormat.choices.bind(null, { json_str })'

    def __init__(self,
                 key: Optional[str] = None,
                 name: str = None,
                 sort_keys: list[str] = None,
                 search: Optional[Union[bool, list[str]]] = None,
                 label: str = None,
                 orderable: bool = None,
                 enabled: bool = True,
                 renderer: Optional[Callable[[CellData], JsonDataType]] = None,
                 default_sort: Optional[SortOrder] = None,
                 order_sequence: Optional[list[SortOrder]] = None,
                 client_renderer: Optional[str] = None,
                 client_renderer_td: Optional[str] = None,
                 visible: bool = True,
                 detail: bool = False,
                 css_class: str = None,
                 extra_columns: Optional[list[str]] = None,
                 include_in_csv: Optional[bool] = None,
                 width: Optional[int] = None,
                 header_title: Optional[str] = None,
                 server_side_formatter: Optional[Callable[[dict, str], Any]] = None,
                 data_type: Optional[str] = None,
                 filterable: Optional[bool] = None,
                 filter_choices: Optional[dict] = None,
                 queryset_field: bool = True,
                 sort_annotation: Optional[Any] = None,
                 sort_nulls_first: bool = False):
        """
        #TODO consolidate, orderable, default_sort, sort_order_sequence
        :param key: A column name to be retrieved and returned and sorted on
        :param name: A name to be shared between both client and server for this value
        :param sort_keys: If provided, use this array to order_by when ordering by this column
        :param search: If true and search_box_enabled, sort the key
        :param label: The label of the column
        :param orderable: Can this column be sorted on
        :param enabled: Is this column enabled for this environment/user
        :param renderer: Optional server renderer for the value
        :param default_sort: If this column should be the default sort order, provide ascending or descending here
        :param order_sequence: How the column should sort when clicked on (can change to be DESC first, can only allow ASC etc)
        :param client_renderer: JavaScript function to render the client
        :param client_renderer_td: JavaScript function that maps to DataTables createdCell https://datatables.net/reference/option/columns.createdCell
        :param visible: If false column would be hidden (useful for sending data we don't want to display)
        :param detail: If True, the column will be shown in the expand section of the table only (requires responsive)
        :param css_class: css class to apply to the column
        :param extra_columns: other columns that need to be selected out for the server renderer
        :param include_in_csv: whether the raw value goes in the server side CSV (defaults to having a key)
        :param width: column width in pixels - tables laid out 'table-layout: fixed' need one per column
        :param header_title: tooltip on the column header
        :param server_side_formatter: func(row, field) applied to the raw value before it reaches either
                                      the client or the CSV, so the two can't disagree
        :param data_type: 'int'/'float'/'date' - what the filter builder offers, and whether a column
                          summarises quantitatively
        :param filterable: whether the filter builder offers this column (defaults to having a key)
        :param filter_choices: {value: label} the filter builder offers instead of a free text input
        :param queryset_field: False for columns that aren't selected from the queryset
        :param sort_annotation: Django expression to annotate and order by when sorting on this column,
                                for columns with nothing sortable at their key
        :param sort_nulls_first: sort nulls as the lowest value (ascending puts them first) rather than
                                 always at the bottom
        """
        self.key = key
        self.sort_keys = sort_keys
        if orderable is None:
            orderable = bool(sort_keys) or bool(order_sequence) or bool(default_sort)
        if not order_sequence:
            order_sequence = [SortOrder.ASC, SortOrder.DESC]
        self.search = []
        if (search is None or search is True) and key:
            self.search = [key]
        elif search is True and not key:
            raise ValueError("Cannot have search = True if no key provided, provide list of searchable columns instead")
        elif search:
            self.search = search

        if orderable and not self.key and not self.sort_keys and sort_annotation is None:
            raise ValueError("Cannot create an 'orderable' RichColumn without key, sort_keys or sort_annotation")
        self.name = name or key
        if not self.name:
            raise ValueError("Cannot create a RichColumn without key or name must be provided")
        if "." in self.name:
            # This will be treated as nested objects in JS, which is not what we want to do passing literal strings
            # @see https://datatables.net/reference/option/columns.data#string
            raise ValueError("Cannot create a RichColumn with '.' (dot) in name")
        self.label = label
        if not label:
            self.label = self.name
            self.label = pretty_label(self.label)
        if not key and not renderer and not client_renderer:
            raise ValueError("Cannot create a RichColumn without a key, server or client renderer")

        self.orderable = orderable
        self.order_sequence = order_sequence
        self.renderer = renderer
        self.client_renderer = client_renderer
        self.client_renderer_td = client_renderer_td
        self.default_sort = default_sort
        self.enabled = enabled
        self.detail = detail
        self.visible = visible
        self.css_class = css_class
        self.extra_columns = extra_columns
        self.include_in_csv = bool(key) if include_in_csv is None else include_in_csv
        self.width = width
        self.header_title = header_title
        self.server_side_formatter = server_side_formatter
        self.data_type = data_type
        self.filterable = bool(key) if filterable is None else filterable
        self.filter_choices = filter_choices
        self.queryset_field = queryset_field
        self.sort_annotation = sort_annotation
        self.sort_nulls_first = sort_nulls_first

    @property
    def sort_alias(self) -> str:
        """ Alias sort_annotation is annotated under - the name itself may be a real queryset field """
        return f"{self.name}_sort"

    @property
    def css_classes(self) -> str:
        return ' '.join([css for css in [
            self.css_class,
            'none' if self.detail else 'all',
            'dt-' + self.name.replace(' ', '-')
        ] if css is not None]).strip()

    def sort_key(self, key: str, desc: bool) -> str:
        if key.startswith('-'):
            if desc:
                key = key[1:]
        elif desc:
            key = f'-{key}'
        return key

    # the below seems to break special sort keys
    def sort_string(self, desc: bool) -> list[OrderBy]:
        def as_order_by(key: str):
            use_desc = desc
            if key.startswith('-'):
                key = key[1:]
                use_desc = not use_desc

            if use_desc:
                return F(key).desc(nulls_last=True)
            if self.sort_nulls_first:
                return F(key).asc(nulls_first=True)
            return F(key).asc(nulls_last=True)
        if self.sort_keys:
            use_keys = self.sort_keys
        elif self.sort_annotation is not None:
            use_keys = [self.sort_alias]
        else:
            use_keys = [self.key]
        return [as_order_by(key) for key in use_keys]

    # def sort_string(self, desc: bool) -> list[str]:
    #     use_keys = self.sort_keys or [self.key]
    #     return [self.sort_key(key, desc) for key in use_keys]

    @property
    def value_columns(self) -> list[str]:
        columns = []
        if (key := self.key) and self.queryset_field:
            columns.append(key)
        if self.extra_columns:
            columns += self.extra_columns
        if self.search:
            columns += self.search
        return columns

    def __eq__(self, other):
        if isinstance(other, RichColumn):
            return self.name == other.name
        return False


DC = TypeVar('DC', bound=models.Model)  # Data Class


class DatatableConfig(Generic[DC]):
    """
    This class both determines how the client side table should be defined (via tags)
    and how the server will send data to it via ajax (via BaseDatatableView)
    """
    search_box_enabled = False
    download_csv_button_enabled = False
    # Streams every row's raw values from the server, rather than the client side button which pulls
    # the rendered rows back through the ajax endpoint - use it on anything that can grow large
    server_csv_download = False
    csv_name: Optional[str] = None  # filename for CSV download (date suffix added automatically)
    rich_columns: list[RichColumn]  # columns for display
    expand_client_renderer: Optional[str] = None  # if provided, will expand rows and render content with this JavaScript method
    scroll_x = False
    # Set to opt in to rows per page persisting per user (in UserGridConfig, keyed on this name)
    # rather than in the browser's localStorage
    grid_name: Optional[str] = None
    # DataTables' own cap on rows per page. None serves whatever the client asks for
    max_page_length: Optional[int] = 100
    # Whether the unfiltered row count is worth a query of its own, so the pager can say "filtered
    # from N". Off where that count costs as much as the page it's counting
    count_unfiltered = True
    # CSS class added to the <table>, eg to lay it out 'table-layout: fixed'
    table_class: Optional[str] = None
    # Build the table now and fetch rows when the page asks (table.ajax.reload())
    defer_loading = False
    # Send a minimal, stably ordered param set so a cached data endpoint keeps its key, and strip
    # 'draw' from the request (the client restores it from the response)
    cache_stable_params = False
    # Every pixel of chrome under a big grid is a row the user can't see
    compact_controls = False
    # Show the record count as an estimate ('~1.2M') where the config supplies one
    approximate_count = False
    ajax_type: Optional[str] = None  # 'GET' keeps @cache_page on the data endpoint working
    # The column filter builder - rules come back as a 'filters' param (@see library.django_utils.filter_rules)
    filter_builder = False
    filter_builder_toolbar = True

    def row_css(self, row: CellData) -> Optional[str]:
        """
        Override to provide a CSS class to add to the row
        :param row: The row data to be styled
        :return: A css class or None
        """
        return None

    def csv_columns(self) -> list[dict[str, str]]:
        """ Columns for the server side CSV, as {name, label} - de-duplicated as action columns
            (delete etc) share their key with the column they act on """
        columns = {}
        for rc in self.enabled_columns:
            if rc.include_in_csv and rc.key not in columns:
                columns[rc.key] = {"name": rc.key, "label": rc.label}
        return list(columns.values())

    def row_columns(self) -> list[str]:
        """
        Override to include columns regardless of which rows are rendered
        Typically used to provide data to row_css
        :return: A list of column names
        """
        return []

    def value_columns(self) -> list[str]:
        column_names = list(itertools.chain(*[rc.value_columns for rc in self.rich_columns if rc.enabled]))
        if row_columns := self.row_columns():
            column_names += row_columns
        return list(set(column_names))

    def __init__(self, request: HttpRequest):
        self.request: HttpRequest = request
        self.user: User = request.user
        # Set by get_known_count() where the count it returned was an estimate rather than a COUNT(*)
        self.approximate_records: Optional[str] = None

    def get_csv_name(self) -> str:
        if csv_name := self.csv_name:
            return csv_name
        try:
            return nice_class_name(self.get_initial_queryset().model)
        except Exception:
            return "export"

    def get_extra(self) -> JsonObjType:
        """ Grid wide metadata handed to the client renderers, which get no per column config """
        return {}

    def download_url(self) -> str:
        """ Where the toolbar's CSV button points - the table's own URL, asking for the CSV """
        return f"{self.request.path}?{DATATABLE_CSV_PARAM}=1"

    def get_known_count(self, qs: QuerySet[DC]) -> Optional[int]:
        """ Row count for this request when it's already known (a stored count, a planner estimate),
            so the view can skip the COUNT(*). None means count the queryset """
        return None

    @cached_property
    def filter_rules(self) -> Optional[JsonObjType]:
        """ The filter builder's rules - @see library.django_utils.filter_rules """
        return filter_rules_from_params(self._querydict)

    def row_values(self, qs: QuerySet[DC]) -> Iterable[dict]:
        """ The row dicts the columns are rendered from """
        return qs.values(*self.value_columns())

    def csv_rows(self, qs: QuerySet[DC]) -> Iterable[dict]:
        """ The row dicts the server side CSV is written from """
        return qs.values(*[c["name"] for c in self.csv_columns()]).iterator()

    def filter_builder_fields(self) -> list[JsonObjType]:
        """ The fields the filter builder offers. 'field' goes into a rule and straight into a Django
            lookup, so columns the queryset can't filter on are left out """
        fields = []
        for rc in self.enabled_columns:
            if not (rc.filterable and rc.key):
                continue
            field: JsonObjType = {"field": rc.key, "label": rc.label}
            if rc.filter_choices:
                field["type"] = "select"
                field["choices"] = rc.filter_choices
            else:
                field["type"] = rc.data_type or "text"
            fields.append(field)
        return fields

    @cached_property
    def default_sort_order_column(self) -> RichColumn:
        rcs = [rc for rc in self.enabled_columns if rc.default_sort and rc.visible]
        return rcs[0] if rcs else self.enabled_columns[0]

    def column_index(self, rc: RichColumn) -> int:
        return self.enabled_columns.index(rc)

    @cached_property
    def enabled_columns(self) -> list[RichColumn]:
        return [rc for rc in self.rich_columns if rc.enabled]

    def get_initial_queryset(self) -> QuerySet[DC]:
        raise NotImplementedError("Need to provide get_initial_queryset!")

    def power_search(self, qs: QuerySet[DC], search_string: str) -> QuerySet[DC]:
        search_cols = set()
        rich_col: RichColumn
        for rich_col in self.enabled_columns:  # TODO do we want to check not enabled columns too?
            search_cols = search_cols.union(rich_col.search)

        filters: list[Q] = []
        for search_col in search_cols:
            filters.append(Q(**{f'{search_col}__icontains': search_string}))
        or_filter = reduce(operator.or_, filters)
        qs = qs.filter(or_filter)
        return qs

    def filter_queryset(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        """
        Override to apply extra GET/POST params to filter what data is returned
        :param qs: The default QuerySet
        :return: A filtered QuerySet
        """
        return qs

    @staticmethod
    def _row_expand_ajax(expand_view: str, id_field: str = 'id', expected_height: Optional[int] = None) -> str:
        """
        :expand_view Name of the django view - as found in urls.py
        :id_field Parameter to go to the view (assumes that the view takes 1 parameter)
        :expected_height Expected height in pixels of an expanded view, doesn't need to be accurate, just looks better if it's roughly correct
        """
        expected_height_str = f"{expected_height if expected_height else 100}px"
        return f"TableFormat.expandAjax.bind(null, '{expand_view}', '{id_field}', '{expected_height_str}')"

    @cached_property
    def _querydict(self):
        if self.request.method == 'POST':
            return self.request.POST
        else:
            return self.request.GET

    def get_query_param(self, param: str) -> Optional[str]:
        """
        Returns a param value from the GET or POST
        :param param: the key of the param
        :return: the value of the param
        """
        if (resolver_match := self.request.resolver_match) and param in resolver_match.kwargs:
            return resolver_match.kwargs[param]

        return self._querydict.get(param)

    def get_query_json(self, param: str) -> Optional[dict]:
        """
        like get_query_param but parses the value as JSON
        :param param: the key of the param
        :return: the value of the param as JSON (or None)
        """
        value = self.get_query_param(param)
        if value:
            return json.loads(value)
        return None

    def _get_sort_tiebreaker(self) -> OrderBy:
        """ Ensure we always have a 'tie breaker' and thus consistent sort order for paging.
            May need to overwrite if you use a group by/count in queryset thus no PK """
        return F("pk").desc()

    def ordering(self, qs: QuerySet) -> QuerySet[DC]:
        """ Get parameters from the request and prepare order by clause """
        #  'order[0][column]': ['0'], 'order[0][dir]': ['asc']

        sort_columns: list[tuple[RichColumn, bool]] = []
        sorted_set = set()
        for index in range(len(self.enabled_columns)):
            order_key = f'order[{index}][column]'
            column_index_str = self.get_query_param(order_key)
            if column_index_str:
                column_index = int(column_index_str)
                sort_order = self.get_query_param(f'order[{index}][dir]')
                rich_column = self.enabled_columns[column_index]

                sorted_set.add(rich_column.name)
                sort_columns.append((rich_column, sort_order == 'desc'))
            else:
                break

        for col in self.rich_columns:
            if col.default_sort:
                if col.name not in sorted_set:
                    sort_columns.append((col, col.default_sort == SortOrder.DESC))

        sort_annotations = {rc.sort_alias: rc.sort_annotation
                            for rc, _desc in sort_columns if rc.sort_annotation is not None}
        if sort_annotations:
            qs = qs.annotate(**sort_annotations)

        sort_by_list: list[OrderBy] = []
        for rich_column, desc in sort_columns:
            sort_by_list += rich_column.sort_string(desc)
        sort_by_list.append(self._get_sort_tiebreaker())
        return qs.order_by(*sort_by_list)

    def pre_render(self, qs: QuerySet[DC]):
        """
        Last method called before we start rendering
        qs: The QuerySet with all filtering, ordering applied
        """
        pass

    @cached_property
    def _model(self) -> type[DC]:
        return self.get_initial_queryset().model

    def view_primary_key(self, row: CellData) -> JsonDataType:
        """ Relies on being 'id' and object defining get_absolute_url  """
        primary_key_name = self._model._meta.pk.name
        if primary_key_name not in row:
            raise ValueError(f"Need to include primary key ('{primary_key_name}') in columns")
        pk = row[primary_key_name]
        text = row.value
        if pk is None or pk == "":
            # Legacy data can have blank text primary keys (eg Experiment) - no URL to link to
            return {"text": text}
        obj = self._model(pk=pk)
        return {
            "text": text,
            "url": obj.get_absolute_url(),
        }

    def render_delete(self, cell: CellData) -> Optional[str]:
        try:
            obj = self._model.get_instance_for_permission_check(cell.value)
        except self._model.DoesNotExist:
            return None
        if not obj.can_write(self.user):
            return None
        return reverse('group_permissions_object_delete',
                       kwargs={'class_name': full_class_name(self._model), 'primary_key': cell.value})


def sanitize_value(value: Any) -> Any:
    if isinstance(value, datetime):
        value = value.timestamp()
    return value


def limit_value_size(value: Any) -> Any:
    """
    Limits the amount of data that can be returned in one cell
    Will duplicate dicts into dicts with limited text
    """
    LIMIT = 100000
    if isinstance(value, str):
        if (value_len := len(value)) and value_len > LIMIT:
            return value[:LIMIT] + f"... (data is too large to display, full data is {value_len} characters long)"
    elif isinstance(value, dict):
        cloned = value.copy()
        for key, sub_value in value.items():
            cloned[key] = limit_value_size(sub_value)
        value = cloned
    return value


def render_cell(row: dict, column: RichColumn) -> JsonDataType:
    """ Renders a column on a row. column can be given in a module notation e.g. document.invoice.type """
    if column.renderer:
        render_data = CellData(all_data=row, key=column.key)
        return limit_value_size(column.renderer(render_data))
    if column.extra_columns:
        data_dict = {}
        for col in column.value_columns:
            data_dict[col] = limit_value_size(sanitize_value(row.get(col)))
        return data_dict
    if column.key:
        return limit_value_size(sanitize_value(row.get(column.key)))
    return None


def prepare_rows(config: DatatableConfig[DC], qs: QuerySet[DC]) -> list[JsonObjType]:
    config.pre_render(qs)

    data = []
    for row in config.row_values(qs):
        row_json = {}
        for rc in config.enabled_columns:
            row_json[rc.name] = render_cell(row=row, column=rc)
        if row_css := config.row_css(row):
            row_json["row_css"] = row_css
        data.append(row_json)
    return data


def _paging(config: DatatableConfig[DC], qs: QuerySet[DC]) -> QuerySet[DC]:
    limit = int(config.get_query_param('length') or 10)
    if limit == -1:  # pagination disabled ("paging": false)
        return qs
    if config.max_page_length is not None:
        limit = min(limit, config.max_page_length)
    start = int(config.get_query_param('start') or 0)
    return qs[start:start + limit]


def datatable_data(config: DatatableConfig[DC]) -> JsonObjType:
    """ One page of rows, as the DataTables data envelope """
    qs = config.get_initial_queryset()
    filtered_qs = config.filter_queryset(qs)
    if filtered_qs is None:
        raise NotImplementedError("filter_queryset returned None")
    if (search_text := config.get_query_param('search[value]')) and (search_text := search_text.strip()):
        filtered_qs = config.power_search(filtered_qs, search_text)

    if (known_count := config.get_known_count(filtered_qs)) is not None:
        # A known count is what this config holds - there's no unfiltered total worth a second query
        total_records = records_filtered = known_count
    elif filtered_qs is qs:
        # filter_queryset hands back the same queryset when no filters were supplied, and the count is
        # often expensive enough to be worth not repeating
        total_records = records_filtered = qs.count()
    elif config.count_unfiltered:
        total_records = qs.count()
        records_filtered = filtered_qs.count()
    else:
        total_records = records_filtered = filtered_qs.count()

    page_qs = _paging(config, config.ordering(filtered_qs))
    data: JsonObjType = {
        'recordsTotal': total_records,
        'recordsFiltered': records_filtered,
        'data': prepare_rows(config, page_qs),
    }
    # An estimate rather than a COUNT(*) - the pager shows it as "~N"
    if approximate_records := config.approximate_records:
        data["approximateRecords"] = approximate_records
    # A cache stable client strips 'draw' from the request and restores it on the response. Echoing a
    # 0 here would look stale to DataTables and the draw would be discarded
    if (draw := config.get_query_param('draw')) is not None:
        data['draw'] = int(draw)
    return data


def datatable_definition(config: DatatableConfig[DC]) -> JsonObjType:
    """ The table definition DataTableDefinition builds the table from. Computed per request - the
        columns and the UserGridConfig rows are both per user """
    data: JsonObjType = {
        "responsive": any(col.detail for col in config.enabled_columns),
        "searchBoxEnabled": config.search_box_enabled,
        "downloadCsvButtonEnabled": config.download_csv_button_enabled,
        "csvName": config.get_csv_name(),
        "expandClientRenderer": config.expand_client_renderer,
        "scrollX": config.scroll_x,
        "tableClass": config.table_class,
        "ajaxType": config.ajax_type,
        "cacheStableParams": config.cache_stable_params,
        "deferLoading": config.defer_loading,
        "approximateCount": config.approximate_count,
        "compactControls": config.compact_controls,
        "extra": config.get_extra(),
    }
    if config.server_csv_download:
        data["downloadUrl"] = config.download_url()
    if grid_name := config.grid_name:
        rows, row_selections = UserGridConfig.get_rows_and_selections(config.user, grid_name)
        data["gridName"] = grid_name
        data["pageLength"] = rows
        data["lengthMenu"] = row_selections
    if config.filter_builder:
        data["filterBuilder"] = {
            "fields": config.filter_builder_fields(),
            "operations": filter_operations(),
        }
        data["filterBuilderToolbar"] = config.filter_builder_toolbar

    if default_sort_column := config.default_sort_order_column:
        sort_order = "desc" if default_sort_column.default_sort == SortOrder.DESC else "asc"
        data["order"] = [[config.column_index(default_sort_column), sort_order]]

    columns: list[JsonObjType] = []
    for rc in config.enabled_columns:
        column: JsonObjType = {
            "data": rc.name,
            "label": rc.label,
            "render": rc.client_renderer,
            "createdCell": rc.client_renderer_td,
            "orderable": rc.orderable,
            "orderSequence": [x.value for x in rc.order_sequence],
            "className": rc.css_classes,
            "visible": rc.visible,
        }
        if rc.width:
            column["width"] = f"{rc.width}px"
        if rc.header_title:
            column["headerTitle"] = rc.header_title
        columns.append(column)
    data["columns"] = columns

    return data


def datatable_csv_response(config: DatatableConfig[DC]) -> StreamingHttpResponse:
    """ Streams every row of the config's (filtered, ordered) queryset as CSV """
    if not config.server_csv_download:
        raise PermissionDenied(f"CSV download requested but 'server_csv_download' not set on "
                               f"{nice_class_name(config)}")
    columns = config.csv_columns()
    qs = config.ordering(config.filter_queryset(config.get_initial_queryset()))
    return csv_streaming_response(config.get_csv_name(), grid_export_csv(columns, config.csv_rows(qs)))


class DatabaseTableView(Generic[DC], JSONResponseView):
    """
    Wraps a column_class to give it functionality for a view to provide data to a DataTables view
    """
    config: DatatableConfig

    column_class: type[DC] = None

    def config_for_request(self, request: HttpRequest) -> DatatableConfig[DC]:
        return self.column_class(request)

    def get(self, request: HttpRequest, *args, **kwargs):
        self.request = request
        self.config = self.config_for_request(request)
        if self._querydict.get(DATATABLE_CSV_PARAM):
            return datatable_csv_response(self.config)
        return super().get(request, *args, **kwargs)

    @property
    def _querydict(self) -> QueryDict:
        if self.request.method == 'POST':
            return self.request.POST
        else:
            return self.request.GET

    def initialize(self, *args, **kwargs):
        pass

    def handle_exception(self, e: BaseException):
        report_exc_info()
        logger.exception(str(e))
        raise e

    def get_context_data(self, *args, **kwargs):
        if self._querydict.get("dataTableDefinition"):
            return datatable_definition(self.config)

        try:
            self.initialize(*args, **kwargs)
            return datatable_data(self.config)
        except Exception as e:
            return self.handle_exception(e)
