import enum
import itertools
import logging
import math
import operator
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
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

from library.django_utils.filter_rules import filter_operations_json, parse_filters, rules_to_q
from library.django_utils.grid_export import csv_streaming_response, grid_export_csv
from library.django_utils.major_operation import MajorOperationViewMixin
from library.log_utils import report_exc_info
from library.utils import JsonDataType, JsonObjType, full_class_name, nice_class_name, pretty_label
from snpdb.models import AvatarDetails, UserGridConfig, UserSettings
from snpdb.views.datatable_mixins import JSONResponseView

logger = logging.getLogger(__name__)

# Client asks for the server side CSV by adding this to the table's own ajax params
DATATABLE_CSV_PARAM = "dataTableCsv"
# The column filter rules the client sends up (@see library.django_utils.filter_rules)
DATATABLE_FILTERS_PARAM = "filters"


class SortOrder(enum.Enum):
    ASC = 'asc'
    DESC = 'desc'


class NullOrder(enum.Enum):
    """ Where NULLs go when sorting on a column """
    LAST = 'last'  # last whichever direction - the default
    FIRST_ON_ASC = 'first_on_asc'  # first ascending, last descending (what the variant grids do)


@dataclass(frozen=True)
class FilterField:
    """ How the client filter builder should offer a column, and what a rule on it means server side.
        'field' in an emitted rule is the column's key, which goes straight into a Django lookup """
    type: str = 'text'  # 'text' / 'int' / 'float' / 'date' / 'select'
    choices: Optional[dict] = field(default=None)

    def as_json(self, name: str, label: str) -> JsonObjType:
        data: JsonObjType = {"field": name, "label": label, "type": self.type}
        if self.choices:
            data["choices"] = self.choices
        return data


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


def sanitize_value(value: Any) -> Any:
    if isinstance(value, datetime):
        value = value.timestamp()
    elif isinstance(value, float) and math.isnan(value):
        # JSON goes out strict - a NaN in a float annotation column renders blank rather than
        # breaking JSON.parse in the browser
        value = None
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
                 csv_rendered: bool = False,
                 width: Optional[int] = None,
                 header_title: Optional[str] = None,
                 client_renderer_kwargs: Optional[dict] = None,
                 sort_menu: Optional[list[dict]] = None,
                 column_filter: Optional[FilterField] = None,
                 null_order: NullOrder = NullOrder.LAST):
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
        :param include_in_csv: whether the column goes in the server side CSV (defaults to having a key)
        :param csv_rendered: write render_cell() output to the CSV rather than the raw value - set it
                             wherever the server renderer converts the value (AF as percent, choices
                             expanded, packed genotype unpacked) so the CSV matches the grid
        :param width: column width in pixels (the variant grids lay out table-layout: fixed)
        :param header_title: tooltip on the column header
        :param client_renderer_kwargs: per column settings handed to the client renderer
        :param sort_menu: [{label, column}] alternative sort keys offered on a composite cell's header
        :param column_filter: offer this column in the client's filter builder, and accept rules on it
        :param null_order: where NULLs sort
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

        if orderable and not self.key and not self.sort_keys:
            raise ValueError("Cannot create an 'orderable' RichColumn without key or sort_keys")
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
        self.csv_rendered = csv_rendered
        self.width = width
        self.header_title = header_title
        self.client_renderer_kwargs = client_renderer_kwargs
        self.sort_menu = sort_menu
        self.column_filter = column_filter
        self.null_order = null_order

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
        nulls_first_on_asc = self.null_order == NullOrder.FIRST_ON_ASC

        def as_order_by(key: str):
            use_desc = desc
            if key.startswith('-'):
                key = key[1:]
                use_desc = not use_desc

            if use_desc:
                return F(key).desc(nulls_last=True)
            if nulls_first_on_asc:
                return F(key).asc(nulls_first=True)
            return F(key).asc(nulls_last=True)
        use_keys = self.sort_keys or [self.key]
        return [as_order_by(key) for key in use_keys]

    # def sort_string(self, desc: bool) -> list[str]:
    #     use_keys = self.sort_keys or [self.key]
    #     return [self.sort_key(key, desc) for key in use_keys]

    @property
    def value_columns(self) -> list[str]:
        columns = []
        if key := self.key:
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
    # recordsTotal only feeds DataTables' "(filtered from N total)" text. Turn this off where the
    # unfiltered queryset is expensive to count - the filtered count is then reported for both.
    count_unfiltered = True
    # Extra css classes on the <table> - the variant grids lay out table-layout: fixed
    table_class: Optional[str] = None
    # Fills in for any column that doesn't set its own width (px). table-layout: fixed needs one per column
    default_column_width: Optional[int] = None
    # Trim the chrome around the table - every pixel of it is a row the user can't see
    compact_controls = False
    # Build the table but hold the first row request back until the page asks for it
    defer_loading = False
    # Send a minimal, stably ordered param set so identical grid state keeps a cached response's key
    cache_stable_params = False
    ajax_type = 'POST'  # 'GET' keeps @cache_page on a data endpoint working
    # Offer the column filter builder (columns carrying a column_filter), and accept its rules
    filter_builder = False
    # The page can mount its own builder off the definition without the grid's "Filter grid..." button
    filter_builder_toolbar = True
    # Whether the pager may show a "~N" estimate instead of a count (@see approximate_count)
    approximate_count_enabled = False
    # Hovering a row for 500ms fetches its expanded content early
    expand_prefetch = True
    max_page_length = 100
    # CSV header labels that differ from the column's own (@see csv_columns)
    csv_label_overrides: Optional[dict[str, str]] = None

    def row_css(self, row: CellData) -> Optional[str]:
        """
        Override to provide a CSS class to add to the row
        :param row: The row data to be styled
        :return: A css class or None
        """
        return None

    def csv_columns(self) -> list[dict[str, str]]:
        """ Columns for the server side CSV, as {name, label} - one per export column """
        label_overrides = self.csv_label_overrides or {}
        return [{"name": rc.name, "label": label_overrides.get(rc.name, rc.label)}
                for rc in self.export_columns()]

    def get_extra(self) -> JsonObjType:
        """ Grid wide metadata for the client renderers - things that belong to the grid rather than
            to any one column """
        return {}

    def get_table_classes(self) -> list[str]:
        return [self.table_class] if self.table_class else []

    def post_data(self) -> Optional[JsonObjType]:
        """ Per request state the page sends back as the ajax params of every row request """
        return None

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
        self._page_rows: list[dict] = []
        self._page_writable_pks: Optional[set] = None
        self._user_labels: dict[int, str] = {}

    @cached_property
    def viewer_settings(self) -> UserSettings:
        return UserSettings.get_for_user(self.user)

    def render_user(self, cell: CellData) -> JsonDataType:
        """ For a "<fk>__username" column with extra_columns=["<fk>__id"]: renders the name through
            AvatarDetails so a title holder gets their crown (#1819). Sort/search/CSV
            stay on the username """
        user_id = cell.get(cell.key.removesuffix("__username") + "__id")
        if user_id is None:
            return ""
        if (label := self._user_labels.get(user_id)) is None:
            user = User.objects.filter(pk=user_id).first()
            label = str(AvatarDetails.avatar_for(user).grid_label_html(self.viewer_settings)) if user else ""
            self._user_labels[user_id] = label
        return label

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

    @cached_property
    def filter_rules(self) -> Optional[dict]:
        """ The column filter rules this request carries (@see library.django_utils.filter_rules) """
        if not self.filter_builder:
            return None
        return parse_filters(self.get_query_param(DATATABLE_FILTERS_PARAM))

    @property
    def filter_rules_supplied(self) -> bool:
        """ known_count/approximate_count implementations decline when column filters narrow the rows """
        return self.filter_rules is not None

    def filter_fields(self) -> list[JsonObjType]:
        """ The fields the filter builder offers. 'field' goes straight into a Django lookup, so only
            columns that named a column_filter are here """
        fields = []
        for rc in self.enabled_columns:
            if column_filter := rc.column_filter:
                fields.append(column_filter.as_json(rc.name, rc.label))
        return fields

    def apply_filter_rules(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        if rules := self.filter_rules:
            if q := rules_to_q(rules):
                qs = qs.filter(q)
        return qs

    def known_count(self, qs: QuerySet[DC]) -> Optional[int]:
        """ Row count for this request when it's already known (a stored node count, a planner
            estimate), so nothing has to run a COUNT(*). None means count the queryset """
        return None

    def approximate_count(self, qs: QuerySet[DC]) -> Optional[str]:
        """ Called after known_count - the "~N" the pager shows in place of an exact count """
        return None

    def render_cell(self, row: dict, column: RichColumn) -> JsonDataType:
        """ Renders a column on a row """
        if column.renderer:
            return limit_value_size(column.renderer(CellData(all_data=row, key=column.key)))
        if column.extra_columns:
            return {col: limit_value_size(sanitize_value(row.get(col))) for col in column.value_columns}
        if column.key:
            return limit_value_size(sanitize_value(row.get(column.key)))
        return None

    def render_rows(self, rows: Iterable[dict]) -> Iterator[dict]:
        """ Raw .values() rows -> {column name: rendered value} """
        for row in rows:
            row_json = {rc.name: self.render_cell(row, rc) for rc in self.enabled_columns}
            if row_css := self.row_css(row):
                row_json["row_css"] = row_css
            yield row_json

    def export_columns(self) -> list[RichColumn]:
        """ The columns the CSV/VCF export writes. Two columns writing the same raw value collapse to
            one - an action column (delete etc) shares its key with the column it acts on """
        columns = []
        raw_keys = set()
        for rc in self.enabled_columns:
            if not rc.include_in_csv:
                continue
            if not rc.csv_rendered and rc.key is not None:
                if rc.key in raw_keys:
                    continue
                raw_keys.add(rc.key)
            columns.append(rc)
        return columns

    def render_export_rows(self, rows: Iterable[dict],
                           columns: Optional[list[RichColumn]] = None) -> Iterator[dict]:
        """ Rows for the CSV/VCF export - the server renderer's output where the column asked for it
            (csv_rendered), the raw value everywhere else """
        if columns is None:
            columns = self.export_columns()
        for row in rows:
            yield {rc.name: self.render_cell(row, rc) if rc.csv_rendered else sanitize_value(row.get(rc.key))
                   for rc in columns}

    def iter_export_rows(self, qs: QuerySet[DC]) -> Iterator[dict]:
        return self.render_export_rows(qs.values(*self.value_columns()).iterator())

    def apply_filters(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        """ Everything that narrows the rows: the config's own params, the search box, then the
            client's column filter rules """
        filtered_qs = self.filter_queryset(qs)
        if filtered_qs is None:
            raise NotImplementedError("filter_queryset returned None")
        if (search_text := self.get_query_param('search[value]')) and (search_text := search_text.strip()):
            filtered_qs = self.power_search(filtered_qs, search_text)
        return self.apply_filter_rules(filtered_qs)

    def paging(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        limit = min(int(self.get_query_param('length') or 10), self.max_page_length)
        start = int(self.get_query_param('start') or 0)
        if limit == -1:  # if pagination is disabled ("paging": false)
            return qs
        return qs[start:start + limit]

    def get_csv_name(self) -> str:
        if csv_name := self.csv_name:
            return csv_name
        try:
            return nice_class_name(self._model)
        except Exception:
            # The definition names the download before any request filtering, and a config whose
            # queryset needs those params raises here - a generic filename is fine
            return "export"

    def initial_order(self) -> Optional[list]:
        """ The client's initial sort - None leaves the table unsorted """
        if rc := self.default_sort_order_column:
            return [[self.column_index(rc), "desc" if rc.default_sort == SortOrder.DESC else "asc"]]
        return None

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

    def requested_ordering(self) -> list[tuple[RichColumn, bool]]:
        """ (column, descending) for each order[i] the request carries """
        #  'order[0][column]': ['0'], 'order[0][dir]': ['asc']
        ordering = []
        for index in range(len(self.enabled_columns)):
            column_index_str = self.get_query_param(f'order[{index}][column]')
            if not column_index_str:
                break
            try:
                rich_column = self.enabled_columns[int(column_index_str)]
            except (ValueError, IndexError):
                logger.warning("%s: order column '%s' is not one of its columns",
                               nice_class_name(self), column_index_str)
                break
            ordering.append((rich_column, self.get_query_param(f'order[{index}][dir]') == 'desc'))
        return ordering

    def ordering(self, qs: QuerySet) -> QuerySet[DC]:
        """ Get parameters from the request and prepare order by clause """
        sort_by_list: list[OrderBy] = []
        sorted_set = set()
        for rich_column, desc in self.requested_ordering():
            sorted_set.add(rich_column.name)
            sort_by_list += rich_column.sort_string(desc)

        for col in self.rich_columns:
            if col.default_sort:
                if col.name not in sorted_set:
                    sort_by_list += col.sort_string(col.default_sort == SortOrder.DESC)

        sort_by_list.append(self._get_sort_tiebreaker())
        return qs.order_by(*sort_by_list)

    def pre_render(self, qs: QuerySet[DC], rows: list[dict]):
        """
        Last method called before we start rendering
        qs: The QuerySet with all filtering, ordering applied
        rows: The page's raw values - use it to resolve in one query what would otherwise be per-row
        Overrides call super() - render_delete resolves the page's write permissions off these rows
        """
        self._page_rows = rows
        self._page_writable_pks = None

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
        """ The delete link, or None where the user can't write the row """
        if cell.value not in self._writable_pks_for_page(cell.key):
            return None
        return reverse('group_permissions_object_delete',
                       kwargs={'class_name': full_class_name(self._model), 'primary_key': cell.value})

    def _writable_pks_for_page(self, pk_column: str) -> set:
        """ Which of the page's rows the user can write, resolved on the first row that asks -
            a pair of Guardian lookups per row is the single most expensive thing a grid can do """
        if self._page_writable_pks is None:
            pks = {pk for row in self._page_rows if (pk := row.get(pk_column)) is not None}
            writable_qs = self._model.filter_writable_for_user(self.user).filter(pk__in=pks)
            self._page_writable_pks = set(writable_qs.values_list("pk", flat=True))
        return self._page_writable_pks


def datatable_response(config: DatatableConfig, draw: Optional[str] = None) -> JsonObjType:
    """ One page of a config's rows in the DataTables envelope. Shared by DatabaseTableView and the
        analysis node grid handler, which sits on its own view mixin for error/lock handling """
    qs = config.get_initial_queryset()
    filtered_qs = config.apply_filters(qs)

    if (known_count := config.known_count(filtered_qs)) is not None:
        total_display_records = known_count
    else:
        total_display_records = filtered_qs.count()

    # The total before filtering - apply_filters hands back the same queryset when nothing was
    # supplied, and counting the unfiltered queryset can be expensive enough to skip
    if config.count_unfiltered and filtered_qs is not qs and known_count is None:
        total_records = qs.count()
    else:
        total_records = total_display_records

    page_qs = config.paging(config.ordering(filtered_qs))
    rows = list(page_qs.values(*config.value_columns()))
    config.pre_render(page_qs, rows)

    data: JsonObjType = {
        'recordsTotal': total_records,
        'recordsFiltered': total_display_records,
        'data': list(config.render_rows(rows)),
    }
    if approximate_records := config.approximate_count(filtered_qs):
        # An estimate rather than a COUNT(*) - the pager shows it as "~N"
        data["approximateRecords"] = approximate_records
    # The client can strip 'draw' so the request URL stays cacheable, and restore it on the response.
    # Echoing a 0 here would look stale to DataTables and the draw would be discarded
    if draw is not None:
        data['draw'] = int(draw)
    return data


def rich_column_json(rc: RichColumn, default_column_width: Optional[int] = None) -> JsonObjType:
    """ One column of the table definition - what DataTableDefinition builds a column (and its cell
        renderer) from. @see the annotation descriptions page, which draws example cells from these """
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
    if width := (rc.width or default_column_width):
        column["width"] = f"{width}px"
    if rc.header_title:
        column["headerTitle"] = rc.header_title
    if rc.client_renderer_kwargs:
        column["renderKwargs"] = rc.client_renderer_kwargs
    if rc.sort_menu:
        # Alternative sort keys for a composite cell - each names another column whose own
        # definition already carries the sort key. @see DataTableDefinition.setupSortMenus
        column["sortMenu"] = rc.sort_menu
    return column


def datatable_definition(config: DatatableConfig, download_url: Optional[str] = None) -> JsonObjType:
    """ The table definition DataTableDefinition builds the table from. Computed per request -
        column visibility and UserGridConfig rows are both per user """
    data: JsonObjType = {
        "responsive": any(col.detail for col in config.enabled_columns),
        "searchBoxEnabled": config.search_box_enabled,
        "downloadCsvButtonEnabled": config.download_csv_button_enabled,
        "csvName": config.get_csv_name(),
        "expandClientRenderer": config.expand_client_renderer,
        "scrollX": config.scroll_x,
    }
    if download_url:
        data["downloadUrl"] = download_url
    if table_classes := config.get_table_classes():
        data["tableClass"] = " ".join(table_classes)
    if config.compact_controls:
        data["compactControls"] = True
    if config.defer_loading:
        data["deferLoading"] = True
    if config.cache_stable_params:
        data["cacheStableParams"] = True
    if config.ajax_type != 'POST':
        data["ajaxType"] = config.ajax_type
    if config.approximate_count_enabled:
        data["approximateCount"] = True
    if config.expand_client_renderer and not config.expand_prefetch:
        data["expandPrefetch"] = False
    if extra := config.get_extra():
        data["extra"] = extra
    if post_data := config.post_data():
        data["postData"] = post_data
    if grid_name := config.grid_name:
        rows, row_selections = UserGridConfig.get_rows_and_selections(config.user, grid_name)
        data["gridName"] = grid_name
        data["pageLength"] = rows
        data["lengthMenu"] = row_selections
    if config.filter_builder:
        data["filterBuilder"] = {
            "fields": config.filter_fields(),
            "operations": filter_operations_json(),
        }
        data["filterBuilderToolbar"] = config.filter_builder_toolbar

    if (order := config.initial_order()) is not None:
        data["order"] = order

    data["columns"] = [rich_column_json(rc, config.default_column_width) for rc in config.enabled_columns]
    return data


class DatabaseTableView(Generic[DC], MajorOperationViewMixin, JSONResponseView):
    """
    Wraps a column_class to give it functionality for a view to provide data to a DataTables view
    """
    # Strict JSON so a NaN in a float annotation column renders blank rather than going down as a
    # bare NaN token and breaking JSON.parse in the browser
    json_allow_nan = False
    config: DatatableConfig

    column_class: type[DC] = None

    def config_for_request(self, request: HttpRequest, **kwargs) -> DatatableConfig[DC]:
        """ kwargs are the URL kwargs, for a config that needs them to build (node, genome build).
            DatatableConfig.get_query_param already falls back to resolver_match.kwargs """
        return self.column_class(request)

    def get(self, request: HttpRequest, *args, **kwargs):
        self.request = request
        self.config = self.config_for_request(request, **kwargs)
        if self._querydict.get(DATATABLE_CSV_PARAM):
            return self.download_csv()
        return super().get(request, *args, **kwargs)

    def download_csv(self) -> StreamingHttpResponse:
        config = self.config
        if not config.server_csv_download:
            raise PermissionDenied(f"CSV download requested but 'server_csv_download' not set on "
                                   f"{nice_class_name(config)}")
        qs = config.ordering(config.apply_filters(config.get_initial_queryset()))
        return csv_streaming_response(self._csv_name(),
                                      grid_export_csv(config.csv_columns(), config.iter_export_rows(qs)))

    def _csv_name(self) -> str:
        return self.config.get_csv_name()

    @property
    def _querydict(self) -> QueryDict:
        if self.request.method == 'POST':
            return self.request.POST
        else:
            return self.request.GET

    def initialize(self, *args, **kwargs):
        pass

    @staticmethod
    def sanitize_value(value: Any) -> Any:
        return sanitize_value(value)

    @staticmethod
    def limit_value_size(value: Any) -> Any:
        return limit_value_size(value)

    def render_cell(self, row: dict, column: RichColumn) -> JsonDataType:
        return self.config.render_cell(row, column)

    def ordering(self, qs: QuerySet[DC]):
        return self.config.ordering(qs)

    def paging(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        return self.config.paging(qs)

    def get_initial_queryset(self) -> QuerySet[DC]:
        return self.config.get_initial_queryset()

    def get_query_param(self, param: str) -> Any:
        return self._querydict.get(param)

    def get_query_json(self, param: str) -> JsonDataType:
        value = self.get_query_param(param)
        if value:
            return json.loads(value)
        return None

    def filter_queryset(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        return self.config.apply_filters(qs)

    def prepare_results(self, qs: QuerySet[DC]):
        # select out all columns but only send down data for enabled columns
        rows = list(qs.values(*self.config.value_columns()))
        self.config.pre_render(qs, rows)
        return list(self.config.render_rows(rows))

    def handle_exception(self, e: BaseException):
        report_exc_info()
        logger.exception(str(e))
        raise e

    def _download_url(self) -> Optional[str]:
        if not self.config.server_csv_download:
            return None
        return f"{self.request.path}?{DATATABLE_CSV_PARAM}=1"

    def json_definition(self) -> JsonObjType:
        return datatable_definition(self.config, download_url=self._download_url())

    def get_context_data(self, *args, **kwargs):
        if self.get_query_param("dataTableDefinition"):
            return self.json_definition()

        try:
            self.initialize(*args, **kwargs)
            return datatable_response(self.config, draw=self._querydict.get('draw'))
        except Exception as e:
            return self.handle_exception(e)
