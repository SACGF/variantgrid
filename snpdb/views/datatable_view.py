# -*- coding: utf-8 -*-
import enum
import itertools
import logging
import operator
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property, reduce
from typing import Optional, Any, Callable, Union, TypeVar, Generic, Type, List

from django.contrib.auth.models import User
from django.db import models
from django.db.models import QuerySet, Q, F, OrderBy
from django.http import HttpRequest, QueryDict
from kombu.utils import json

from library.log_utils import report_exc_info
from library.utils import pretty_label, JsonDataType, JsonObjType
from snpdb.views.datatable_mixins import JSONResponseView

logger = logging.getLogger(__name__)


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
                 order_sequence: Optional[List[SortOrder]] = None,
                 client_renderer: Optional[str] = None,
                 client_renderer_td: Optional[str] = None,
                 visible: bool = True,
                 detail: bool = False,
                 css_class: str = None,
                 extra_columns: Optional[list[str]] = None):
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
            else:
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
    rich_columns: list[RichColumn]  # columns for display
    expand_client_renderer: Optional[str] = None  # if provided, will expand rows and render content with this JavaScript method
    scroll_x = False

    def row_css(self, row: CellData) -> Optional[str]:
        """
        Override to provide a CSS class to add to the row
        :param row: The row data to be styled
        :return: A css class or None
        """
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
        if param in self.request.resolver_match.kwargs:
            return self.request.resolver_match.kwargs[param]

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

        sort_by_list: list[OrderBy] = []
        sorted_set = set()
        for index in range(len(self.enabled_columns)):
            order_key = f'order[{index}][column]'
            column_index_str = self.get_query_param(order_key)
            if column_index_str:
                column_index = int(column_index_str)
                sort_order = self.get_query_param(f'order[{index}][dir]')
                rich_column = self.enabled_columns[column_index]

                sorted_set.add(rich_column.name)
                sort_by_list += rich_column.sort_string(sort_order == 'desc')
            else:
                break

        for col in self.rich_columns:
            if col.default_sort:
                if col.name not in sorted_set:
                    sort_by_list += col.sort_string(col.default_sort == SortOrder.DESC)

        sort_by_list.append(self._get_sort_tiebreaker())
        return qs.order_by(*sort_by_list)

    def pre_render(self, qs: QuerySet[DC]):
        """
        Last method called before we start rendering
        qs: The QuerySet with all filtering, ordering applied
        """
        pass

    def view_primary_key(self, row: CellData) -> JsonDataType:
        """ Relies on being 'id' and object defining get_absolute_url  """
        qs = self.get_initial_queryset()
        primary_key_name = qs.model._meta.pk.name
        pk = row.get(primary_key_name)
        if not pk:
            raise ValueError(f"Need to include primary key ('{primary_key_name}') in columns")
        text = row.value
        obj = qs.model(pk=pk)
        return {
            "text": text,
            "url": obj.get_absolute_url(),
        }


class DatabaseTableView(Generic[DC], JSONResponseView):
    """
    Wraps a column_class to give it functionality for a view to provide data to a DataTables view
    """
    config: DatatableConfig
    max_display_length = 100

    column_class: Type[DC] = None

    def config_for_request(self, request: HttpRequest) -> DatatableConfig[DC]:
        return self.column_class(request)

    def get(self, request: HttpRequest, *args, **kwargs):
        self.config = self.config_for_request(request)
        return super().get(request, *args, **kwargs)

    @property
    def _querydict(self) -> QueryDict:
        if self.request.method == 'POST':
            return self.request.POST
        else:
            return self.request.GET

    def initialize(self, *args, **kwargs):
        pass
        # can we set config here? how do we get request back out?
        # if not self.config:
        #    raise ValueError('DatatableMixin must set self.config in initialize')

    @staticmethod
    def sanitize_value(value: Any) -> Any:
        if isinstance(value, datetime):
            value = value.timestamp()
        return value

    @staticmethod
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
                cloned[key] = DatabaseTableView.limit_value_size(sub_value)
            value = cloned
        return value

    def render_cell(self, row: dict, column: RichColumn) -> JsonDataType:
        """ Renders a column on a row. column can be given in a module notation e.g. document.invoice.type """
        data: Any
        if column.renderer:
            render_data = CellData(all_data=row, key=column.key)
            return DatabaseTableView.limit_value_size(column.renderer(render_data))
        elif column.extra_columns:
            data_dict = {}
            for col in column.value_columns:
                data_dict[col] = DatabaseTableView.limit_value_size(DatabaseTableView.sanitize_value(row.get(col)))
            return data_dict

        elif column.key:
            return DatabaseTableView.limit_value_size(DatabaseTableView.sanitize_value(row.get(column.key)))
        else:
            return None

    def ordering(self, qs: QuerySet[DC]):
        return self.config.ordering(qs)

    def paging(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        limit = min(int(self._querydict.get('length', 10)), self.max_display_length)
        start = int(self._querydict.get('start', 0))

        # if pagination is disabled ("paging": false)
        if limit == -1:
            return qs

        offset = start + limit

        return qs[start:offset]

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
        qs = self.config.filter_queryset(qs)
        if qs is not None:
            if (search_text := self.get_query_param('search[value]')) and (search_text := search_text.strip()):
                qs = self.config.power_search(qs, search_text)
            return qs

        raise NotImplementedError("filter_queryset returned None")

    def prepare_results(self, qs: QuerySet[DC]):
        self.config.pre_render(qs)

        data = []
        # select out all columns but only send down data for enabled columns
        all_columns = self.config.value_columns()

        for row in qs.values(*all_columns):
            # fix me, do server side rendering
            row_json = {}
            for rc in self.config.enabled_columns:
                value = self.render_cell(row=row, column=rc)
                row_json[rc.name] = value
            if row_css := self.config.row_css(row):
                row_json["row_css"] = row_css;
            data.append(row_json)
        return data

    def handle_exception(self, e: BaseException):
        report_exc_info()
        logger.exception(str(e))
        raise e

    def json_definition(self) -> JsonObjType:
        config = self.config

        data: JsonObjType = {
            "responsive": any(col.detail for col in config.enabled_columns),
            "searchBoxEnabled": config.search_box_enabled,
            "downloadCsvButtonEnabled": config.download_csv_button_enabled,
            "expandClientRenderer": config.expand_client_renderer,
            "scrollX": config.scroll_x,
        }
        if config.default_sort_order_column:
            data["order"] = [[config.column_index(config.default_sort_order_column), "asc" if config.default_sort_order_column.default_sort != SortOrder.DESC else "desc"]]

        columns: list[JsonObjType] = []
        for rc in config.enabled_columns:
            columns.append({
                "data": rc.name,
                "label": rc.label,
                "render": rc.client_renderer,
                "createdCell": rc.client_renderer_td,
                "orderable": rc.orderable,
                "orderSequence": [x.value for x in rc.order_sequence],
                "className": rc.css_classes,
                "visible": rc.visible,
            })
        data["columns"] = columns

        return data

    def get_context_data(self, *args, **kwargs):
        if definition_request := self.get_query_param("dataTableDefinition"):
            return self.json_definition()

        try:
            self.initialize(*args, **kwargs)

            # prepare initial queryset
            qs = self.get_initial_queryset()

            # store the total number of records (before filtering)
            total_records = qs.count()

            # apply filters
            qs = self.filter_queryset(qs)

            # number of records after filtering
            total_display_records = qs.count()

            # apply ordering
            qs = self.ordering(qs)

            # apply pagintion
            qs = self.paging(qs)

            # prepare output data
            data = self.prepare_results(qs)

            ret = {'draw': int(self._querydict.get('draw', 0)),
                   'recordsTotal': total_records,
                   'recordsFiltered': total_display_records,
                   'data': data
                   }
            return ret
        except Exception as e:
            return self.handle_exception(e)
