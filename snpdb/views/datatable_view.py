# -*- coding: utf-8 -*-
import enum
import logging
from datetime import datetime
from typing import List, Dict, Optional

from django.db.models import QuerySet
from kombu.utils import json
from lazy import lazy

from library.log_utils import report_exc_info
from snpdb.views.datatable_mixins import JSONResponseView

logger = logging.getLogger(__name__)


class SortOrder(enum.Enum):
    ASC = 'asc'
    DESC = 'desc'


class RichColumn:
    """
    A column to be presented on a DataTable
    """

    def __init__(self,
                 key: str,
                 name: str = None,
                 sort_keys: List[str] = None,
                 label: str = None,
                 orderable: bool = False,
                 enabled: bool = True,
                 renderer=None,
                 default_sort: Optional[SortOrder] = None,
                 client_renderer: Optional[str] = None,
                 detail: bool = False,
                 css_class: str = None):
        """
        :param key: A column name to be retrieved and returned and sorted on
        :param name: A name to be shared between both client and server for this value
        :param sort_keys: If provided, use this array to order_by when ordering by this column
        :param label: The label of the column
        :param orderable: Can this column be sorted on
        :param enabled: Is this column enabled for this environment/user
        :param renderer: Optional server renderer for the value
        :param default_sort: If this column should be the default sort order, provide ascending or descending here
        :param client_renderer: JavaScript function to render the client
        :param detail: If True, the column will be shown in the expand section of the table only (requires responsive)
        :param css_class: css class to apply to the column
        """
        self.key = key
        self.sort_keys = sort_keys
        self.name = name or key
        self.label = label
        if not label:
            self.label = name or key
            self.label = self.label[0].upper() + self.label[1:]
        self.orderable = orderable
        self.renderer = renderer
        self.client_renderer = client_renderer
        self.default_sort = default_sort
        self.enabled = enabled
        self.detail = detail
        self.css_class = css_class

    @property
    def css_classes(self) -> str:
        return ' '.join([css for css in [
            self.css_class,
            'none' if self.detail else 'all',
            'dt-' + self.name.replace(' ', '-')
        ] if css is not None]).strip()

    def sort_string(self, desc: bool) -> List[str]:
        use_keys = self.sort_keys or [self.key]
        sdir = '-' if desc else ''
        return [f'{sdir}{key}' for key in use_keys]

    def __eq__(self, other):
        if isinstance(other, RichColumn):
            return self.name == other.name
        return False


class DatatableConfig:
    """
    This class both determines how the client side table should be defined (via tags)
    and how the server will send data to it via ajax (via BaseDatatableView)

    Most functionality from BaseDatatableView has been slowly migrated into this class
    as the ability to initate this with just a request comes in handy
    """

    rich_columns: List[RichColumn]  # columns for display
    extra_columns: List[str] = []  # bonus columns to retrieve from the database

    def value_columns(self) -> List[str]:
        column_names = [rc.key for rc in self.rich_columns if rc.key]
        all_columns = list(set(column_names + self.extra_columns))
        return all_columns

    def __init__(self, request):
        self.request = request
        self.user = request.user

    @lazy
    def default_sort_order_column(self) -> RichColumn:
        rcs = [rc for rc in self.enabled_columns if rc.default_sort]
        return rcs[0] if rcs else self.enabled_columns[0]

    def column_index(self, rc: RichColumn) -> int:
        return self.enabled_columns.index(rc)

    @lazy
    def enabled_columns(self) -> List[RichColumn]:
        return [rc for rc in self.rich_columns if rc.enabled]

    def get_initial_queryset(self):
        raise NotImplementedError("Need to provide a model or implement get_initial_queryset!")

    def filter_queryset(self, qs: QuerySet) -> QuerySet:
        """
        Override to apply extra GET/POST params to filter what data is returned
        :param qs: The default QuerySet
        :return: A filtered QuerySet
        """
        return qs

    @lazy
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
        return self._querydict.get(param)

    def get_query_json(self, param: str) -> Optional[Dict]:
        """
        like get_query_param but parses the value as JSON
        :param param: the key of the param
        :return: the value of the param as JSON (or None)
        """
        value = self.get_query_param(param)
        if value:
            return json.loads(value)
        return None

    def ordering(self, qs: QuerySet) -> QuerySet:
        """ Get parameters from the request and prepare order by clause
        """
        #'order[0][column]': ['0'], 'order[0][dir]': ['asc']

        sort_by_list = list()
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
            break

        return qs.order_by(*sort_by_list)


class DatatableMixin(object):
    """ JSON data for datatables """
    config: DatatableConfig
    max_display_length = 100

    def get(self, request, *args, **kwargs):
        self.config = self.config(request)
        return super().get(request, *args, **kwargs)

    @property
    def _querydict(self):
        if self.request.method == 'POST':
            return self.request.POST
        else:
            return self.request.GET

    def initialize(self, *args, **kwargs):
        pass
        # can we set config here? how do we get request back out?
        #if not self.config:
        #    raise ValueError('DatatableMixin must set self.config in initialize')

    def render_column(self, row: Dict, column: RichColumn):
        """ Renders a column on a row. column can be given in a module notation eg. document.invoice.type
        """
        if column.renderer:
            return column.renderer(row)
        if column.key:
            value = row.get(column.key)
            if isinstance(value, datetime):
                value = value.timestamp()
            return value
        else:
            return None

    def ordering(self, qs: QuerySet):
        return self.config.ordering(qs)

    def paging(self, qs):
        limit = min(int(self._querydict.get('length', 10)), self.max_display_length)
        start = int(self._querydict.get('start', 0))

        # if pagination is disabled ("paging": false)
        if limit == -1:
            return qs

        offset = start + limit

        return qs[start:offset]

    def get_initial_queryset(self):
        return self.config.get_initial_queryset()

    def get_query_param(self, param: str):
        return self._querydict.get(param)

    def get_query_json(self, param: str):
        value = self.get_query_param(param)
        if value:
            return json.loads(value)
        return None

    def filter_queryset(self, qs):
        return self.config.filter_queryset(qs)

    def prepare_results(self, qs: QuerySet):
        data = []
        # select out all columns but only send down data for enabled columns
        all_columns = self.config.value_columns()

        for row in qs.values(*all_columns):
            # fix me, do server side rendering
            row_json = {}
            for rc in self.config.enabled_columns:
                value = self.render_column(row=row, column=rc)
                row_json[rc.name] = value
            data.append(row_json)
        return data

    def handle_exception(self, e):
        report_exc_info()
        logger.exception(str(e))
        raise e

    def get_context_data(self, *args, **kwargs):
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


class BaseDatatableView(DatatableMixin, JSONResponseView):
    pass
