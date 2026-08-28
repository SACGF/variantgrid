""" Serves a pandas DataFrame to the DataTables client.

The DataFrame equivalent of snpdb.views.datatable_view - same definition/data/CSV protocol, but the
rows come from get_dataframe() rather than a queryset, and the columns are whatever the DataFrame
turned out to have. Column keys are positional ('c0', 'c1', ...) so column names holding dots (which
DataTables reads as nested object paths) or repeated names can't collide.
"""
import logging
import math
from typing import Any, Optional

import pandas as pd
from django.http import HttpRequest, QueryDict, StreamingHttpResponse

from library.django_utils.grid_export import csv_streaming_response, grid_export_csv
from library.pandas_utils import df_nan_to_none
from library.utils import JsonObjType, nice_class_name
from snpdb.models import UserGridConfig
from snpdb.views.datatable_mixins import JSONResponseView
from snpdb.views.datatable_view import DATATABLE_CSV_PARAM, SortOrder

logger = logging.getLogger(__name__)

INDEX_COLUMN_KEY = "c0"


def _column_key(index: int) -> str:
    return f"c{index}"


class DataFrameDatatableConfig:
    """ Subclasses provide get_dataframe() - the index becomes the first column """
    grid_name: Optional[str] = None
    csv_name: Optional[str] = None
    search_box_enabled = False
    server_csv_download = True
    scroll_x = False

    index_label = "ID"
    index_visible = True
    index_client_renderer: Optional[str] = None
    # Column name to sort on initially - None sorts on the index
    default_sort_column: Optional[str] = None
    default_sort_order = SortOrder.ASC

    def __init__(self, request: HttpRequest):
        self.request = request
        self.user = request.user

    def get_dataframe(self) -> pd.DataFrame:
        raise NotImplementedError()

    def get_column_label(self, column_name: Any) -> str:
        return str(column_name)

    def get_column_client_renderer(self, column_name: Any) -> Optional[str]:
        return None

    def get_column_width(self, column_name: Any) -> Optional[int]:
        return None

    def get_extra(self) -> JsonObjType:
        """ Grid wide metadata handed to the client renderers """
        return {}

    def get_query_param(self, param: str) -> Optional[str]:
        if (resolver_match := self.request.resolver_match) and param in resolver_match.kwargs:
            return resolver_match.kwargs[param]
        return self.request.GET.get(param)


class DataFrameTableView(JSONResponseView):
    """ Wraps a DataFrameDatatableConfig to serve the DataTables client """
    column_class: Optional[type[DataFrameDatatableConfig]] = None
    max_display_length = 100

    config: DataFrameDatatableConfig

    def config_for_request(self, request: HttpRequest) -> DataFrameDatatableConfig:
        return self.column_class(request)

    def get(self, request: HttpRequest, *args, **kwargs):
        self.request = request
        self.config = self.config_for_request(request)
        if self._querydict.get(DATATABLE_CSV_PARAM):
            return self.download_csv()
        return super().get(request, *args, **kwargs)

    @property
    def _querydict(self) -> QueryDict:
        if self.request.method == 'POST':
            return self.request.POST
        return self.request.GET

    def _dataframe(self) -> pd.DataFrame:
        return df_nan_to_none(self.config.get_dataframe())  # JSON can't handle NaN

    def _columns(self, df: pd.DataFrame) -> list[tuple[str, Any]]:
        """ (client column key, DataFrame column name) - the index is the first column """
        columns = [(INDEX_COLUMN_KEY, None)]
        columns += [(_column_key(i), name) for i, name in enumerate(df.columns, start=1)]
        return columns

    def _sort(self, df: pd.DataFrame) -> pd.DataFrame:
        config = self.config
        column_key = self._querydict.get('order[0][column]')
        if column_key is not None:
            ascending = self._querydict.get('order[0][dir]') != 'desc'
            column_index = int(column_key)
        else:
            if config.default_sort_column is not None and config.default_sort_column in df.columns:
                column_index = list(df.columns).index(config.default_sort_column) + 1
            else:
                column_index = 0
            ascending = config.default_sort_order == SortOrder.ASC

        if column_index == 0:
            return df.sort_index(ascending=ascending)
        try:
            column_name = df.columns[column_index - 1]
        except IndexError:
            logger.error("DataFrame grid asked to sort on column %s which it doesn't have", column_index)
            return df
        return df.sort_values(column_name, ascending=ascending)

    def _rows(self, df: pd.DataFrame) -> list[JsonObjType]:
        columns = self._columns(df)
        rows = []
        for index, row in df.iterrows():
            if isinstance(index, float) and math.isnan(index):
                index = None
            data = {INDEX_COLUMN_KEY: index}
            for column_key, column_name in columns[1:]:
                data[column_key] = row[column_name]
            rows.append(data)
        return rows

    def get_context_data(self, *args, **kwargs):
        if self._querydict.get("dataTableDefinition"):
            return self.json_definition()

        df = self._sort(self._dataframe())
        records_total = len(df)

        length = min(int(self._querydict.get('length', 10)), self.max_display_length)
        if length == -1:
            df_page = df
        else:
            start = int(self._querydict.get('start', 0))
            df_page = df[start:start + length]

        return {
            'draw': int(self._querydict.get('draw', 0)),
            'recordsTotal': records_total,
            'recordsFiltered': records_total,
            'data': self._rows(df_page),
        }

    def _csv_name(self) -> str:
        return self.config.csv_name or nice_class_name(self.config)

    def download_csv(self) -> StreamingHttpResponse:
        df = self._sort(self._dataframe())
        colmodels = [{"name": key, "label": label} for key, label in self._column_labels(df)]
        return csv_streaming_response(self._csv_name(), grid_export_csv(colmodels, self._rows(df)))

    def _column_labels(self, df: pd.DataFrame) -> list[tuple[str, str]]:
        config = self.config
        labels = [(INDEX_COLUMN_KEY, config.index_label)]
        for column_key, column_name in self._columns(df)[1:]:
            labels.append((column_key, config.get_column_label(column_name)))
        return labels

    def json_definition(self) -> JsonObjType:
        config = self.config
        df = self._dataframe()

        data: JsonObjType = {
            "responsive": False,
            "searchBoxEnabled": config.search_box_enabled,
            "downloadCsvButtonEnabled": False,
            "csvName": self._csv_name(),
            "expandClientRenderer": None,
            "scrollX": config.scroll_x,
            "extra": config.get_extra(),
        }
        if config.server_csv_download:
            data["downloadUrl"] = f"{self.request.path}?{DATATABLE_CSV_PARAM}=1"
        if grid_name := config.grid_name:
            rows, row_selections = UserGridConfig.get_rows_and_selections(self.request.user, grid_name)
            data["gridName"] = grid_name
            data["pageLength"] = rows
            data["lengthMenu"] = row_selections

        columns: list[JsonObjType] = []
        for i, (column_key, column_name) in enumerate(self._columns(df)):
            is_index = i == 0
            columns.append({
                "data": column_key,
                "label": config.index_label if is_index else config.get_column_label(column_name),
                "render": config.index_client_renderer if is_index else config.get_column_client_renderer(column_name),
                "orderable": True,
                "orderSequence": [SortOrder.ASC.value, SortOrder.DESC.value],
                "className": f"dt-{column_key}",
                "visible": config.index_visible if is_index else True,
                "width": config.get_column_width(column_name),
            })
        data["columns"] = columns
        data["order"] = [[self._default_order_index(df), config.default_sort_order.value]]
        return data

    def _default_order_index(self, df: pd.DataFrame) -> int:
        if (column_name := self.config.default_sort_column) is not None and column_name in df.columns:
            return list(df.columns).index(column_name) + 1
        return 0
