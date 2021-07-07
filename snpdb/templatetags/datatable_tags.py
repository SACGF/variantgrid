from typing import Union

from django.template import Library
from snpdb.views.datatable_view import DatatableConfig, SortOrder

register = Library()


@register.inclusion_tag("datatables/datatable.html")
def datatable(table_config: DatatableConfig, table_id: str, class_name: str = ''):
    return {"rich_columns": table_config.enabled_columns, "table_id": table_id, "class_name": class_name}


@register.inclusion_tag("datatables/datatable_definition.html")
def datatable_definition(
        table_config: DatatableConfig,
        table_id: str,
        url: str,
        data: str = None,
        hide_filter_count: bool = False,
        responsive: Union[bool, str] = False):
    if not isinstance(table_config, DatatableConfig):
        raise ValueError(f"Expected DatatableConfig but got '{table_config}'")
    sort_order = None
    default_sort_order_column = table_config.default_sort_order_column
    try:
        index = table_config.column_index(default_sort_order_column)
        sort_direction = default_sort_order_column.default_sort or SortOrder.ASC
        sort_order = [[index, sort_direction.value]]
    except ValueError:
        pass
    return {
        "rich_columns": table_config.enabled_columns,
        "table_id": table_id,
        "url": url, "data": data,
        "hide_filter_count": hide_filter_count,
        "sort_order": sort_order,
        "responsive": responsive}
