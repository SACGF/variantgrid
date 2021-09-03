from typing import Union, Optional

from django.template import Library
from snpdb.views.datatable_view import DatatableConfig, SortOrder

register = Library()


@register.inclusion_tag("datatables/datatable_complete.html")
def datatable_complete(
        table_config: DatatableConfig,
        table_id: str,
        class_name: str = '',
        data: Optional[str] = None,
        hide_filter_count: bool = False,
        responsive: Union[bool, str] = False,
        wait_for: Optional[str] = None):
    """
    :param table_config: The DatatableConfig provided in context, defines the column data and URL and other attributes
    :param table_id: the HTML ID that will be assigned to the table
    :param class_name: Any css classes to assign to the table (some classes will be automatic)
    :param data: JavaScript function that will put parameters on the request (this is how we can filter data with custom filters)
    :param hide_filter_count: Set to true if we don't want "Showing 1 of 10 of 200 results"
    :param responsive: If columns should be automatically put in an expand spot if too many columns to fit - note that this
    is incompatible with ajax expand
    :param wait_for: JavaScript promise, don't start rendering data until this is complete
    """

    if not isinstance(table_config, DatatableConfig):
        raise ValueError(f"Expected DatatableConfig but got '{table_config}'")
    sort_order = None
    default_sort_order_column = table_config.default_sort_order_column
    if responsive:
        class_name += ' responsive'
    try:
        index = table_config.column_index(default_sort_order_column)
        sort_direction = default_sort_order_column.default_sort or SortOrder.ASC
        sort_order = [[index, sort_direction.value]]
    except ValueError:
        pass

    return {
        "rich_columns": table_config.enabled_columns,
        "search_box_enabled": table_config.search_box_enabled,
        "expand_client_renderer": table_config.expand_client_renderer,
        "table_id": table_id,
        "class_name": class_name,
        "url": table_config.data_view,
        "data": data,
        "hide_filter_count": hide_filter_count,
        "sort_order": sort_order,
        "responsive": responsive,
        "wait_for": wait_for
    }


@register.inclusion_tag("datatables/datatable.html")
def datatable(table_config: DatatableConfig, table_id: str, class_name: str = ''):
    return {
        "rich_columns": table_config.enabled_columns,
        "search_box_enabled": table_config.search_box_enabled,
        "table_id": table_id,
        "class_name": class_name,
        "expand_client_renderer": table_config.expand_client_renderer
    }


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
        "search_box_enabled": table_config.search_box_enabled,
        "expand_client_renderer": table_config.expand_client_renderer,
        "table_id": table_id,
        "url": url, "data": data,
        "hide_filter_count": hide_filter_count,
        "sort_order": sort_order,
        "responsive": responsive}
