from typing import Union

from django.template import Library

from snpdb.views.datatable_view import DatatableConfig, SortOrder

register = Library()


@register.inclusion_tag("datatables/datatable.html")
def datatable(table_config: DatatableConfig, table_id: str, class_name: str = ''):
    """
    Deprecated, use instead
    <table id="{table_id}" data-datatable-url="{url}" class="{class_name}"></table>
    """
    return {
        "rich_columns": table_config.enabled_columns,
        "search_box_enabled": table_config.search_box_enabled,
        "download_csv_button_enabled": table_config.download_csv_button_enabled,
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
    """
    Deprecated, use table tag as seen in datatable tag instead
    """
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
        "download_csv_button_enabled": table_config.download_csv_button_enabled,
        "expand_client_renderer": table_config.expand_client_renderer,
        "table_id": table_id,
        "url": url, "data": data,
        "hide_filter_count": hide_filter_count,
        "sort_order": sort_order,
        "responsive": responsive}
