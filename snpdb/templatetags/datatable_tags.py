from django.template import Library

from snpdb.models import UserGridConfig

register = Library()


@register.inclusion_tag("datatables/user_data_grid_filter.html", takes_context=True)
def datatable_user_data_filter(context, table_id: str, grid_name: str,
                               group_data: bool = True, hidden_data: bool = False):
    """ Toggles stored per user in UserGridConfig, see filter_queryset in the table config """
    user_grid_config = UserGridConfig.get(context["user"], grid_name)
    return {
        "table_id": table_id,
        "grid_name": grid_name,
        "group_data": group_data,
        "hidden_data": hidden_data,
        "user_grid_config": user_grid_config,
    }
