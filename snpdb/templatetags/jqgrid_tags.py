import uuid

from django.template import Library, loader
from django.urls.base import reverse
from django.urls.exceptions import NoReverseMatch

from snpdb.models import UserSettings, UserGridConfig

register = Library()


def get_grid_id(grid_name):
    return f"{grid_name}-grid"


@register.simple_tag(takes_context=True)
def jqgrid(context, jqgrid_config_url, name=None, search=True, gbox_id=None,
           template_name='jqgrid/jqgrid.html', init_func=None, jqgrid_config_get_parameters_func=None,
           delete=False, download_grid_json_as_csv=False, *args, **kwargs):
    grid_complete = kwargs.pop("grid_complete", None)
    modify_export_url = kwargs.pop("modify_export_url", None)

    grid_id = 'grid'
    pager_id = 'pager'
    if name is not None:
        name = str(name).replace(" ", "_")
        grid_id = get_grid_id(name)
        pager_id = name + '-' + pager_id

    if gbox_id is None:
        gbox_id = 'gbox_' + grid_id

    unique_code = '_%s' % uuid.uuid4().hex.replace('-', '')

    csrf_token = context["csrf_token"]
    try:
        # New way - passed general JQGrid url and set op=config
        url_kwargs = {"op": "config"}
        url_kwargs.update(kwargs)
        jqgrid_config_url = reverse(jqgrid_config_url, args=args, kwargs=url_kwargs)
    except NoReverseMatch:
        # Old way - passed special config url
        jqgrid_config_url = reverse(jqgrid_config_url, args=args, kwargs=kwargs)

    context = {
        'jqgrid_config_url': jqgrid_config_url,
        'grid_id': grid_id,
        'unique_code': unique_code,
        'pager_id': pager_id,
        'gbox_id': gbox_id,
        'search': search,
        'download_grid_json_as_csv': download_grid_json_as_csv,
        'init_func': init_func,
        'init_kwargs': kwargs,
        'jqgrid_config_get_parameters_func': jqgrid_config_get_parameters_func,
        "grid_complete": grid_complete,
        "modify_export_url": modify_export_url,
        'delete': delete,
        'csrf_token': csrf_token,
        'user': context["user"],
        'user_settings': UserSettings.get_for_user(context["user"])
    }
    t = loader.get_template(template_name)
    return t.render(context=context)


@register.inclusion_tag('jqgrid/user_data_grid_filter.html', takes_context=True)
def user_data_grid_filter(context, grid_id_suffix, caption, group_data=True, incomplete_data=False, hidden_data=False, filter_name_choices=None):
    user = context["user"]
    user_grid_config = UserGridConfig.get(user, caption)

    return {
        'group_data': group_data,
        'incomplete_data': incomplete_data,
        'hidden_data': hidden_data,
        'filter_name_choices': filter_name_choices,
        'caption': caption,
        'grid_id': get_grid_id(grid_id_suffix),
        'user_grid_config': user_grid_config
    }
