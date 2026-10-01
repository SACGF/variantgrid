"""
Renders the menus from settings.MENUS (uicore/menus.py:get_menus): menu_bar_main (the top bar) and menu_bar_sub (the side bar for the
current page), both picked from the request's url name. Also the page chrome's site_messages and the absolute_url tag.
"""
from typing import Optional

from django.http import HttpRequest
from django.template.library import Library
from django.urls import reverse

from library.django_utils import get_url_from_view_path
from uicore.menus import MenuItem, current_menu, get_menus

register = Library()


def _current_url_name(request: HttpRequest) -> Optional[str]:
    if rm := request.resolver_match:
        return rm.url_name
    return None


def _item_context(item: MenuItem, active: bool) -> dict:
    return {
        'url': item.href or reverse(item.url_name),
        'css_class': item.css_class,
        'title': item.display_title,
        'icon': item.icon,
        'admin_only': item.admin_only,
        'active': active,
        'type': 'side',
        'id': f'submenu-{item.url_name}',
        'external': item.external,
        'method': item.method,
    }


@register.inclusion_tag("uicore/menus/menu_bar_main.html", takes_context=True)
def menu_bar_main(context):
    active_menu = current_menu(_current_url_name(context.request))
    top_items = [{
        'url': reverse(menu.url_name),
        'title': menu.title,
        'active': menu == active_menu,
        'type': 'top',
        'id': f'menu-top-{menu.key}',
        'method': 'get',
    } for menu in get_menus() if menu.in_top_bar]
    return {
        'top_items': top_items,
        'help_url': context.get('help_url'),
        'user': context.get('user'),
    }


@register.inclusion_tag("uicore/menus/menu_bar_sub.html", takes_context=True)
def menu_bar_sub(context):
    request = context.request
    url_name = _current_url_name(request)
    if not (menu := current_menu(url_name)):
        return {}

    items = [_item_context(item, item.owns(url_name)) for item in menu.visible_items()
             if request.user.is_superuser or not item.admin_only]
    return {
        'items': items,
        'footer_template': menu.footer_template,
    }


@register.simple_tag()
def absolute_url(name, *args, **kwargs) -> str:
    return get_url_from_view_path(reverse(name, args=args, kwargs=kwargs))


@register.inclusion_tag("uicore/site_messages/site_messages.html")
def site_messages(site_messages):
    return {"site_messages": site_messages}
