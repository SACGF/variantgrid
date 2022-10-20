from typing import Any, Optional

from django.template.library import Library
from django.urls import reverse

from library.django_utils import get_url_from_view_path
from variantgrid.perm_path import get_visible_url_names

register = Library()
UNSET = '!@#$%^&*()'


@register.inclusion_tag("uicore/tags/menu_item.html", takes_context=True)
def menu_top(context,
             url_name: str,
             app_name: str,
             title: str = None):
    url = None
    for url_name_part in url_name.split('|'):
        if get_visible_url_names().get(url_name_part):
            url = reverse(url_name_part)
            break

    if not url:
        return {'invalid': True}

    is_active = False
    request = context.request
    app_names = app_name.split('|')
    for app_name_part in app_names:
        if request.get_full_path().startswith('/' + app_name_part):
            is_active = True
            break

    if not title:
        title = url_name.replace('_', ' ')
    title = title[0].upper() + title[1::]

    return {
        'url': url,
        'active': is_active,
        'title': title,
        'id': f'menu-top-{app_names[0]}'
    }


@register.inclusion_tag("uicore/tags/menu_item.html", takes_context=True)
def menu_item(
        context,
        url_name: str,
        css_class: str = '',
        arg1: Any = UNSET,
        arg2: Any = UNSET,
        badge_count: Optional[int] = None,
        title: Optional[str] = None,
        icon: Optional[str] = None,
        href: Optional[str] = None,
        admin_only=False,
        external=False,
        other_urls: Optional[str] = None):

    request = context.request
    if admin_only and not request.user.is_superuser:
        return {'invalid': True}

    args = []
    if arg1 != UNSET:
        args.append(arg1)
        if arg2 != UNSET:
            args.append(arg2)

    url = href
    if not url:
        if not get_visible_url_names().get(url_name):
            return {'invalid': True}
        url = reverse(url_name, args=args)

    current_url_name: str = context.request.resolver_match.url_name
    is_active = url_name == current_url_name
    if not is_active and other_urls:
        other_url_parts = other_urls.split(',')
        current_url = context.request.path
        for other_url_part in other_url_parts:
            if current_url.startswith(other_url_part) or current_url_name == other_url_part:
                is_active = True
                break

    if not title:
        title = url_name.replace('_', ' ')
        title_case = ''
        next_capital = True
        for letter in title:
            if next_capital:
                letter = letter.upper()
                next_capital = False
            if letter == ' ':
                next_capital = True
            title_case += letter
        title = title_case

    return {
        'url': url,
        'css_class': css_class,
        'title': title,
        'icon': icon,
        'admin_only': admin_only,
        'active': is_active,
        'badge_count': badge_count,
        'id': f'submenu-{url_name}',
        'external': external
    }


@register.simple_tag()
def absolute_url(name, *args, **kwargs) -> str:
    return get_url_from_view_path(reverse(name, args=args, kwargs=kwargs))


@register.inclusion_tag("uicore/menus/current_record.html", takes_context=True)
def current_record(context, current_record: Any):
    return {"current_record": current_record}
