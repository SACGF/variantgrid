from typing import Any

from django.template.library import Library
from django.urls import reverse

from library.django_utils import get_url_from_view_path
from variantgrid.perm_path import get_visible_url_names

register = Library()
UNSET = '!@#$%^&*()'

@register.inclusion_tag("snpdb/tags/menu_item.html", takes_context=True)
def menu_top(
        context,
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

    css_classes = []
    request = context.request
    for app_name_part in app_name.split('|'):
        if app_name_part in request.get_full_path():
            css_classes.append('current')
            break

    css_class = ' '.join(css_classes).strip()

    if not title:
        title = url_name.replace('_', ' ')

    return {
        'url': url,
        'css_class': css_class,
        'title': title
    }

@register.inclusion_tag("snpdb/tags/menu_item.html", takes_context=True)
def menu_item(
        context,
        url_name:str,
        css_class: str = '',
        arg1: Any = UNSET,
        arg2: Any = UNSET,
        badge_count: int = None,
        title: str = None,
        href: str = None,
        admin_only=False):

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

    current_url_name = context.request.resolver_match.url_name
    is_current = url_name == current_url_name

    css_classes = [css_class]
    if is_current:
        css_classes.append('current')
    if admin_only:
        css_classes.append('admin')
    css_class = ' '.join(css_classes).strip()

    if not title:
        title = url_name.replace('_', ' ')
        title_case = ''
        next_capital = True
        for index in range(len(title)):
            letter = title[index]
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
        'badge_count': badge_count
    }


@register.simple_tag()
def absolute_url(name, *args, **kwargs) -> str:
    return get_url_from_view_path(reverse(name, args=args, kwargs=kwargs))
