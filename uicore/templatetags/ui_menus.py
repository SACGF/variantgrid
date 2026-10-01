"""
absolute_url: a full URL (scheme and host) for a url name, for emails. The menus themselves are rendered by
uicore/chrome.py and fetched by the page after it loads.
"""
from django.template.library import Library
from django.urls import reverse

from library.django_utils import get_url_from_view_path

register = Library()


@register.simple_tag()
def absolute_url(name, *args, **kwargs) -> str:
    return get_url_from_view_path(reverse(name, args=args, kwargs=kwargs))
