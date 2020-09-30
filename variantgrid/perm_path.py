"""

@see https://github.com/SACGF/variantgrid/wiki/URL---Menu-configuration

Replace path() with perm_path in your urls.urlpatterns[]

You can't return None to urlpatterns so need to redirect or something...

"""
import logging
from collections import defaultdict
from typing import Mapping

from django.conf import settings
from django.urls.conf import path, re_path
from django.urls.resolvers import get_resolver

from library.cache import timed_cache
from library.django_utils import require_superuser


def _perm_path(route, view, path_func, **kwargs):
    name = kwargs.get('name')
    if name is not None:
        if not settings.URLS_NAME_REGISTER[name]:
            view = require_superuser(view)
    else:
        logging.warning("url: '%s' has no name, so is not tested via URLS_NAME_REGISTER", route)
    return path_func(route, view, **kwargs)


def perm_path(route, view, **kwargs):
    return _perm_path(route, view, path, **kwargs)


def re_perm_path(route, view, **kwargs):
    return _perm_path(route, view, re_path, **kwargs)


@timed_cache()
def get_visible_url_names() -> Mapping[str, bool]:
    # Only include loaded URLs, then use what's configured in URLS_NAME_REGISTER
    url_name_visible = defaultdict(lambda: False)
    for url_name in (k for k in get_resolver(None).reverse_dict.keys() if isinstance(k, str)):
        url_name_visible[url_name] = settings.URLS_NAME_REGISTER[url_name]
    return url_name_visible
