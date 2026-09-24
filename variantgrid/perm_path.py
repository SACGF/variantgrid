"""

@see https://github.com/SACGF/variantgrid/wiki/URL---Menu-configuration

Replace path() with perm_path in your urls.urlpatterns[]

You can't return None to urlpatterns so need to redirect or something...

A named view goes through path(); a DRF router's patterns, which Django builds itself and
so never pass through here, go through router_urls().

"""
import logging
from collections import defaultdict
from collections.abc import Mapping

from django.conf import settings
from django.urls.conf import path as django_path
from django.urls.resolvers import URLPattern, get_resolver

from library.cache import timed_cache
from library.django_utils import require_superuser
from library.utils.django_utils import view_to_string


def _perm_path(route, view, path_func, **kwargs):
    name = kwargs.get('name')
    if name is not None:
        if not settings.URLS_NAME_REGISTER[name]:
            view = require_superuser(view)
    else:
        logging.warning("url: route='%s' view=%s, has no name, so is not tested via URLS_NAME_REGISTER",
                        route, view_to_string(view))
    return path_func(route, view, **kwargs)


def path(route, view, **kwargs):
    return _perm_path(route, view, django_path, **kwargs)


def router_urls(router) -> list[URLPattern]:
    """ A DRF router's patterns with URLS_NAME_REGISTER applied, as path() does for a named view """
    urls = []
    for url in router.urls:
        if url.name is not None and not settings.URLS_NAME_REGISTER[url.name]:
            url = URLPattern(url.pattern, require_superuser(url.callback), url.default_args, url.name)
        urls.append(url)
    return urls


@timed_cache()
def get_visible_url_names() -> Mapping[str, bool]:
    # Only include loaded URLs, then use what's configured in URLS_NAME_REGISTER
    url_name_visible = defaultdict(lambda: False)
    for url_name in (k for k in get_resolver(None).reverse_dict.keys() if isinstance(k, str)):
        url_name_visible[url_name] = settings.URLS_NAME_REGISTER[url_name]
    return url_name_visible
