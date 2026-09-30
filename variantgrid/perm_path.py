"""

@see https://github.com/SACGF/variantgrid/wiki/URL---Menu-configuration

Replace path() with perm_path in your urls.urlpatterns[]

You can't return None to urlpatterns so need to redirect or something...

A named view goes through path(); a DRF router's patterns, which Django builds itself and
so never pass through here, go through router_urls(). deprecated_path() is path() for a URL we want
to remove but an external client might still call: each hit is reported to Rollbar.

"""
import logging
from collections import defaultdict
from collections.abc import Mapping
from functools import wraps

from django.conf import settings
from django.urls.conf import path as django_path
from django.urls.resolvers import URLPattern, get_resolver

from library.cache import timed_cache
from library.django_utils import require_superuser
from library.django_utils.view_utils import view_to_string
from library.log_utils import report_message


def url_name_enabled(name: str) -> bool:
    """ URLS_NAME_REGISTER, and for the names in URLS_NAME_REGISTER_SETTINGS the feature setting as it stands now """
    if not settings.URLS_NAME_REGISTER[name]:
        return False
    if setting_name := settings.URLS_NAME_REGISTER_SETTINGS.get(name):
        return bool(getattr(settings, setting_name))
    return True


def _perm_path(route, view, path_func, **kwargs):
    name = kwargs.get('name')
    if name is not None:
        if not url_name_enabled(name):
            view = require_superuser(view)
    else:
        logging.warning("url: route='%s' view=%s, has no name, so is not tested via URLS_NAME_REGISTER",
                        route, view_to_string(view))
    return path_func(route, view, **kwargs)


def path(route, view, **kwargs):
    return _perm_path(route, view, django_path, **kwargs)


def deprecated_path(route, view, **kwargs):
    """ A URL with no callers we know of, kept in case an external client uses it (#1475).
        Remove once Rollbar has gone ~6 months without reporting it """
    name = kwargs.get('name')

    @wraps(view)
    def report_deprecated_view(request, *args, **view_kwargs):
        report_message(f"Deprecated URL '{name}' accessed", level='warning', request=request,
                       extra_data={"target": request.get_full_path()})
        return view(request, *args, **view_kwargs)

    return path(route, report_deprecated_view, **kwargs)


def router_urls(router) -> list[URLPattern]:
    """ A DRF router's patterns with URLS_NAME_REGISTER applied, as path() does for a named view """
    urls = []
    for url in router.urls:
        if url.name is not None and not url_name_enabled(url.name):
            url = URLPattern(url.pattern, require_superuser(url.callback), url.default_args, url.name)
        urls.append(url)
    return urls


@timed_cache()
def get_visible_url_names() -> Mapping[str, bool]:
    # Only include loaded URLs, then use what's configured in URLS_NAME_REGISTER
    url_name_visible = defaultdict(lambda: False)
    for url_name in (k for k in get_resolver(None).reverse_dict.keys() if isinstance(k, str)):
        url_name_visible[url_name] = url_name_enabled(url_name)
    return url_name_visible
