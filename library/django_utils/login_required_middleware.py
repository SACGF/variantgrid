"""
Site-wide login: Django's LoginRequiredMiddleware, redirecting to settings.LOGIN_URL. A view opts out with
django.contrib.auth.decorators.login_not_required, and a third-party URLconf with login_not_required_include().
DRF marks every APIView and ViewSet login_not_required, so REST endpoints authenticate at view time (DRF's
DEFAULT_PERMISSION_CLASSES) and a token-authenticated request gets DRF's answer rather than a login redirect.
"""
from importlib import import_module

from django.conf import settings
from django.contrib.admin import AdminSite
from django.contrib.auth.decorators import login_not_required
from django.contrib.auth.middleware import LoginRequiredMiddleware
from django.urls import URLPattern, URLResolver, include


class SiteLoginRequiredMiddleware(LoginRequiredMiddleware):
    """ Anonymous requests always go to settings.LOGIN_URL (OIDC on Shariant). Django would follow a view's own
        login_url instead, and the admin, staff_member_required and user_passes_test views name the admin's password
        login form """

    def process_view(self, request, view_func, view_args, view_kwargs):
        # Django marks the admin login login_not_required - keep it behind the site login
        if isinstance(getattr(view_func, "__self__", None), AdminSite) and not request.user.is_authenticated:
            return self.handle_no_permission(request, view_func)
        return super().process_view(request, view_func, view_args, view_kwargs)

    def get_login_url(self, view_func):
        return settings.LOGIN_URL

    def get_redirect_field_name(self, view_func):
        return self.redirect_field_name


def _mark_login_not_required(url_patterns):
    for entry in url_patterns:
        if isinstance(entry, URLResolver):
            _mark_login_not_required(entry.url_patterns)
        elif isinstance(entry, URLPattern):
            login_not_required(entry.callback)


def login_not_required_include(urlconf_module: str):
    """ include() for a third-party URLconf (registration, OIDC) whose views must be reachable before login """
    urlconf = import_module(urlconf_module)
    _mark_login_not_required(urlconf.urlpatterns)
    return include(urlconf)
